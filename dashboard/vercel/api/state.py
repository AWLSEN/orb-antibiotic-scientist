"""
Vercel serverless function — aggregates Orb Cloud state into JSON.

Secrets (set via `vercel env add`):
  ORB_API_KEY        bearer token for api.orbcloud.dev

No per-request secret escapes to the client; the public frontend only
ever sees the aggregated state object.

Response shape:
  {
    "timestamp": "<iso>",
    "usage": { runtime_gb_hours, disk_gb_hours, checkpoint_cycles, ... },
    "computers": [
      { id, name, status, runtime_mb, disk_mb, agent_state, agent_pid,
        findings_candidates, session_started_at, last_event_at }
    ],
    "detail": null | { /* full detail for ?computer=<id> */ }
  }

Route: GET /api/state               → overview (usage + computers summary)
       GET /api/state?computer=ID   → overview + full detail for that id
"""

from __future__ import annotations

import json
import os
import re
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from http.server import BaseHTTPRequestHandler

ORB_BASE = os.environ.get("ORB_BASE_URL", "https://api.orbcloud.dev").rstrip("/")
ORB_API_KEY = os.environ.get("ORB_API_KEY", "").strip()

# Short in-process cache so repeated polls from many clients don't hammer Orb.
# Separate cache keys for overview vs detail-per-id.
_CACHE: dict[str, dict] = {}
_CACHE_TTL = 30.0   # bumped 4s → 30s; agent produces ~1 cand/min, staleness is fine

# Orb file requests issued in parallel.
_MAX_PARALLEL = 16


def _orb_get(path: str, *, text: bool = False, timeout: int = 8):
    if not ORB_API_KEY:
        raise RuntimeError("ORB_API_KEY not configured on this deployment")
    url = f"{ORB_BASE}{path}"
    req = urllib.request.Request(
        url, headers={"Authorization": f"Bearer {ORB_API_KEY}"}
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        raw = resp.read()
    return raw.decode("utf-8", errors="replace") if text else json.loads(raw)


def _file_text(cid: str, path: str) -> str:
    try:
        return _orb_get(
            f"/v1/computers/{cid}/files/{path.lstrip('/')}", text=True
        )
    except Exception:
        return ""


def _file_dir(cid: str, path: str) -> list[dict]:
    try:
        data = _orb_get(f"/v1/computers/{cid}/files/{path.lstrip('/')}")
        return data.get("files", []) if isinstance(data, dict) else []
    except Exception:
        return []


def _file_json(cid: str, path: str):
    raw = _file_text(cid, path)
    if not raw:
        return None
    try:
        return json.loads(raw)
    except Exception:
        return None


def _file_jsonl(cid: str, path: str) -> list:
    raw = _file_text(cid, path)
    if not raw:
        return []
    out = []
    for line in raw.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            out.append(json.loads(line))
        except Exception:
            continue
    return out


LOG_RUN_START_RE = re.compile(
    r"^\[(?P<ts>\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})\] \[agent-sdk\] === RUN #(?P<n>\d+)"
)
LOG_RUN_END_RE = re.compile(
    r"^\[(?P<ts>\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})\] \[agent-sdk\] Run complete\. "
    r"Session: (?P<sid>[\w\-]+), Error: (?P<err>True|False), Turns: (?P<turns>\d+), "
    r"Cost: \$(?P<cost>[\d.]+)"
)
LOG_RESTART_RE = re.compile(
    r"^\[(?P<ts>\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})\] \[agent-sdk\] Restarting in "
)
LOG_STARTUP_RE = re.compile(
    r"^\[(?P<ts>\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})\] \[agent-sdk\] "
    r"orb-antibiotic-scientist -- Claude Agent SDK"
)
LOG_TS_RE = re.compile(r"^\[(?P<ts>\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})\]")


def _ts_to_iso(ts: str) -> str:
    return ts.replace(" ", "T") + "Z"


def _parse_log_events(log_text: str) -> list[dict]:
    """Extract structured events from agent-sdk.log."""
    events: list[dict] = []
    for line in log_text.splitlines():
        m = LOG_RUN_START_RE.match(line)
        if m:
            events.append({
                "kind": "run_start",
                "ts": _ts_to_iso(m.group("ts")),
                "run_num": int(m.group("n")),
            })
            continue
        m = LOG_RUN_END_RE.match(line)
        if m:
            events.append({
                "kind": "run_end",
                "ts": _ts_to_iso(m.group("ts")),
                "session_id": m.group("sid"),
                "error": m.group("err") == "True",
                "turns": int(m.group("turns")),
                "cost_usd": float(m.group("cost")),
            })
            continue
        m = LOG_RESTART_RE.match(line)
        if m:
            events.append({
                "kind": "restart",
                "ts": _ts_to_iso(m.group("ts")),
            })
            continue
        m = LOG_STARTUP_RE.match(line)
        if m:
            events.append({
                "kind": "startup",
                "ts": _ts_to_iso(m.group("ts")),
            })
    return events


def _session_started_at(events: list[dict]) -> str | None:
    for e in events:
        if e["kind"] in ("startup", "run_start"):
            return e["ts"]
    return None


def _last_event_ts(log_text: str) -> str | None:
    for line in reversed(log_text.splitlines()):
        m = LOG_TS_RE.match(line)
        if m:
            return _ts_to_iso(m.group("ts"))
    return None


def _list_computers() -> list[dict]:
    try:
        data = _orb_get("/v1/computers")
        return data.get("computers", []) if isinstance(data, dict) else []
    except Exception:
        return []


def _computer_agent_state(cid: str) -> dict:
    try:
        data = _orb_get(f"/v1/computers/{cid}/agents")
        a = (data.get("agents") or [{}])[0] if data.get("agents") else {}
        return {
            "state": a.get("state"),
            "pid": a.get("pid"),
            "port": a.get("port"),
        }
    except Exception:
        return {"state": None, "pid": None, "port": None}


def _candidate_count(cid: str) -> int:
    items = _file_dir(cid, "agent/code/findings/candidates")
    return sum(
        1 for f in items
        if f.get("type") == "directory" and f.get("name") not in (".", "..")
    )


def _summary_for_computer(c: dict) -> dict:
    cid = c.get("id")
    agent = _computer_agent_state(cid)
    log_text = _file_text(cid, "agent/code/logs/agent-sdk.log")
    events = _parse_log_events(log_text)
    return {
        "id": cid,
        "name": c.get("name"),
        "status": c.get("status"),
        "runtime_mb": c.get("runtime_mb"),
        "disk_mb": c.get("disk_mb"),
        "agent_state": agent["state"],
        "agent_pid": agent["pid"],
        "agent_port": agent["port"],
        "findings_candidates": _candidate_count(cid),
        "session_started_at": _session_started_at(events),
        "last_event_at": _last_event_ts(log_text),
        "events_count": len(events),
    }


def _load_one_candidate(cid: str, name: str) -> dict:
    """Fetch the files needed to derive a candidate's experiment status.

    Issues 2–5 Orb file requests per candidate. Called in parallel via
    a thread pool to keep total latency bounded by the slowest request
    rather than the sum.
    """
    base = f"agent/code/findings/candidates/{name}"
    inner = _file_dir(cid, base)
    artifacts = [f.get("name") for f in inner if f.get("type") == "file"]
    has = set(artifacts)

    meta = _file_json(cid, f"{base}/candidate.json") or {}
    scored = _file_json(cid, f"{base}/scored.json") if "scored.json" in has else None
    docking = _file_json(cid, f"{base}/docking.json") if "docking.json" in has else None
    redteam = _file_json(cid, f"{base}/redteam.json") if "redteam.json" in has else None

    status = "in_progress"
    reason = None
    rigor = None
    dg = None
    passed_gates = 1
    total_gates = 12

    if "validate.json" in has: passed_gates = max(passed_gates, 4)
    if "docking.json" in has: passed_gates = max(passed_gates, 5)
    if "docking-secondary.json" in has: passed_gates = max(passed_gates, 6)
    if "mechanism.json" in has: passed_gates = max(passed_gates, 7)
    if "admet.json" in has: passed_gates = max(passed_gates, 8)
    if "novelty.json" in has: passed_gates = max(passed_gates, 9)
    if "retrosynthesis.json" in has: passed_gates = max(passed_gates, 11)
    if "redteam.json" in has: passed_gates = max(passed_gates, 12)

    if docking and isinstance(docking, dict):
        dg = docking.get("best_energy_kcalmol")
        thr = docking.get("threshold_kcalmol") or -8.0
        if dg is not None and dg > thr:
            status = "failed"
            reason = f"docking too weak — {dg:.2f} kcal/mol (needed ≤ {thr:.1f})"

    if redteam and isinstance(redteam, dict):
        if redteam.get("substantive_flaw"):
            status = "failed"
            flaws = redteam.get("flaws") or []
            sub = next((f for f in flaws if f.get("substantive")), None)
            reason = sub.get("description") if sub else "red-team vetoed"
            if reason and len(reason) > 140:
                reason = reason[:138] + "…"
        else:
            if scored and scored.get("above_threshold"):
                status = "passed"
                rigor = scored.get("rigor_score")

    if scored and isinstance(scored, dict):
        rigor = scored.get("rigor_score")
        if status == "in_progress" and scored.get("above_threshold") is True:
            status = "passed"
        elif status == "in_progress" and scored.get("veto_applied"):
            status = "failed"
            reason = scored.get("veto_reason") or "vetoed by pipeline"

    return {
        "candidate_id": meta.get("candidate_id") or name,
        "name": meta.get("name"),
        "designed_at": meta.get("designed_at"),
        "smiles": meta.get("smiles"),
        "scaffold_class": (meta.get("design_rationale") or {}).get("scaffold_class"),
        "status": status,
        "reason": reason,
        "rigor_score": rigor,
        "docking_dg": dg,
        "passed_gates": passed_gates,
        "total_gates": total_gates,
    }


def _load_candidate_list(cid: str, limit: int = 18) -> list[dict]:
    """Parallelised candidate fetch. Issues up to _MAX_PARALLEL simultaneous
    Orb file requests so total latency is bounded by the slowest request."""
    dirs = [
        f for f in _file_dir(cid, "agent/code/findings/candidates")
        if f.get("type") == "directory" and f.get("name") not in (".", "..", ".gitkeep")
    ]
    names = [d["name"] for d in sorted(dirs, key=lambda x: x.get("name", ""), reverse=True)[:limit]]

    if not names:
        return []

    results: list[dict | None] = [None] * len(names)
    with ThreadPoolExecutor(max_workers=_MAX_PARALLEL) as pool:
        futures = {pool.submit(_load_one_candidate, cid, n): i for i, n in enumerate(names)}
        for fut in as_completed(futures):
            i = futures[fut]
            try:
                results[i] = fut.result()
            except Exception:
                results[i] = {"candidate_id": names[i], "status": "in_progress"}
    return [r for r in results if r is not None]


def _detail_for_computer(cid: str) -> dict:
    log_text = _file_text(cid, "agent/code/logs/agent-sdk.log")
    events = _parse_log_events(log_text)
    snap = _file_json(cid, "agent/code/logs/env-snapshot.json") or {}
    env_out = {}
    for k in (
        "ANTHROPIC_BASE_URL", "ANTHROPIC_AUTH_TOKEN", "ANTHROPIC_API_KEY",
        "ORB_PROXY_URL", "API_TIMEOUT_MS",
        "ANTHROPIC_DEFAULT_OPUS_MODEL",
        "ANTHROPIC_DEFAULT_SONNET_MODEL",
        "ANTHROPIC_DEFAULT_HAIKU_MODEL",
    ):
        v = snap.get(k)
        if v is None or v == "<MISSING>":
            env_out[k] = "<MISSING>"
        elif k in ("ANTHROPIC_AUTH_TOKEN", "ANTHROPIC_API_KEY"):
            env_out[k] = "<set>" if v not in ("", "<redacted>") else "<MISSING>"
        else:
            env_out[k] = str(v)

    findings = {}
    for name in ("candidates", "docking", "admet", "novelty", "mechanism",
                 "red-team", "retrosynthesis", "weekly-reports"):
        items = _file_dir(cid, f"agent/code/findings/{name}")
        if name == "candidates":
            n = sum(1 for f in items if f.get("type") == "directory"
                    and f.get("name") not in (".gitkeep", ".", ".."))
        else:
            n = sum(1 for f in items if f.get("type") == "file"
                    and f.get("name") != ".gitkeep")
        findings[name.replace("-", "_")] = n

    return {
        "computer_id": cid,
        "events": events,
        "env_snapshot": env_out,
        "findings": findings,
        "candidates_list": _load_candidate_list(cid, limit=30),
        "leaderboard": _file_json(cid, "agent/code/findings/leaderboard.json") or {
            "candidates": [], "updated_at": None,
        },
        "loop_health": {
            "positive_control": _file_jsonl(cid, "agent/code/findings/loop-health/positive-control.jsonl"),
            "negative_control": _file_jsonl(cid, "agent/code/findings/loop-health/negative-control.jsonl"),
            "positive_alert": _file_json(cid, "agent/code/findings/loop-health/positive-control-alert.json"),
            "negative_alert": _file_json(cid, "agent/code/findings/loop-health/negative-control-alert.json"),
        },
        "recent_log": "\n".join(log_text.splitlines()[-250:]),
    }


def _usage() -> dict:
    try:
        u = _orb_get("/v1/usage")
        # Rough $ estimate based on typical cloud-compute GB-hour rates.
        # Orb pricing is not public at this level of granularity, so label
        # clearly as an estimate. $0.03/GB-hr runtime, $0.001/GB-hr disk.
        rt_usd = float(u.get("runtime_gb_hours", 0)) * 0.03
        disk_usd = float(u.get("disk_gb_hours", 0)) * 0.001
        u["estimated_cost_usd"] = round(rt_usd + disk_usd, 4)
        u["estimate_note"] = "rough estimate at $0.03/GB-hr runtime + $0.001/GB-hr disk"
        return u
    except Exception as exc:
        return {"error": str(exc)}


def _aggregate(computer_id: str | None) -> dict:
    ts = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    out: dict = {"timestamp": ts, "error": None}

    if not ORB_API_KEY:
        out["error"] = "ORB_API_KEY not set on Vercel"
        return out

    out["usage"] = _usage()

    computers = _list_computers()
    out["computers"] = [_summary_for_computer(c) for c in computers]

    if computer_id:
        # validate it exists in this org
        if any(c.get("id") == computer_id for c in computers):
            out["detail"] = _detail_for_computer(computer_id)
        else:
            out["detail"] = {"error": f"unknown computer id {computer_id}"}
    else:
        out["detail"] = None

    return out


def _get_cached(cache_key: str, computer_id: str | None) -> dict:
    now = time.time()
    cached = _CACHE.get(cache_key)
    if cached and (now - cached["t"]) < _CACHE_TTL:
        return cached["data"]
    try:
        data = _aggregate(computer_id)
    except Exception as exc:
        data = {"error": f"aggregate crashed: {exc}",
                "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    _CACHE[cache_key] = {"t": now, "data": data}
    return data


class handler(BaseHTTPRequestHandler):  # noqa: N801 — Vercel entrypoint
    def do_GET(self):  # noqa: N802
        parsed = urllib.parse.urlparse(self.path)
        qs = urllib.parse.parse_qs(parsed.query)
        computer_id = (qs.get("computer") or [None])[0]
        cache_key = f"computer={computer_id or '_'}"
        body = json.dumps(_get_cached(cache_key, computer_id)).encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Cache-Control", "no-store, max-age=0")
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)
