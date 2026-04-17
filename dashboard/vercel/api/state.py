"""
Vercel serverless function — aggregates Orb Cloud state into JSON.

Secrets (set via `vercel env add`):
  ORB_API_KEY        bearer token for api.orbcloud.dev
  ORB_COMPUTER_ID    UUID of the computer to observe

No per-request secret escapes to the client; the public frontend only
ever sees the aggregated state object.

Route: GET /api/state   (via vercel.json)
"""

from __future__ import annotations

import json
import os
import re
import time
import urllib.parse
import urllib.request
from http.server import BaseHTTPRequestHandler

ORB_BASE = os.environ.get("ORB_BASE_URL", "https://api.orbcloud.dev").rstrip("/")
ORB_API_KEY = os.environ.get("ORB_API_KEY", "").strip()
COMPUTER_ID = os.environ.get("ORB_COMPUTER_ID", "").strip()

# Lightweight in-function TTL cache so repeated polls don't hammer Orb.
_CACHE: dict = {"t": 0.0, "data": None}
_CACHE_TTL = 4.0  # seconds


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


def _file_text(path: str) -> str:
    try:
        return _orb_get(f"/v1/computers/{COMPUTER_ID}/files/{path.lstrip('/')}", text=True)
    except Exception:
        return ""


def _file_dir(path: str) -> list[dict]:
    try:
        data = _orb_get(f"/v1/computers/{COMPUTER_ID}/files/{path.lstrip('/')}")
        return data.get("files", []) if isinstance(data, dict) else []
    except Exception:
        return []


def _file_json(path: str):
    raw = _file_text(path)
    if not raw:
        return None
    try:
        return json.loads(raw)
    except Exception:
        return None


def _file_jsonl(path: str) -> list:
    raw = _file_text(path)
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


LOG_RUN_RE = re.compile(
    r"\[(?P<ts>[\d\-: ]+)\] \[agent-sdk\] Run complete\. "
    r"Session: (?P<sid>[\w\-]+), Error: (?P<err>True|False), "
    r"Turns: (?P<turns>\d+), Cost: \$(?P<cost>[\d.]+)"
)
LOG_RUN_START_RE = re.compile(
    r"\[(?P<ts>[\d\-: ]+)\] \[agent-sdk\] === RUN #(?P<n>\d+) "
)


def _parse_runs(log_text: str, max_runs: int = 80) -> list[dict]:
    if not log_text:
        return []
    lines = log_text.splitlines()
    runs: list[dict] = []
    current = None
    for line in lines:
        m0 = LOG_RUN_START_RE.match(line)
        if m0:
            current = int(m0.group("n"))
            continue
        m1 = LOG_RUN_RE.match(line)
        if m1:
            runs.append({
                "run_num": current if current is not None else len(runs) + 1,
                "timestamp": m1.group("ts"),
                "session_id": m1.group("sid"),
                "error": m1.group("err") == "True",
                "turns": int(m1.group("turns")),
                "cost_usd": float(m1.group("cost")),
                "error_message": None,
            })
    return runs[-max_runs:]


def _load_candidate_list() -> list[dict]:
    """Return a chronological list of candidates by scanning
    findings/candidates/*/candidate.json."""
    dirs = [
        f for f in _file_dir("agent/code/findings/candidates")
        if f.get("type") == "directory" and f.get("name") not in (".", "..")
    ]
    out: list[dict] = []
    # Cap at 15 most recently named (candidate id encodes the date)
    for d in sorted(dirs, key=lambda x: x.get("name", ""), reverse=True)[:15]:
        meta = _file_json(f"agent/code/findings/candidates/{d['name']}/candidate.json")
        if meta:
            out.append({
                "candidate_id": meta.get("candidate_id") or d["name"],
                "name": meta.get("name"),
                "designed_at": meta.get("designed_at"),
                "smiles": meta.get("smiles"),
            })
        else:
            out.append({"candidate_id": d["name"], "name": None, "designed_at": None})
    return out


def _aggregate() -> dict:
    ts = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    out: dict = {"timestamp": ts, "error": None}

    if not (ORB_API_KEY and COMPUTER_ID):
        out["error"] = "ORB_API_KEY or ORB_COMPUTER_ID not set on Vercel"
        return out

    try:
        out["computer"] = _orb_get(f"/v1/computers/{COMPUTER_ID}")
    except Exception as exc:
        out["error"] = f"computer fetch failed: {exc}"
        out["computer"] = None

    try:
        out["agents"] = (_orb_get(f"/v1/computers/{COMPUTER_ID}/agents") or {}).get("agents", [])
    except Exception:
        out["agents"] = []

    log_text = _file_text("agent/code/logs/agent-sdk.log")
    out["recent_log"] = "\n".join(log_text.splitlines()[-300:])
    out["runs"] = _parse_runs(log_text)

    snap = _file_json("agent/code/logs/env-snapshot.json") or {}
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
    out["env_snapshot"] = env_out

    # findings counts
    findings = {}
    for name in ("candidates", "docking", "admet", "novelty", "mechanism",
                 "red-team", "retrosynthesis", "weekly-reports"):
        items = _file_dir(f"agent/code/findings/{name}")
        if name == "candidates":
            n = sum(1 for f in items if f.get("type") == "directory"
                    and f.get("name") not in (".gitkeep", ".", ".."))
        else:
            n = sum(1 for f in items if f.get("type") == "file"
                    and f.get("name") != ".gitkeep")
        findings[name.replace("-", "_")] = n
    out["findings"] = findings

    out["candidates_list"] = _load_candidate_list()

    out["leaderboard"] = _file_json("agent/code/findings/leaderboard.json") or {
        "candidates": [], "updated_at": None,
    }

    pc = _file_jsonl("agent/code/findings/loop-health/positive-control.jsonl")
    nc = _file_jsonl("agent/code/findings/loop-health/negative-control.jsonl")
    pa = _file_json("agent/code/findings/loop-health/positive-control-alert.json")
    na = _file_json("agent/code/findings/loop-health/negative-control-alert.json")
    out["loop_health"] = {
        "positive_control": pc, "negative_control": nc,
        "positive_alert": pa, "negative_alert": na,
    }

    return out


def _get_cached():
    now = time.time()
    if _CACHE["data"] is not None and (now - _CACHE["t"]) < _CACHE_TTL:
        return _CACHE["data"]
    try:
        data = _aggregate()
    except Exception as exc:
        data = {"error": f"aggregate crashed: {exc}", "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    _CACHE["t"] = now
    _CACHE["data"] = data
    return data


class handler(BaseHTTPRequestHandler):  # noqa: N801 — Vercel Python runtime looks for class handler
    def do_GET(self):  # noqa: N802
        body = json.dumps(_get_cached()).encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Cache-Control", "no-store, max-age=0")
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)
