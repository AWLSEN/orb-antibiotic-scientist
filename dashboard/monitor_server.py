#!/usr/bin/env python3
"""
Local HTTP proxy for dashboard/monitor.html.

Serves:
  GET /                  → dashboard/monitor.html
  GET /api/state         → aggregated JSON snapshot of the Orb computer

The aggregator pulls, via the Orb Cloud REST API:
  - computer status + agent state
  - recent agent-sdk.log (tail-parsed into per-run records)
  - logs/env-snapshot.json (written by agent.py)
  - findings/leaderboard.json
  - findings/loop-health/positive-control.jsonl + alerts
  - findings/loop-health/negative-control.jsonl + alerts
  - findings/candidates, docking, admet, ... counts (directory listings)

Required env vars:
  ORB_API_KEY       Bearer token (reads ~/.orb/credentials as a fallback).
  ORB_COMPUTER_ID   Computer UUID (reads .orb-state/computer-id as fallback).

Optional:
  PORT              default 7777
  HOST              default 127.0.0.1

Usage:
  python dashboard/monitor_server.py
  open http://127.0.0.1:7777/
"""

from __future__ import annotations

import json
import os
import re
import sys
import time
import traceback
import urllib.parse
import urllib.request
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

REPO = Path(__file__).parent.parent
DASH = Path(__file__).parent
HTML = DASH / "monitor.html"

ORB_BASE = os.environ.get("ORB_BASE_URL", "https://api.orbcloud.dev").rstrip("/")
DEFAULT_PORT = int(os.environ.get("PORT", "7777"))
DEFAULT_HOST = os.environ.get("HOST", "127.0.0.1")


def _load_orb_api_key() -> str | None:
    env = os.environ.get("ORB_API_KEY")
    if env:
        return env.strip()
    creds = Path.home() / ".orb" / "credentials"
    if creds.exists():
        try:
            data = json.loads(creds.read_text())
            return data.get("token")
        except Exception:
            return None
    return None


def _load_computer_id() -> str | None:
    env = os.environ.get("ORB_COMPUTER_ID")
    if env:
        return env.strip()
    state = REPO / ".orb-state" / "computer-id"
    if state.exists():
        return state.read_text().strip()
    return None


ORB_API_KEY = _load_orb_api_key()
COMPUTER_ID = _load_computer_id()


# ----------------------------------------------------------------------
# Orb API helpers (synchronous — called at /api/state request time)
# ----------------------------------------------------------------------


def _orb_get(path: str, *, text: bool = False, timeout: int = 8):
    if not ORB_API_KEY:
        raise RuntimeError("ORB_API_KEY not set and no ~/.orb/credentials found")
    url = f"{ORB_BASE}{path}"
    req = urllib.request.Request(
        url, headers={"Authorization": f"Bearer {ORB_API_KEY}"}
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        raw = resp.read()
    if text:
        return raw.decode("utf-8", errors="replace")
    return json.loads(raw)


def _orb_file_text(container_path: str) -> str:
    """Read a file from inside the Orb computer via the files API."""
    path = container_path.lstrip("/")
    try:
        return _orb_get(f"/v1/computers/{COMPUTER_ID}/files/{path}", text=True)
    except Exception:
        return ""


def _orb_dir(container_path: str) -> list[dict]:
    path = container_path.lstrip("/")
    try:
        data = _orb_get(f"/v1/computers/{COMPUTER_ID}/files/{path}")
        return data.get("files", []) if isinstance(data, dict) else []
    except Exception:
        return []


def _orb_file_json(container_path: str):
    raw = _orb_file_text(container_path)
    if not raw:
        return None
    try:
        return json.loads(raw)
    except Exception:
        return None


def _orb_file_jsonl(container_path: str) -> list:
    raw = _orb_file_text(container_path)
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


# ----------------------------------------------------------------------
# Log parsing — turn agent-sdk.log into run records
# ----------------------------------------------------------------------


LOG_RUN_RE = re.compile(
    r"\[(?P<ts>[\d\-: ]+)\] \[agent-sdk\] Run complete\. "
    r"Session: (?P<sid>[\w\-]+), Error: (?P<err>True|False), "
    r"Turns: (?P<turns>\d+), Cost: \$(?P<cost>[\d.]+)"
)
LOG_RUN_START_RE = re.compile(
    r"\[(?P<ts>[\d\-: ]+)\] \[agent-sdk\] === RUN #(?P<n>\d+) "
)


def parse_runs(log_text: str, max_runs: int = 80) -> list[dict]:
    if not log_text:
        return []
    runs: list[dict] = []
    lines = log_text.splitlines()
    current_run_num = None
    for line in lines:
        m0 = LOG_RUN_START_RE.match(line)
        if m0:
            current_run_num = int(m0.group("n"))
            continue
        m1 = LOG_RUN_RE.match(line)
        if m1:
            # find error message on a subsequent "Run #N error:" line — approximate
            runs.append({
                "run_num": current_run_num if current_run_num is not None else len(runs) + 1,
                "timestamp": m1.group("ts"),
                "session_id": m1.group("sid"),
                "error": m1.group("err") == "True",
                "turns": int(m1.group("turns")),
                "cost_usd": float(m1.group("cost")),
                "error_message": None,
            })
    # Scan for explicit error lines to attach to the most recent run
    for i, line in enumerate(lines):
        if "] [agent-sdk] Run #" in line and " error:" in line:
            if runs:
                msg = line.split(" error:", 1)[-1].strip()
                # find which run this belonged to
                m = re.search(r"Run #(\d+) error", line)
                if m:
                    n = int(m.group(1))
                    for r in reversed(runs):
                        if r.get("run_num") == n:
                            r["error_message"] = msg[:240]
                            break
    return runs[-max_runs:]


# ----------------------------------------------------------------------
# Aggregate state — called per request
# ----------------------------------------------------------------------


def aggregate_state() -> dict:
    ts = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    out: dict = {"timestamp": ts, "error": None}

    if not COMPUTER_ID or not ORB_API_KEY:
        out["error"] = "missing ORB_API_KEY or ORB_COMPUTER_ID"
        return out

    try:
        out["computer"] = _orb_get(f"/v1/computers/{COMPUTER_ID}")
    except Exception as exc:
        out["error"] = f"computer fetch failed: {exc}"
        out["computer"] = None

    try:
        out["agents"] = (
            _orb_get(f"/v1/computers/{COMPUTER_ID}/agents") or {}
        ).get("agents", [])
    except Exception:
        out["agents"] = []

    # Logs
    log_text = _orb_file_text("agent/code/logs/agent-sdk.log")
    out["recent_log"] = "\n".join(log_text.splitlines()[-200:])
    out["runs"] = parse_runs(log_text)

    # Env snapshot written by agent.py diag
    snap = _orb_file_json("agent/code/logs/env-snapshot.json") or {}
    # Normalize: we only need the keys the dashboard cares about + "<set>"/"<MISSING>"
    env_out: dict[str, str] = {}
    for k in (
        "ANTHROPIC_BASE_URL", "ANTHROPIC_AUTH_TOKEN", "ANTHROPIC_API_KEY",
        "API_TIMEOUT_MS",
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

    # Findings counts (directory listing)
    findings = {}
    for name in ("candidates", "docking", "admet", "novelty", "mechanism",
                 "red-team", "retrosynthesis", "weekly-reports"):
        items = _orb_dir(f"agent/code/findings/{name}")
        # Count directories for candidates (one dir per candidate),
        # files otherwise. Skip .gitkeep.
        if name == "candidates":
            n = sum(1 for f in items if f.get("type") == "directory" and f.get("name") != "candidates")
        else:
            n = sum(1 for f in items if f.get("type") == "file" and f.get("name") != ".gitkeep")
        key = name.replace("-", "_")
        findings[key] = n
    out["findings"] = findings

    # Leaderboard
    out["leaderboard"] = _orb_file_json("agent/code/findings/leaderboard.json") or {
        "candidates": [], "updated_at": None
    }

    # Loop health
    pc = _orb_file_jsonl("agent/code/findings/loop-health/positive-control.jsonl")
    nc = _orb_file_jsonl("agent/code/findings/loop-health/negative-control.jsonl")
    pa = _orb_file_json("agent/code/findings/loop-health/positive-control-alert.json")
    na = _orb_file_json("agent/code/findings/loop-health/negative-control-alert.json")
    out["loop_health"] = {
        "positive_control": pc,
        "negative_control": nc,
        "positive_alert": pa,
        "negative_alert": na,
    }

    return out


# ----------------------------------------------------------------------
# HTTP server
# ----------------------------------------------------------------------


class Handler(BaseHTTPRequestHandler):
    def _send_json(self, obj, code: int = 200) -> None:
        body = json.dumps(obj).encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.send_header("Cache-Control", "no-store")
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _send_file(self, path: Path, mime: str = "text/html; charset=utf-8") -> None:
        try:
            body = path.read_bytes()
        except FileNotFoundError:
            self.send_error(404, "not found")
            return
        self.send_response(200)
        self.send_header("Content-Type", mime)
        self.send_header("Cache-Control", "no-store")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self) -> None:  # noqa: N802
        parsed = urllib.parse.urlparse(self.path)
        if parsed.path in ("/", "/index.html", "/monitor.html"):
            self._send_file(HTML)
            return
        if parsed.path == "/api/state":
            try:
                state = aggregate_state()
                self._send_json(state)
            except Exception as exc:
                self._send_json(
                    {"error": f"{type(exc).__name__}: {exc}",
                     "trace": traceback.format_exc()},
                    code=500,
                )
            return
        if parsed.path == "/healthz":
            self._send_json({"ok": True})
            return
        self.send_error(404, "not found")

    def log_message(self, fmt: str, *args) -> None:  # noqa: N802
        # quieter default logs
        sys.stderr.write(f"[monitor] {self.address_string()} {fmt % args}\n")


def main() -> int:
    if not ORB_API_KEY:
        print("[monitor] warning: no ORB_API_KEY found (set ORB_API_KEY or place creds in ~/.orb/credentials)", file=sys.stderr)
    if not COMPUTER_ID:
        print("[monitor] warning: no computer id found (set ORB_COMPUTER_ID or place one in .orb-state/computer-id)", file=sys.stderr)

    httpd = ThreadingHTTPServer((DEFAULT_HOST, DEFAULT_PORT), Handler)
    print(f"[monitor] dashboard at http://{DEFAULT_HOST}:{DEFAULT_PORT}/")
    print(f"[monitor] watching computer {COMPUTER_ID} via {ORB_BASE}")
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass
    return 0


if __name__ == "__main__":
    sys.exit(main())
