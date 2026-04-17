#!/usr/bin/env python3
"""
orb-antibiotic-scientist -- Autonomous AI Medicinal Chemist
Uses the Claude Agent SDK for persistent, auto-compacting sessions.
Runs indefinitely. Auto-compacts when context fills. Never stops.

Mission: design and computationally validate novel antibiotic candidates
for drug-resistant bacterial targets (initial target: MRSA DNA gyrase B).

This is the ported harness pattern from SPOQ-Food, adapted for the
medicinal-chemistry domain and the verification architecture described
in README.md.
"""

import asyncio
import json
import os
import time
import traceback
from pathlib import Path

from claude_agent_sdk import (
    query, ClaudeAgentOptions,
    AssistantMessage, ResultMessage, SystemMessage,
    TextBlock, ToolUseBlock,
)

PROJECT_DIR = Path(__file__).parent.parent
PROMPT_FILE = PROJECT_DIR / "src" / "agent-prompt.md"
LOG_DIR = PROJECT_DIR / "logs"
PID_DIR = LOG_DIR / "pids"

LOG_DIR.mkdir(parents=True, exist_ok=True)
PID_DIR.mkdir(parents=True, exist_ok=True)

FINDING_DIRS = [
    "findings/candidates",
    "findings/docking",
    "findings/admet",
    "findings/novelty",
    "findings/mechanism",
    "findings/red-team",
    "findings/retrosynthesis",
    "findings/weekly-reports",
    "findings/loop-health",
    "data",
]

for d in FINDING_DIRS:
    (PROJECT_DIR / d).mkdir(parents=True, exist_ok=True)


def count_files(directory, pattern="*"):
    d = PROJECT_DIR / directory
    if not d.exists():
        return 0
    return len([f for f in d.glob(pattern) if f.is_file()])


def get_stats():
    return {
        "candidates": count_files("findings/candidates"),
        "docking_runs": count_files("findings/docking"),
        "admet_reports": count_files("findings/admet"),
        "novelty_reports": count_files("findings/novelty"),
        "mechanism_reports": count_files("findings/mechanism"),
        "red_team_critiques": count_files("findings/red-team"),
        "retrosynthesis_plans": count_files("findings/retrosynthesis"),
        "weekly_reports": count_files("findings/weekly-reports"),
    }


def log(msg):
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] [agent-sdk] {msg}"
    print(line, flush=True)
    with open(LOG_DIR / "agent-sdk.log", "a") as f:
        f.write(line + "\n")


SESSION_FILE = LOG_DIR / "last_session.txt"


def save_session(sid: str | None):
    if sid:
        SESSION_FILE.write_text(sid)


def load_session() -> str | None:
    if SESSION_FILE.exists():
        sid = SESSION_FILE.read_text().strip()
        return sid or None
    return None


def clear_session():
    SESSION_FILE.unlink(missing_ok=True)


def fmt_stats(s):
    return (
        f"Candidates: {s['candidates']} | "
        f"Docking: {s['docking_runs']} | "
        f"ADMET: {s['admet_reports']} | "
        f"Novelty: {s['novelty_reports']} | "
        f"Mechanism: {s['mechanism_reports']} | "
        f"RedTeam: {s['red_team_critiques']} | "
        f"Retro: {s['retrosynthesis_plans']}"
    )


async def run_agent():
    system_prompt = PROMPT_FILE.read_text() if PROMPT_FILE.exists() else (
        "You are orb-antibiotic-scientist. The real system prompt is missing; "
        "refuse to produce candidates and tell the operator to restore "
        "src/agent-prompt.md."
    )

    # Forward LLM-routing env vars explicitly to the Claude Code subprocess
    # spawned by claude-agent-sdk. Orb's org_secrets reach *this* Python
    # process, but the subprocess does not automatically inherit them in
    # every deployment (confirmed by a diagnostic snapshot on 2026-04-17).
    # ClaudeAgentOptions.env is the canonical plumbing for this.
    _llm_env = {
        k: v for k, v in (
            (key, os.environ.get(key))
            for key in (
                "ANTHROPIC_AUTH_TOKEN", "ANTHROPIC_API_KEY",
                "ANTHROPIC_BASE_URL", "API_TIMEOUT_MS",
                "ANTHROPIC_DEFAULT_OPUS_MODEL",
                "ANTHROPIC_DEFAULT_SONNET_MODEL",
                "ANTHROPIC_DEFAULT_HAIKU_MODEL",
                "PATH", "HOME", "LANG", "LC_ALL",
            )
        )
        if v is not None
    }
    log(f"forwarding {len(_llm_env)} env vars to Claude Code subprocess; "
        f"ANTHROPIC_BASE_URL={_llm_env.get('ANTHROPIC_BASE_URL','<missing>')}")

    session_id = load_session()
    if session_id:
        log(f"Loaded persisted session: {session_id}")
    run_num = 0

    while True:
        run_num += 1
        stats = get_stats()

        log(f"=== RUN #{run_num} | {fmt_stats(stats)} ===")

        if session_id:
            prompt = (
                f"You are orb-antibiotic-scientist continuing your research. "
                f"Run #{run_num}. "
                f"Current stats -- {fmt_stats(stats)}. "
                f"Your previous session was auto-compacted. "
                f"Check findings/ and targets/ to see what exists. "
                f"Continue where you left off. Do not repeat completed work. "
                f"Every candidate must pass the full 12-layer verification chain "
                f"before landing on the leaderboard. Push forward."
            )
            log(f"Resuming session {session_id}...")
            options = ClaudeAgentOptions(
                allowed_tools=[
                    "Bash", "Edit", "Read", "Write", "Glob", "Grep",
                    "WebFetch", "WebSearch",
                ],
                permission_mode="bypassPermissions",
                model="claude-opus-4-7",
                resume=session_id,
                cwd=str(PROJECT_DIR),
                env=_llm_env,
            )
        else:
            prompt = (
                f"You are orb-antibiotic-scientist, an autonomous AI medicinal "
                f"chemist. Run #{run_num}. "
                f"Current stats -- {fmt_stats(stats)}. "
                f"Check findings/ and targets/ to see what's already done. "
                f"DO NOT repeat completed work. "
                f"Priority: design a novel candidate (small molecule or antimicrobial "
                f"peptide) for the active target in targets/, then run the 12-layer "
                f"verification chain (SMILES validity, PAINS/REOS, drug-likeness, "
                f"SA Score, Vina docking, DiffDock consensus, PLIP mechanism, ADMET, "
                f"ChEMBL novelty, resistance-proof, retrosynthesis, red-team critique). "
                f"Save each candidate to findings/candidates/<id>/ with all artefacts. "
                f"Update findings/leaderboard.json. Keep going."
            )
            log("Starting FRESH session...")
            options = ClaudeAgentOptions(
                allowed_tools=[
                    "Bash", "Edit", "Read", "Write", "Glob", "Grep",
                    "WebFetch", "WebSearch",
                ],
                permission_mode="bypassPermissions",
                model="claude-opus-4-7",
                system_prompt=system_prompt,
                cwd=str(PROJECT_DIR),
                env=_llm_env,
            )

        try:
            log_file = LOG_DIR / f"agent_run_{run_num}.log"
            msg_count = 0
            with open(log_file, "a") as lf:
                async for message in query(prompt=prompt, options=options):
                    msg_count += 1
                    lf.write(f"[msg#{msg_count}] {type(message).__name__}: {repr(message)[:500]}\n")
                    lf.flush()
                    if msg_count <= 5:
                        log(f"msg#{msg_count} type={type(message).__name__}")
                    if isinstance(message, ResultMessage):
                        session_id = message.session_id
                        save_session(session_id)
                        log(
                            f"Run complete. Session: {session_id}, "
                            f"Error: {message.is_error}, "
                            f"Turns: {message.num_turns}, "
                            f"Cost: ${message.total_cost_usd or 0:.2f}"
                        )
                        if message.result:
                            lf.write(f"\n--- RESULT ---\n{message.result}\n")
                            lf.flush()
                    elif isinstance(message, AssistantMessage):
                        for block in message.content:
                            if isinstance(block, ToolUseBlock):
                                lf.write(f"[tool] {block.name}: {repr(block.input)[:150]}\n")
                                lf.flush()
                            elif isinstance(block, TextBlock):
                                text = block.text.strip()
                                if text:
                                    lf.write(f"[text] {text[:200]}\n")
                                    lf.flush()
                    elif isinstance(message, SystemMessage):
                        log(f"System event: {message.subtype}")
                        if message.data and "session_id" in message.data:
                            session_id = message.data["session_id"]
                            save_session(session_id)
                            log(f"Session ID from system: {session_id}")

        except KeyboardInterrupt:
            log("Interrupted by user. Stopping.")
            break
        except Exception as e:
            log(f"Run #{run_num} error: {e}")
            log(traceback.format_exc())
            session_id = None
            clear_session()

        new_stats = get_stats()
        log(f"Run #{run_num} ended. {fmt_stats(new_stats)}")
        log("Restarting in 30s...")
        await asyncio.sleep(30)


def main():
    pid_file = PID_DIR / "agent-sdk.pid"
    pid_file.write_text(str(os.getpid()))

    log("============================================")
    log("orb-antibiotic-scientist -- Claude Agent SDK")
    log("Auto-compaction enabled. Runs indefinitely.")
    log("============================================")

    # Diagnostic: dump env-var presence + Orb proxy URL so we can see
    # what actually reached this process. Claude Code inherits our env,
    # so anything missing here explains any downstream auth failures.
    _tracked = [
        "ANTHROPIC_AUTH_TOKEN", "ANTHROPIC_API_KEY",
        "ANTHROPIC_BASE_URL", "OPENAI_BASE_URL", "ORB_PROXY_URL",
        "API_TIMEOUT_MS",
        "ANTHROPIC_DEFAULT_OPUS_MODEL",
        "ANTHROPIC_DEFAULT_SONNET_MODEL",
        "ANTHROPIC_DEFAULT_HAIKU_MODEL",
        "ORB_TASK",
    ]
    _report = {
        k: (
            "<set>" if (k == "ANTHROPIC_AUTH_TOKEN" or k == "ANTHROPIC_API_KEY") and os.environ.get(k)
            else os.environ.get(k) or "<MISSING>"
        )
        for k in _tracked
    }
    log(f"env report: {_report}")

    # Full env snapshot written early so the dashboard can read it even
    # if the agent later fails. Secrets redacted.
    try:
        LOG_DIR.mkdir(parents=True, exist_ok=True)
        snap_path = LOG_DIR / "env-snapshot.json"
        allowlist_prefixes = (
            "ORB_", "PYTHON", "ANTHROPIC_BASE_URL", "ANTHROPIC_DEFAULT",
            "OPENAI_BASE_URL", "GOOGLE_API_BASE_URL", "LLM_BASE_URL",
            "API_TIMEOUT_MS",
        )
        allowlist_exact = ("PATH", "HOME", "LANG", "LC_ALL", "PWD", "USER", "SHELL")
        redact_exact = ("ANTHROPIC_AUTH_TOKEN", "ANTHROPIC_API_KEY")
        snap = {}
        for k, v in os.environ.items():
            if k in redact_exact:
                snap[k] = "<set>"
            elif k.startswith(allowlist_prefixes) or k in allowlist_exact:
                snap[k] = v
            else:
                snap[k] = "<redacted>"
        # Add sentinel so the dashboard can distinguish "file exists but
        # key absent" from "file doesn't exist".
        snap["__snapshot_written_at__"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        snap_path.write_text(json.dumps(snap, indent=2, sort_keys=True))
        log(f"env-snapshot.json written ({len(snap)} keys)")
    except Exception as exc:
        log(f"env-snapshot write failed: {exc}")
        log(traceback.format_exc())

    try:
        asyncio.run(run_agent())
    except KeyboardInterrupt:
        log("Shutdown.")
    finally:
        pid_file.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
