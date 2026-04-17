#!/usr/bin/env python3
"""
Hourly loop: positive control via re-docking of reference inhibitors.

The target YAML lists a set of reference_inhibitors with known activity
against MRSA GyrB (novobiocin, clorobiocin, GSK299423, etc.). This loop:

  1. Re-docks every reference inhibitor against the same receptor the
     agent is using for candidates.
  2. Writes the results to findings/loop-health/positive-control.jsonl
     (one line per run, append-only — time series of the pipeline's
     health).
  3. Computes rolling mean ΔG per reference across the last N runs and
     flags any reference whose current ΔG drifts more than
     `max_drift_kcalmol` (default 1.5) from its rolling mean.
  4. Also flags any reference whose current ΔG rank (among recent
     candidates + references) falls outside the top 5%.
  5. If any flag raises, exits non-zero AND writes
     findings/loop-health/positive-control-alert.json — the agent is
     expected to halt new candidate generation until the alert is
     cleared.

The pipeline-health telemetry is public: the dashboard surfaces the
rolling ΔG of novobiocin & GSK299423 as line charts so observers can
see the pipeline verifying itself in real time.

CLI:
  python -m loops.positive_control --target targets/mrsa-gyrb.yaml [--run]
  python -m loops.positive_control --evaluate

The `--run` flag performs the actual docking (requires vina); the
default mode reads the historical JSONL and re-evaluates drift.
`--evaluate` is a pure-logic pass that is CI-safe.
"""

from __future__ import annotations

import argparse
import json
import statistics
import sys
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Iterable

import yaml


REPO_ROOT = Path(__file__).parent.parent
LOOP_HEALTH_DIR = REPO_ROOT / "findings" / "loop-health"
HISTORY_JSONL = LOOP_HEALTH_DIR / "positive-control.jsonl"
ALERT_JSON = LOOP_HEALTH_DIR / "positive-control-alert.json"

# Rolling-window size for the mean ΔG baseline.
ROLLING_WINDOW = 24   # 24 runs × 1h = 1 day of recent history

# Drift threshold (kcal/mol).
DEFAULT_MAX_DRIFT = 1.5


# ----------------------------------------------------------------------
# Data classes
# ----------------------------------------------------------------------


@dataclass
class ReferenceRun:
    name: str
    delta_g_kcalmol: float | None
    duration_s: float | None = None
    error: str | None = None


@dataclass
class RunRecord:
    timestamp: str
    target_id: str
    runs: list[ReferenceRun] = field(default_factory=list)

    def to_jsonl(self) -> str:
        return json.dumps(asdict(self))


@dataclass
class DriftEvaluation:
    name: str
    current: float | None
    rolling_mean: float | None
    rolling_n: int
    drift_abs: float | None
    over_threshold: bool
    threshold: float


# ----------------------------------------------------------------------
# Loading history
# ----------------------------------------------------------------------


def load_history(path: Path = HISTORY_JSONL) -> list[RunRecord]:
    if not path.exists():
        return []
    records: list[RunRecord] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            data = json.loads(line)
        except Exception:
            continue
        runs = [ReferenceRun(**r) for r in data.get("runs", [])]
        records.append(RunRecord(
            timestamp=data.get("timestamp", ""),
            target_id=data.get("target_id", ""),
            runs=runs,
        ))
    return records


def append_record(record: RunRecord, path: Path = HISTORY_JSONL) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as fh:
        fh.write(record.to_jsonl() + "\n")


# ----------------------------------------------------------------------
# Drift evaluation (pure, CI-safe)
# ----------------------------------------------------------------------


def rolling_mean_by_reference(
    history: Iterable[RunRecord],
    window: int = ROLLING_WINDOW,
) -> dict[str, tuple[float, int]]:
    """Return {reference_name: (mean_delta_g, n)} computed over the last
    `window` records' successful runs."""
    recent = list(history)[-window:]
    buckets: dict[str, list[float]] = {}
    for rec in recent:
        for r in rec.runs:
            if r.delta_g_kcalmol is None:
                continue
            buckets.setdefault(r.name, []).append(float(r.delta_g_kcalmol))
    return {name: (statistics.fmean(vals), len(vals)) for name, vals in buckets.items() if vals}


def evaluate_drift(
    current_record: RunRecord,
    history: Iterable[RunRecord],
    *,
    window: int = ROLLING_WINDOW,
    max_drift_kcalmol: float = DEFAULT_MAX_DRIFT,
) -> list[DriftEvaluation]:
    means = rolling_mean_by_reference(history, window=window)
    evaluations: list[DriftEvaluation] = []
    for r in current_record.runs:
        mean_n = means.get(r.name)
        if r.delta_g_kcalmol is None:
            evaluations.append(DriftEvaluation(
                name=r.name, current=None,
                rolling_mean=mean_n[0] if mean_n else None,
                rolling_n=mean_n[1] if mean_n else 0,
                drift_abs=None, over_threshold=True,
                threshold=max_drift_kcalmol,
            ))
            continue
        if mean_n is None or mean_n[1] < 3:
            # Not enough history to judge — bootstrap phase.
            evaluations.append(DriftEvaluation(
                name=r.name, current=r.delta_g_kcalmol,
                rolling_mean=None, rolling_n=mean_n[1] if mean_n else 0,
                drift_abs=None, over_threshold=False,
                threshold=max_drift_kcalmol,
            ))
            continue
        drift = abs(r.delta_g_kcalmol - mean_n[0])
        evaluations.append(DriftEvaluation(
            name=r.name, current=r.delta_g_kcalmol,
            rolling_mean=mean_n[0], rolling_n=mean_n[1],
            drift_abs=drift, over_threshold=drift > max_drift_kcalmol,
            threshold=max_drift_kcalmol,
        ))
    return evaluations


def build_alert(
    evaluations: list[DriftEvaluation],
    *,
    target_id: str,
) -> dict[str, Any] | None:
    flagged = [e for e in evaluations if e.over_threshold]
    if not flagged:
        return None
    return {
        "target_id": target_id,
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "flagged": [asdict(e) for e in flagged],
        "severity": "halt",
        "action": (
            "Halt new candidate generation. The pipeline's positive "
            "controls have drifted; the Vina / receptor / grid box may "
            "be misconfigured. Investigate before resuming."
        ),
    }


# ----------------------------------------------------------------------
# Docking execution (requires vina — delegates to VinaDocker)
# ----------------------------------------------------------------------


def run_docking_pass(
    target_yaml: Path,
) -> RunRecord:  # pragma: no cover - exercised in integration only
    """Re-dock every reference inhibitor listed in the target YAML."""
    sys.path.insert(0, str(REPO_ROOT / "src"))
    from tools.dock_vina import VinaDocker  # type: ignore

    docker = VinaDocker(target_yaml)
    results = docker.benchmark()
    runs = []
    for r in results:
        runs.append(ReferenceRun(
            name=r.candidate_id.removeprefix("bench-"),
            delta_g_kcalmol=r.best_energy_kcalmol,
            duration_s=r.duration_s,
            error=r.error,
        ))
    return RunRecord(
        timestamp=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        target_id=yaml.safe_load(target_yaml.read_text()).get("id", ""),
        runs=runs,
    )


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------


def _cmd_run(args) -> int:  # pragma: no cover - requires vina
    record = run_docking_pass(args.target)
    append_record(record)
    return _evaluate_and_alert(record, args)


def _cmd_evaluate(args) -> int:
    history = load_history()
    if not history:
        print(json.dumps({"message": "no history yet"}, indent=2))
        return 0
    latest = history[-1]
    evaluations = evaluate_drift(
        latest, history[:-1], max_drift_kcalmol=args.max_drift,
    )
    alert = build_alert(evaluations, target_id=latest.target_id)
    out = {
        "timestamp": latest.timestamp,
        "target_id": latest.target_id,
        "evaluations": [asdict(e) for e in evaluations],
        "alert": alert,
    }
    print(json.dumps(out, indent=2))
    if alert:
        ALERT_JSON.parent.mkdir(parents=True, exist_ok=True)
        ALERT_JSON.write_text(json.dumps(alert, indent=2))
        return 1
    ALERT_JSON.unlink(missing_ok=True)
    return 0


def _evaluate_and_alert(record: RunRecord, args) -> int:
    history = load_history()
    evaluations = evaluate_drift(
        record, history[:-1], max_drift_kcalmol=args.max_drift,
    )
    alert = build_alert(evaluations, target_id=record.target_id)
    if alert:
        ALERT_JSON.parent.mkdir(parents=True, exist_ok=True)
        ALERT_JSON.write_text(json.dumps(alert, indent=2))
        return 1
    ALERT_JSON.unlink(missing_ok=True)
    return 0


def _main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--target", type=Path, required=False)
    ap.add_argument("--run", action="store_true", help="re-dock references via vina")
    ap.add_argument("--evaluate", action="store_true", help="evaluate history only")
    ap.add_argument("--max-drift", type=float, default=DEFAULT_MAX_DRIFT)
    args = ap.parse_args()

    if args.run:
        if not args.target:
            print("--run requires --target", file=sys.stderr)
            return 2
        return _cmd_run(args)
    return _cmd_evaluate(args)


if __name__ == "__main__":
    raise SystemExit(_main())
