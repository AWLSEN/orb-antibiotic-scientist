#!/usr/bin/env python3
"""
Hourly loop: negative control via DUD-E decoys.

A docking pipeline that *ranks* is only trustworthy when it separates
known actives from property-matched decoys. This loop:

  1. Docks the reference inhibitors (actives) from targets/<id>.yaml and
     a curated decoy set from data/dud-e/<target>-decoys.smi.
  2. Computes two pipeline-quality metrics over the combined ranked list:
       - ROC-AUC (area under the receiver operating curve)
       - EF1%   (enrichment factor at 1% of the ranked list)
  3. Appends a time-series entry to findings/loop-health/negative-
     control.jsonl and writes a halt-alert if either metric drops below
     the target thresholds (defaults EF1% ≥ 5, ROC-AUC ≥ 0.75).

Rationale: a docking setup that scores random decoys as high as the
known inhibitors cannot be trusted to rank candidates. The enrichment
metric is the canonical virtual-screening benchmark (Mysinger et al.,
DUD-E 2012).

The scoring logic (EF1%, ROC-AUC) is implemented locally so CI tests
can exercise it without sklearn installed.

CLI:
  python -m loops.negative_control --target targets/mrsa-gyrb.yaml --run
  python -m loops.negative_control --evaluate
  python -m loops.negative_control --target ... --decoys path/to/decoys.smi --run
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
import sys
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Iterable

import yaml


REPO_ROOT = Path(__file__).parent.parent
LOOP_HEALTH_DIR = REPO_ROOT / "findings" / "loop-health"
HISTORY_JSONL = LOOP_HEALTH_DIR / "negative-control.jsonl"
ALERT_JSON = LOOP_HEALTH_DIR / "negative-control-alert.json"
DECOYS_DIR = REPO_ROOT / "data" / "dud-e"

DEFAULT_EF1_THRESHOLD = 5.0
DEFAULT_AUC_THRESHOLD = 0.75


# ----------------------------------------------------------------------
# Scoring metrics — pure functions
# ----------------------------------------------------------------------


def enrichment_factor(
    ranked: list[tuple[float, bool]],
    top_fraction: float = 0.01,
) -> float:
    """Compute EF at `top_fraction` from a list of (score, is_active) pairs
    where `score` is the docking ΔG in kcal/mol (more negative = better).
    Actives are ranked higher when their score is more negative.

    EF = (actives_in_top_x / size_of_top_x) / (actives_in_full / N)

    Returns 0.0 if there are no actives in the top fraction.
    """
    if not ranked:
        return 0.0
    # Sort ascending by ΔG (best = lowest first)
    ordered = sorted(ranked, key=lambda p: p[0])
    n = len(ordered)
    n_top = max(1, int(math.ceil(n * top_fraction)))
    top = ordered[:n_top]
    actives_top = sum(1 for _s, a in top if a)
    actives_total = sum(1 for _s, a in ordered if a)
    if actives_total == 0 or n_top == 0:
        return 0.0
    top_active_rate = actives_top / n_top
    overall_rate = actives_total / n
    if overall_rate == 0:
        return 0.0
    return top_active_rate / overall_rate


def roc_auc(ranked: list[tuple[float, bool]]) -> float:
    """Compute AUC-ROC given (score, is_active) pairs. More-negative
    ΔG = more-likely-active, so we negate ΔG when feeding Wilcoxon.

    Uses the Mann-Whitney U formulation: AUC = U / (n_pos × n_neg).
    """
    if not ranked:
        return 0.0
    actives = [p for p in ranked if p[1]]
    decoys = [p for p in ranked if not p[1]]
    if not actives or not decoys:
        return 0.0
    u = 0.0
    for a_score, _ in actives:
        for d_score, _ in decoys:
            if a_score < d_score:          # active wins (more negative ΔG)
                u += 1.0
            elif a_score == d_score:
                u += 0.5
    return u / (len(actives) * len(decoys))


@dataclass
class ControlRun:
    name: str
    delta_g_kcalmol: float | None
    is_active: bool
    error: str | None = None


@dataclass
class RunRecord:
    timestamp: str
    target_id: str
    runs: list[ControlRun] = field(default_factory=list)
    ef1: float | None = None
    roc_auc: float | None = None

    def to_jsonl(self) -> str:
        return json.dumps(asdict(self))


@dataclass
class EvaluationResult:
    timestamp: str
    target_id: str
    ef1: float
    roc_auc: float
    ef1_threshold: float
    auc_threshold: float
    ef1_pass: bool
    auc_pass: bool
    passed: bool


# ----------------------------------------------------------------------
# History + alerts
# ----------------------------------------------------------------------


def append_record(record: RunRecord, path: Path = HISTORY_JSONL) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as fh:
        fh.write(record.to_jsonl() + "\n")


def load_history(path: Path = HISTORY_JSONL) -> list[RunRecord]:
    if not path.exists():
        return []
    out: list[RunRecord] = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        try:
            data = json.loads(line)
        except Exception:
            continue
        runs = [ControlRun(**r) for r in data.get("runs", [])]
        out.append(RunRecord(
            timestamp=data.get("timestamp", ""),
            target_id=data.get("target_id", ""),
            runs=runs,
            ef1=data.get("ef1"),
            roc_auc=data.get("roc_auc"),
        ))
    return out


def build_alert(evaluation: EvaluationResult) -> dict[str, Any] | None:
    if evaluation.passed:
        return None
    reasons = []
    if not evaluation.ef1_pass:
        reasons.append(
            f"EF1% {evaluation.ef1:.2f} < threshold {evaluation.ef1_threshold:.2f}"
        )
    if not evaluation.auc_pass:
        reasons.append(
            f"ROC-AUC {evaluation.roc_auc:.3f} < threshold {evaluation.auc_threshold:.3f}"
        )
    return {
        "timestamp": evaluation.timestamp,
        "target_id": evaluation.target_id,
        "severity": "halt",
        "reasons": reasons,
        "ef1": evaluation.ef1,
        "roc_auc": evaluation.roc_auc,
        "action": (
            "Pipeline fails actives-vs-decoys discrimination. Halt "
            "candidate generation. Investigate docking setup before "
            "resuming."
        ),
    }


def evaluate_record(
    record: RunRecord,
    *,
    ef1_threshold: float = DEFAULT_EF1_THRESHOLD,
    auc_threshold: float = DEFAULT_AUC_THRESHOLD,
) -> EvaluationResult:
    scored = [
        (r.delta_g_kcalmol, r.is_active)
        for r in record.runs
        if r.delta_g_kcalmol is not None
    ]
    ef1 = enrichment_factor(scored, top_fraction=0.01)
    auc = roc_auc(scored)
    return EvaluationResult(
        timestamp=record.timestamp,
        target_id=record.target_id,
        ef1=ef1,
        roc_auc=auc,
        ef1_threshold=ef1_threshold,
        auc_threshold=auc_threshold,
        ef1_pass=ef1 >= ef1_threshold,
        auc_pass=auc >= auc_threshold,
        passed=(ef1 >= ef1_threshold) and (auc >= auc_threshold),
    )


# ----------------------------------------------------------------------
# Decoy loading
# ----------------------------------------------------------------------


def load_decoys(target_id: str, explicit_path: Path | None = None) -> list[str]:
    """Load decoy SMILES, one per line, from either --decoys or
    data/dud-e/<target_id>-decoys.smi. Returns [] if not found."""
    path = explicit_path or (DECOYS_DIR / f"{target_id}-decoys.smi")
    if not path.exists():
        return []
    lines = path.read_text().splitlines()
    decoys: list[str] = []
    for line in lines:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        decoys.append(line.split()[0])  # first whitespace-delimited token
    return decoys


# ----------------------------------------------------------------------
# Docking execution (integration only)
# ----------------------------------------------------------------------


def run_docking_pass(
    target_yaml: Path,
    decoys: list[str],
) -> RunRecord:  # pragma: no cover - integration only
    sys.path.insert(0, str(REPO_ROOT / "src"))
    from tools.dock_vina import VinaDocker  # type: ignore

    docker = VinaDocker(target_yaml)
    target = yaml.safe_load(target_yaml.read_text())

    runs: list[ControlRun] = []
    # Actives
    for ref in target.get("reference_inhibitors", []):
        smi = ref.get("smiles", "").strip()
        if not smi:
            continue
        r = docker.dock(smi, f"active-{ref['name'].replace(' ', '_').lower()}")
        runs.append(ControlRun(
            name=ref["name"], delta_g_kcalmol=r.best_energy_kcalmol,
            is_active=True, error=r.error,
        ))
    # Decoys
    for i, smi in enumerate(decoys):
        r = docker.dock(smi, f"decoy-{i:04d}")
        runs.append(ControlRun(
            name=f"decoy-{i:04d}", delta_g_kcalmol=r.best_energy_kcalmol,
            is_active=False, error=r.error,
        ))

    return RunRecord(
        timestamp=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        target_id=target.get("id", ""),
        runs=runs,
    )


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------


def _cmd_run(args) -> int:  # pragma: no cover - requires vina
    target = yaml.safe_load(args.target.read_text())
    decoys = load_decoys(target.get("id", ""), args.decoys)
    if not decoys:
        print(f"no decoys found for {target.get('id')}", file=sys.stderr)
        return 2
    record = run_docking_pass(args.target, decoys)
    evaluation = evaluate_record(
        record, ef1_threshold=args.ef1_threshold, auc_threshold=args.auc_threshold,
    )
    record.ef1 = evaluation.ef1
    record.roc_auc = evaluation.roc_auc
    append_record(record)
    alert = build_alert(evaluation)
    if alert:
        ALERT_JSON.parent.mkdir(parents=True, exist_ok=True)
        ALERT_JSON.write_text(json.dumps(alert, indent=2))
        return 1
    ALERT_JSON.unlink(missing_ok=True)
    return 0


def _cmd_evaluate(args) -> int:
    history = load_history()
    if not history:
        print(json.dumps({"message": "no history yet"}, indent=2))
        return 0
    latest = history[-1]
    evaluation = evaluate_record(
        latest, ef1_threshold=args.ef1_threshold, auc_threshold=args.auc_threshold,
    )
    alert = build_alert(evaluation)
    print(json.dumps({
        "evaluation": asdict(evaluation),
        "alert": alert,
    }, indent=2))
    return 0 if alert is None else 1


def _main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--target", type=Path)
    ap.add_argument("--decoys", type=Path, default=None)
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--evaluate", action="store_true")
    ap.add_argument("--ef1-threshold", type=float, default=DEFAULT_EF1_THRESHOLD)
    ap.add_argument("--auc-threshold", type=float, default=DEFAULT_AUC_THRESHOLD)
    args = ap.parse_args()

    if args.run:
        if not args.target:
            print("--run requires --target", file=sys.stderr)
            return 2
        return _cmd_run(args)
    return _cmd_evaluate(args)


if __name__ == "__main__":
    raise SystemExit(_main())
