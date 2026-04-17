#!/usr/bin/env python3
"""
Weekly loop: cross-method consensus.

Docking scoring is notoriously method-dependent: Vina, DiffDock, GNINA,
and Glide can disagree by several kcal/mol on the same ligand in the
same pocket. A trustworthy candidate should be a consensus winner, not
a single-method artefact.

This loop takes the current top-N leaderboard, re-docks each with a
secondary engine (DiffDock / GNINA / Vina-with-different-seed), and
computes agreement metrics:

  - Spearman rank correlation between primary and secondary ranks.
  - Mean pose RMSD between primary and secondary top poses (Å).
  - Per-candidate disagreement matrix (rank delta > N triggers a flag).

Gate:
  - Spearman correlation ≥ 0.5
  - Mean pose RMSD ≤ 3.0 Å

If a candidate disagrees by more than 5 rank positions AND has RMSD > 4 Å,
it is marked "consensus-unstable" on the leaderboard (not removed, but
flagged).

Pure logic implementation; CI-safe. Integration hook (run_secondary_pass)
is separate.

CLI:
  python -m loops.consensus --target targets/mrsa-gyrb.yaml --evaluate
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
HISTORY_JSONL = LOOP_HEALTH_DIR / "consensus.jsonl"
ALERT_JSON = LOOP_HEALTH_DIR / "consensus-alert.json"

DEFAULT_SPEARMAN_THRESHOLD = 0.5
DEFAULT_MEAN_RMSD_THRESHOLD = 3.0
DEFAULT_RANK_DELTA_UNSTABLE = 5
DEFAULT_RMSD_UNSTABLE = 4.0


@dataclass
class ConsensusEntry:
    candidate_id: str
    primary_score: float
    secondary_score: float
    primary_rank: int
    secondary_rank: int
    rank_delta: int
    pose_rmsd_a: float | None


@dataclass
class ConsensusReport:
    timestamp: str
    target_id: str
    n_candidates: int
    spearman: float
    mean_rmsd_a: float | None
    spearman_threshold: float
    rmsd_threshold: float
    spearman_pass: bool
    rmsd_pass: bool
    passed: bool
    unstable_candidates: list[str] = field(default_factory=list)
    entries: list[dict[str, Any]] = field(default_factory=list)


# ----------------------------------------------------------------------
# Rank correlation (Spearman)
# ----------------------------------------------------------------------


def _rank_asc(values: list[float]) -> list[float]:
    """Average-rank order ascending (ties share mean rank)."""
    indexed = sorted(enumerate(values), key=lambda p: p[1])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(indexed):
        j = i
        while j + 1 < len(indexed) and indexed[j + 1][1] == indexed[i][1]:
            j += 1
        avg = (i + j) / 2 + 1.0
        for k in range(i, j + 1):
            ranks[indexed[k][0]] = avg
        i = j + 1
    return ranks


def spearman_correlation(xs: list[float], ys: list[float]) -> float:
    if len(xs) != len(ys) or len(xs) < 2:
        return 0.0
    rx = _rank_asc(xs)
    ry = _rank_asc(ys)
    mean_rx = statistics.fmean(rx)
    mean_ry = statistics.fmean(ry)
    cov = sum((a - mean_rx) * (b - mean_ry) for a, b in zip(rx, ry))
    var_x = sum((a - mean_rx) ** 2 for a in rx)
    var_y = sum((b - mean_ry) ** 2 for b in ry)
    denom = math.sqrt(var_x * var_y)
    if denom == 0:
        return 0.0
    return cov / denom


# ----------------------------------------------------------------------
# Consensus computation
# ----------------------------------------------------------------------


def compute_consensus(
    *,
    primary: dict[str, float],
    secondary: dict[str, float],
    rmsds: dict[str, float] | None = None,
    target_id: str = "",
    spearman_threshold: float = DEFAULT_SPEARMAN_THRESHOLD,
    rmsd_threshold: float = DEFAULT_MEAN_RMSD_THRESHOLD,
    rank_delta_unstable: int = DEFAULT_RANK_DELTA_UNSTABLE,
    rmsd_unstable: float = DEFAULT_RMSD_UNSTABLE,
) -> ConsensusReport:
    """Compute the weekly consensus report.

    Parameters
    ----------
    primary, secondary : {candidate_id: ΔG (kcal/mol)}
    rmsds              : optional {candidate_id: pose_rmsd_a} between the
                         primary and secondary top poses.
    """
    ids = [c for c in primary.keys() if c in secondary]
    p_scores = [primary[c] for c in ids]
    s_scores = [secondary[c] for c in ids]

    spearman = spearman_correlation(p_scores, s_scores)

    rmsd_list = []
    if rmsds:
        rmsd_list = [rmsds[c] for c in ids if c in rmsds and rmsds[c] is not None]
    mean_rmsd = statistics.fmean(rmsd_list) if rmsd_list else None

    p_ranks = _rank_asc(p_scores)
    s_ranks = _rank_asc(s_scores)

    entries: list[ConsensusEntry] = []
    unstable: list[str] = []
    for idx, cid in enumerate(ids):
        rdelta = int(abs(p_ranks[idx] - s_ranks[idx]))
        e = ConsensusEntry(
            candidate_id=cid,
            primary_score=p_scores[idx],
            secondary_score=s_scores[idx],
            primary_rank=int(p_ranks[idx]),
            secondary_rank=int(s_ranks[idx]),
            rank_delta=rdelta,
            pose_rmsd_a=rmsds.get(cid) if rmsds else None,
        )
        entries.append(e)
        if (e.rank_delta >= rank_delta_unstable
                and e.pose_rmsd_a is not None
                and e.pose_rmsd_a > rmsd_unstable):
            unstable.append(cid)

    report = ConsensusReport(
        timestamp=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        target_id=target_id,
        n_candidates=len(ids),
        spearman=round(spearman, 4),
        mean_rmsd_a=round(mean_rmsd, 3) if mean_rmsd is not None else None,
        spearman_threshold=spearman_threshold,
        rmsd_threshold=rmsd_threshold,
        spearman_pass=spearman >= spearman_threshold,
        rmsd_pass=(mean_rmsd is None) or (mean_rmsd <= rmsd_threshold),
        passed=(spearman >= spearman_threshold)
        and ((mean_rmsd is None) or (mean_rmsd <= rmsd_threshold)),
        unstable_candidates=unstable,
        entries=[asdict(e) for e in entries],
    )
    return report


# ----------------------------------------------------------------------
# History + alert
# ----------------------------------------------------------------------


def append_report(report: ConsensusReport, path: Path = HISTORY_JSONL) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as fh:
        fh.write(json.dumps(asdict(report)) + "\n")


def build_alert(report: ConsensusReport) -> dict[str, Any] | None:
    if report.passed and not report.unstable_candidates:
        return None
    reasons: list[str] = []
    if not report.spearman_pass:
        reasons.append(
            f"Spearman {report.spearman:.3f} < threshold {report.spearman_threshold:.2f}"
        )
    if not report.rmsd_pass:
        reasons.append(
            f"mean pose RMSD {report.mean_rmsd_a:.2f} Å > "
            f"threshold {report.rmsd_threshold:.2f}"
        )
    return {
        "timestamp": report.timestamp,
        "target_id": report.target_id,
        "severity": "warn",  # not a hard halt — re-evaluate before acting
        "reasons": reasons,
        "unstable_candidates": report.unstable_candidates,
        "action": (
            "Cross-method docking disagrees. Verify the secondary engine's "
            "configuration. If primary vs secondary disagree systematically, "
            "candidates in the top-20 should be tagged 'consensus-unstable' "
            "on the leaderboard until resolved."
        ),
    }


# ----------------------------------------------------------------------
# Integration hook (not exercised in CI)
# ----------------------------------------------------------------------


def run_secondary_pass(
    target_yaml: Path,
    leaderboard_path: Path,
    top_n: int = 20,
    *,
    secondary_engine: str = "vina_seed2",
) -> ConsensusReport:  # pragma: no cover - integration only
    """Integration-only: re-dock the current leaderboard top-N with a
    secondary engine. For production wire DiffDock or GNINA here."""
    sys.path.insert(0, str(REPO_ROOT / "src"))
    from tools.dock_vina import VinaDocker  # type: ignore

    board = json.loads(Path(leaderboard_path).read_text())
    candidates = board.get("candidates", [])[:top_n]

    primary = {c["candidate_id"]: c["sub_scores"]["docking"] for c in candidates}
    secondary: dict[str, float] = {}

    docker = VinaDocker(target_yaml)
    docker.seed = 2  # deterministic "secondary" via different seed
    for c in candidates:
        smi = c["smiles"]
        r = docker.dock(smi, f"consensus-{c['candidate_id']}")
        secondary[c["candidate_id"]] = r.best_energy_kcalmol or 0.0

    return compute_consensus(
        primary=primary, secondary=secondary, rmsds=None,
        target_id=yaml.safe_load(target_yaml.read_text()).get("id", ""),
    )


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------


def _main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--target", type=Path, required=False)
    ap.add_argument("--evaluate", action="store_true")
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--leaderboard", type=Path,
                    default=REPO_ROOT / "findings" / "leaderboard.json")
    ap.add_argument("--top-n", type=int, default=20)
    args = ap.parse_args()

    if args.run:  # pragma: no cover
        if not args.target:
            print("--run requires --target", file=sys.stderr)
            return 2
        report = run_secondary_pass(args.target, args.leaderboard, top_n=args.top_n)
        append_report(report)
        alert = build_alert(report)
        if alert:
            ALERT_JSON.parent.mkdir(parents=True, exist_ok=True)
            ALERT_JSON.write_text(json.dumps(alert, indent=2))
        print(json.dumps(asdict(report), indent=2))
        return 0 if (alert is None) else 1

    # Evaluate mode: read last history entry
    if not HISTORY_JSONL.exists():
        print(json.dumps({"message": "no consensus history yet"}, indent=2))
        return 0
    lines = [ln for ln in HISTORY_JSONL.read_text().splitlines() if ln.strip()]
    if not lines:
        print(json.dumps({"message": "no consensus history yet"}, indent=2))
        return 0
    data = json.loads(lines[-1])
    print(json.dumps(data, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
