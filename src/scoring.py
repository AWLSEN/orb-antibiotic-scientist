#!/usr/bin/env python3
"""
Composite rigor score + leaderboard updater.

Each candidate produces per-layer artefacts (docking, druglikeness,
admet, novelty, mechanism, retrosynthesis, red-team). This module:

  1. Reduces each layer's artefact to a sub-score in [0, 1].
  2. Combines them with weights from the target YAML (`scoring.weights`)
     into a composite rigor score, also in [0, 1].
  3. Merges the candidate into findings/leaderboard.json, keeping the
     top-N by rigor score and filtering below `scoring.leaderboard_threshold`.

Scoring is intentionally simple and deterministic so that `scoring.py`
is reproducible: same artefacts → same rigor score, byte-for-byte.

CLI:
  python -m scoring update --target targets/mrsa-gyrb.yaml --candidate-dir findings/candidates/cand-0001
  python -m scoring show --target targets/mrsa-gyrb.yaml

Library:
  from scoring import score_candidate, update_leaderboard
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import yaml

REPO_ROOT = Path(__file__).parent.parent
LEADERBOARD_PATH = REPO_ROOT / "findings" / "leaderboard.json"


# ----------------------------------------------------------------------
# Per-layer scoring — each returns a float in [0, 1]
# ----------------------------------------------------------------------


def score_docking(docking: dict[str, Any] | None) -> float:
    """Map Vina ΔG (kcal/mol) to [0, 1]. −5 → 0, −11 → 1, linear."""
    if not docking:
        return 0.0
    e = docking.get("best_energy_kcalmol")
    if e is None:
        return 0.0
    try:
        e = float(e)
    except (TypeError, ValueError):
        return 0.0
    return max(0.0, min(1.0, (-e - 5.0) / 6.0))


def score_novelty(novelty: dict[str, Any] | None) -> float:
    """1 − max_tanimoto, with scaffold hits zeroing the score."""
    if not novelty:
        return 0.0
    if novelty.get("scaffold_hits"):
        return 0.0
    tc = novelty.get("max_tanimoto")
    # NaN (offline mode / query failure) ⇒ partial credit so candidates
    # aren't automatically rejected by network flake.
    try:
        tc_f = float(tc)
        if math.isnan(tc_f):
            return 0.5
    except (TypeError, ValueError):
        return 0.5
    return max(0.0, min(1.0, 1.0 - tc_f))


def score_druglikeness(validate: dict[str, Any] | None) -> float:
    """Fraction of Layer-3 gates (Lipinski, Veber, Ghose) that pass."""
    if not validate:
        return 0.0
    layers = validate.get("layers") or []
    l3 = next((l for l in layers if l.get("name") == "druglikeness"), None)
    if not l3:
        return 0.0
    gates = l3.get("gates") or []
    if not gates:
        return 0.0
    return sum(1.0 for g in gates if g.get("passed")) / len(gates)


def score_admet(admet: dict[str, Any] | None) -> float:
    """1 − 0.15 × red_flag_count, clipped to [0, 1]."""
    if not admet:
        return 0.0
    flags = admet.get("red_flags") or []
    return max(0.0, min(1.0, 1.0 - 0.15 * len(flags)))


def score_synthesizability(sa_score: float | int | None) -> float:
    """SA 1.0 → 1.0, SA 6.0 → 0.0, SA > 6.0 → 0."""
    if sa_score is None:
        return 0.5  # unknown
    try:
        sa = float(sa_score)
    except (TypeError, ValueError):
        return 0.5
    return max(0.0, min(1.0, (6.0 - sa) / 5.0))


def score_mechanism(mechanism: dict[str, Any] | None) -> float:
    """0.0 if gate failed, else 0.5 + 0.125 × catalytic hits, capped at 1.0."""
    if not mechanism:
        return 0.0
    if not mechanism.get("passed"):
        return 0.0
    cat_hits = len(mechanism.get("catalytic_residues_hit") or [])
    return min(1.0, 0.5 + 0.125 * cat_hits)


# Legacy / optional layer: red-team isn't in the target YAML weights but
# acts as a veto. If substantive flaw → composite gets zeroed.


# ----------------------------------------------------------------------
# Composite scoring
# ----------------------------------------------------------------------


@dataclass
class SubScores:
    docking: float
    novelty: float
    druglikeness: float
    admet: float
    synthesizability: float
    mechanism: float


@dataclass
class ScoredCandidate:
    candidate_id: str
    target_id: str
    smiles: str
    rigor_score: float
    sub_scores: dict[str, float]
    veto_applied: bool
    veto_reason: str | None
    above_threshold: bool
    leaderboard_threshold: float
    timestamp: str = field(
        default_factory=lambda: time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    )


def _load_target(target_yaml: Path) -> dict[str, Any]:
    return yaml.safe_load(Path(target_yaml).read_text())


def score_candidate(
    *,
    target_yaml: Path,
    candidate_id: str,
    smiles: str,
    docking: dict[str, Any] | None = None,
    novelty: dict[str, Any] | None = None,
    validate: dict[str, Any] | None = None,
    admet: dict[str, Any] | None = None,
    sa_score: float | None = None,
    mechanism: dict[str, Any] | None = None,
    redteam: dict[str, Any] | None = None,
) -> ScoredCandidate:
    target = _load_target(target_yaml)
    weights = target["scoring"]["weights"]
    threshold = float(target["scoring"].get("leaderboard_threshold", 0.83))

    subs = SubScores(
        docking=score_docking(docking),
        novelty=score_novelty(novelty),
        druglikeness=score_druglikeness(validate),
        admet=score_admet(admet),
        synthesizability=score_synthesizability(sa_score),
        mechanism=score_mechanism(mechanism),
    )

    def w(name: str) -> float:
        return float(weights.get(name, 0.0))

    rigor = (
        w("docking_vina") * subs.docking
        + w("novelty") * subs.novelty
        + w("druglikeness") * subs.druglikeness
        + w("admet") * subs.admet
        + w("synthesizability") * subs.synthesizability
        + w("mechanism") * subs.mechanism
    )

    veto_applied = False
    veto_reason: str | None = None
    if redteam and redteam.get("substantive_flaw"):
        veto_applied = True
        veto_reason = "red-team flagged substantive flaw"
        rigor = 0.0

    rigor = round(rigor, 4)
    above = (not veto_applied) and (rigor >= threshold)

    return ScoredCandidate(
        candidate_id=candidate_id,
        target_id=target.get("id", ""),
        smiles=smiles,
        rigor_score=rigor,
        sub_scores={
            "docking": round(subs.docking, 4),
            "novelty": round(subs.novelty, 4),
            "druglikeness": round(subs.druglikeness, 4),
            "admet": round(subs.admet, 4),
            "synthesizability": round(subs.synthesizability, 4),
            "mechanism": round(subs.mechanism, 4),
        },
        veto_applied=veto_applied,
        veto_reason=veto_reason,
        above_threshold=above,
        leaderboard_threshold=threshold,
    )


# ----------------------------------------------------------------------
# Leaderboard persistence
# ----------------------------------------------------------------------


def _load_leaderboard(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"candidates": [], "updated_at": None, "top_n": 100}
    try:
        return json.loads(path.read_text())
    except Exception:
        return {"candidates": [], "updated_at": None, "top_n": 100}


def update_leaderboard(
    scored: ScoredCandidate,
    *,
    leaderboard_path: Path = LEADERBOARD_PATH,
    top_n: int = 100,
) -> dict[str, Any]:
    """Upsert a candidate into the leaderboard and persist.

    - Replaces existing entry for the same candidate_id.
    - Only candidates with above_threshold=True are kept.
    - Sorted by rigor_score descending, truncated at top_n.
    """
    data = _load_leaderboard(leaderboard_path)
    existing: list[dict[str, Any]] = data.get("candidates", [])

    if scored.above_threshold:
        filtered = [e for e in existing if e.get("candidate_id") != scored.candidate_id]
        filtered.append(asdict(scored))
        filtered.sort(key=lambda x: x.get("rigor_score", 0.0), reverse=True)
        candidates = filtered[:top_n]
    else:
        candidates = [e for e in existing if e.get("candidate_id") != scored.candidate_id]

    leaderboard_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "candidates": candidates,
        "updated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "top_n": top_n,
    }
    leaderboard_path.write_text(json.dumps(payload, indent=2))
    return payload


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------


def _load_json_or_none(path: Path) -> dict[str, Any] | None:
    if path.exists():
        return json.loads(path.read_text())
    return None


def _cmd_update(args) -> int:
    cand_dir: Path = args.candidate_dir
    if not cand_dir.exists():
        print(f"candidate dir not found: {cand_dir}", file=sys.stderr)
        return 2

    meta = _load_json_or_none(cand_dir / "candidate.json") or {}
    docking = _load_json_or_none(cand_dir / "docking.json")
    novelty = _load_json_or_none(cand_dir / "novelty.json")
    validate = _load_json_or_none(cand_dir / "validate.json")
    admet = _load_json_or_none(cand_dir / "admet.json")
    mechanism = _load_json_or_none(cand_dir / "mechanism.json")
    redteam = _load_json_or_none(cand_dir / "redteam.json")
    sa_score = None
    if validate:
        layers = validate.get("layers") or []
        l4 = next((l for l in layers if l.get("name") == "synthesizability"), None)
        if l4:
            for g in l4.get("gates", []):
                if g.get("name") == "sa_score":
                    sa_score = g.get("detail", {}).get("sa_score")

    scored = score_candidate(
        target_yaml=args.target,
        candidate_id=meta.get("candidate_id") or cand_dir.name,
        smiles=meta.get("smiles") or "",
        docking=docking,
        novelty=novelty,
        validate=validate,
        admet=admet,
        sa_score=sa_score,
        mechanism=mechanism,
        redteam=redteam,
    )
    (cand_dir / "scored.json").write_text(json.dumps(asdict(scored), indent=2))
    board = update_leaderboard(scored, leaderboard_path=args.leaderboard, top_n=args.top_n)
    print(json.dumps({
        "rigor_score": scored.rigor_score,
        "above_threshold": scored.above_threshold,
        "veto_applied": scored.veto_applied,
        "leaderboard_size": len(board["candidates"]),
    }, indent=2))
    return 0 if scored.above_threshold else 1


def _cmd_show(args) -> int:
    data = _load_leaderboard(args.leaderboard)
    print(json.dumps(data, indent=2))
    return 0


def _main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    u = sub.add_parser("update", help="score a candidate and upsert into leaderboard")
    u.add_argument("--target", type=Path, required=True)
    u.add_argument("--candidate-dir", type=Path, required=True)
    u.add_argument("--leaderboard", type=Path, default=LEADERBOARD_PATH)
    u.add_argument("--top-n", type=int, default=100)
    u.set_defaults(func=_cmd_update)

    s = sub.add_parser("show", help="print leaderboard.json")
    s.add_argument("--target", type=Path, required=True)
    s.add_argument("--leaderboard", type=Path, default=LEADERBOARD_PATH)
    s.set_defaults(func=_cmd_show)

    args = ap.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(_main())
