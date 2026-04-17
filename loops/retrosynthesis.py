#!/usr/bin/env python3
"""
Layer 11 + per-top-candidate retrosynthesis feasibility.

A brilliant computational candidate that cannot be synthesised is useless
for wet-lab follow-up. This loop gates candidates by retrosynthetic
feasibility.

Two backends:
  - preferred: AiZynthFinder (https://github.com/MolecularAI/aizynthfinder)
    returns actual route graphs with step counts. Production deploys
    should configure it with a template library + expansion model and
    enable this path.
  - fallback: a deterministic heuristic built on RDKit's BRICS bond
    dissection + Ertl SA score. Coarse but always available; used in CI.

Gate (from target YAML `synthesizability`):
  - SA score ≤ max_sa_score (default 6.0)
  - number of BRICS disconnections ≤ retrosynthesis_max_steps (default 4)
    The BRICS count is the count of retrosynthetically relevant disconnections
    returned by rdkit.Chem.BRICS.BRICSDecompose — a reasonable proxy for
    "how many non-trivial bond-forming steps."

When AiZynthFinder is available, the actual number of steps in the best
route replaces the BRICS heuristic.

CLI:
  python -m loops.retrosynthesis --target targets/mrsa-gyrb.yaml --smiles "..."
  python -m loops.retrosynthesis --target ... --smiles "..." --aizynth

Exit codes: 0 pass, 1 fail, 2 setup error.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import yaml

from rdkit import Chem, RDLogger
from rdkit.Chem import BRICS

RDLogger.DisableLog("rdApp.*")

REPO_ROOT = Path(__file__).parent.parent


@dataclass
class RetroResult:
    input_smiles: str
    canonical_smiles: str
    target_id: str
    backend: str
    sa_score: float | None
    sa_score_cap: float
    steps_estimate: int | None
    steps_cap: int
    brics_fragments: list[str] = field(default_factory=list)
    aizynth_routes: list[dict[str, Any]] = field(default_factory=list)
    passed: bool = False
    reasons: list[str] = field(default_factory=list)
    error: str | None = None


# ----------------------------------------------------------------------
# Heuristic backend — BRICS + SA score
# ----------------------------------------------------------------------


def _sa_score(mol: Chem.Mol) -> float | None:
    """Import RDKit's sascorer lazily and compute Ertl SA score."""
    import os, sys as _sys
    try:
        from rdkit.Chem import RDConfig
        contrib = os.path.join(RDConfig.RDContribDir, "SA_Score")
        if contrib not in _sys.path:
            _sys.path.append(contrib)
        import sascorer  # type: ignore[import-not-found]
        return float(sascorer.calculateScore(mol))
    except Exception:
        return None


def _brics_fragments(mol: Chem.Mol) -> list[str]:
    """Return canonical BRICS fragments of the molecule. The count is a
    proxy for how many retrosynthetic disconnections the molecule admits."""
    frags = BRICS.BRICSDecompose(mol, returnMols=False)
    return list(frags) if frags else []


def heuristic_retro(
    smiles: str, target_yaml: Path,
) -> RetroResult:
    target = yaml.safe_load(Path(target_yaml).read_text())
    synth_cfg = target.get("synthesizability", {})
    max_sa = float(synth_cfg.get("max_sa_score", 6.0))
    max_steps = int(synth_cfg.get("retrosynthesis_max_steps", 4))

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return RetroResult(
            input_smiles=smiles,
            canonical_smiles="",
            target_id=target.get("id", ""),
            backend="heuristic",
            sa_score=None, sa_score_cap=max_sa,
            steps_estimate=None, steps_cap=max_steps,
            passed=False,
            error="RDKit could not parse SMILES",
        )

    canonical = Chem.MolToSmiles(mol)
    sa = _sa_score(mol)
    frags = _brics_fragments(mol)
    steps_est = max(0, len(frags) - 1)  # a single "fragment" means zero steps

    reasons: list[str] = []
    if sa is None:
        reasons.append("SA score unavailable (RDKit contrib missing)")
    elif sa > max_sa:
        reasons.append(f"SA score {sa:.2f} > {max_sa}")
    if steps_est > max_steps:
        reasons.append(
            f"BRICS-derived step estimate {steps_est} > max {max_steps}"
        )

    return RetroResult(
        input_smiles=smiles,
        canonical_smiles=canonical,
        target_id=target.get("id", ""),
        backend="heuristic",
        sa_score=round(sa, 3) if sa is not None else None,
        sa_score_cap=max_sa,
        steps_estimate=steps_est,
        steps_cap=max_steps,
        brics_fragments=frags,
        passed=not reasons,
        reasons=reasons,
    )


# ----------------------------------------------------------------------
# AiZynthFinder backend (optional, production)
# ----------------------------------------------------------------------


def aizynth_retro(  # pragma: no cover - requires aizynthfinder install
    smiles: str,
    target_yaml: Path,
    *,
    config_path: Path | None = None,
) -> RetroResult:
    try:
        from aizynthfinder.aizynthfinder import AiZynthFinder  # type: ignore
    except ImportError as exc:
        raise RuntimeError(
            "`pip install aizynthfinder` (plus template + expansion models) "
            "required for --aizynth backend"
        ) from exc

    target = yaml.safe_load(Path(target_yaml).read_text())
    synth_cfg = target.get("synthesizability", {})
    max_steps = int(synth_cfg.get("retrosynthesis_max_steps", 4))

    cfg = config_path or (REPO_ROOT / "data" / "aizynth" / "config.yml")
    finder = AiZynthFinder(configfile=str(cfg))
    finder.target_smiles = smiles
    finder.tree_search()
    finder.build_routes()

    routes_summary = []
    best_steps: int | None = None
    for r in finder.routes.dicts:
        n = r.get("scores", {}).get("number_of_steps")
        if n is not None and (best_steps is None or n < best_steps):
            best_steps = n
        routes_summary.append({
            "score": r.get("scores", {}).get("overall_score"),
            "steps": n,
        })

    passed = (best_steps is not None) and (best_steps <= max_steps)
    reasons = []
    if best_steps is None:
        reasons.append("aizynthfinder found no viable routes")
    elif best_steps > max_steps:
        reasons.append(f"best route needs {best_steps} steps > max {max_steps}")

    return RetroResult(
        input_smiles=smiles,
        canonical_smiles=Chem.MolToSmiles(Chem.MolFromSmiles(smiles)),
        target_id=target.get("id", ""),
        backend="aizynthfinder",
        sa_score=None, sa_score_cap=0.0,
        steps_estimate=best_steps, steps_cap=max_steps,
        aizynth_routes=routes_summary,
        passed=passed, reasons=reasons,
    )


# ----------------------------------------------------------------------
# Public API
# ----------------------------------------------------------------------


def check_retrosynthesis(
    smiles: str,
    target_yaml: Path,
    *,
    use_aizynth: bool = False,
) -> RetroResult:
    if use_aizynth:
        return aizynth_retro(smiles, target_yaml)
    return heuristic_retro(smiles, target_yaml)


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------


def _main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--target", type=Path, required=True)
    ap.add_argument("--smiles", required=True)
    ap.add_argument("--aizynth", action="store_true")
    args = ap.parse_args()

    try:
        result = check_retrosynthesis(
            args.smiles, args.target, use_aizynth=args.aizynth,
        )
    except Exception as exc:
        print(f"setup error: {exc}", file=sys.stderr)
        return 2

    print(json.dumps(asdict(result), indent=2))
    if result.error:
        return 2
    return 0 if result.passed else 1


if __name__ == "__main__":
    raise SystemExit(_main())
