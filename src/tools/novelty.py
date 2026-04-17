#!/usr/bin/env python3
"""
Layer 9: novelty check against ChEMBL + scaffold-class blocklist.

A candidate passes Layer 9 only if:
  1. The maximum Tanimoto similarity vs ChEMBL (Morgan fingerprint,
     radius 2, 2048 bits) is strictly below the target YAML threshold
     (default 0.40).
  2. The candidate does NOT contain any SMARTS scaffold listed in the
     target YAML's novelty.scaffold_class_blocklist — this blocks
     re-derivatives of known antibiotic cores (β-lactams, aminocoumarins,
     fluoroquinolones, tetracyclines, aminoglycosides).

ChEMBL similarity is queried via the public REST API
(`https://www.ebi.ac.uk/chembl/api/data/similarity/{smiles}/{pct_threshold}.json`).
Results are cached under `data/chembl-cache/` keyed by SMILES SHA256 so
repeat checks do not re-query.

CLI:
  python -m tools.novelty --target targets/mrsa-gyrb.yaml --smiles "..."
  python -m tools.novelty --target targets/mrsa-gyrb.yaml --smiles "..." --offline

Exit codes: 0 pass, 1 fail, 2 setup error.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
import urllib.parse
import urllib.request
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import yaml

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

REPO_ROOT = Path(__file__).parent.parent.parent
CACHE_DIR = REPO_ROOT / "data" / "chembl-cache"
CACHE_DIR.mkdir(parents=True, exist_ok=True)

CHEMBL_SIMILARITY_URL = (
    "https://www.ebi.ac.uk/chembl/api/data/similarity/{smiles}/{pct}.json"
    "?limit={limit}&format=json"
)


# ----------------------------------------------------------------------
# Data classes
# ----------------------------------------------------------------------


@dataclass
class ChEMBLHit:
    chembl_id: str
    pref_name: str | None
    similarity: float
    smiles: str


@dataclass
class NoveltyResult:
    input_smiles: str
    canonical_smiles: str
    target_id: str
    threshold: float
    max_tanimoto: float
    nearest_chembl_id: str | None
    nearest_chembl_name: str | None
    nearest_chembl_smiles: str | None
    chembl_hits: list[dict[str, Any]] = field(default_factory=list)
    scaffold_hits: list[str] = field(default_factory=list)
    passed: bool = False
    offline: bool = False
    error: str | None = None


# ----------------------------------------------------------------------
# Fingerprint
# ----------------------------------------------------------------------


def morgan_fingerprint(mol: Chem.Mol, radius: int = 2, n_bits: int = 2048):
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=n_bits)


def tanimoto(fp_a, fp_b) -> float:
    return float(DataStructs.TanimotoSimilarity(fp_a, fp_b))


# ----------------------------------------------------------------------
# ChEMBL lookup
# ----------------------------------------------------------------------


def _cache_path(smiles: str, pct_threshold: int) -> Path:
    key = hashlib.sha256(f"{smiles}|pct={pct_threshold}".encode()).hexdigest()[:16]
    return CACHE_DIR / f"{key}.json"


def query_chembl_similarity(
    smiles: str,
    pct_threshold: int = 70,
    limit: int = 25,
    timeout: int = 30,
) -> list[ChEMBLHit]:
    """Query the ChEMBL similarity endpoint. Returns parsed hits."""
    cache_file = _cache_path(smiles, pct_threshold)
    if cache_file.exists():
        data = json.loads(cache_file.read_text())
    else:
        encoded = urllib.parse.quote(smiles, safe="")
        url = CHEMBL_SIMILARITY_URL.format(
            smiles=encoded, pct=pct_threshold, limit=limit
        )
        req = urllib.request.Request(
            url,
            headers={"User-Agent": "orb-antibiotic-scientist/0.1 (+github)"},
        )
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            data = json.loads(resp.read())
        cache_file.write_text(json.dumps(data))

    hits: list[ChEMBLHit] = []
    for m in data.get("molecules", []):
        sim_raw = m.get("similarity")
        try:
            sim = float(sim_raw) / 100.0 if sim_raw is not None else 0.0
        except (TypeError, ValueError):
            sim = 0.0
        struct = m.get("molecule_structures") or {}
        hit_smiles = struct.get("canonical_smiles") or ""
        hits.append(
            ChEMBLHit(
                chembl_id=m.get("molecule_chembl_id", ""),
                pref_name=m.get("pref_name"),
                similarity=sim,
                smiles=hit_smiles,
            )
        )
    return hits


def _recompute_similarity_via_fingerprint(
    query_mol: Chem.Mol, hits: list[ChEMBLHit]
) -> list[ChEMBLHit]:
    """ChEMBL's returned similarity uses a ChEMBL-specific fingerprint.
    Recompute with our canonical Morgan r=2 2048-bit FP so the gate is
    consistent and reproducible."""
    q_fp = morgan_fingerprint(query_mol)
    out: list[ChEMBLHit] = []
    for h in hits:
        if not h.smiles:
            out.append(h)
            continue
        m = Chem.MolFromSmiles(h.smiles)
        if m is None:
            out.append(h)
            continue
        h.similarity = tanimoto(q_fp, morgan_fingerprint(m))
        out.append(h)
    return out


# ----------------------------------------------------------------------
# Scaffold blocklist (SMARTS)
# ----------------------------------------------------------------------


def check_scaffold_blocklist(
    mol: Chem.Mol, blocklist: list[dict[str, str]]
) -> list[str]:
    """Return the names of blocklist scaffolds that match the molecule."""
    hits: list[str] = []
    for entry in blocklist:
        name = entry.get("name", "unknown")
        smarts = entry.get("smarts", "")
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            continue
        if mol.HasSubstructMatch(patt):
            hits.append(name)
    return hits


# ----------------------------------------------------------------------
# Top-level check
# ----------------------------------------------------------------------


def check_novelty(
    smiles: str, target_yaml: Path, *, offline: bool = False
) -> NoveltyResult:
    target = yaml.safe_load(Path(target_yaml).read_text())
    nv_cfg = target.get("novelty", {})
    threshold = float(nv_cfg.get("tanimoto_threshold", 0.40))
    blocklist = nv_cfg.get("scaffold_class_blocklist", [])

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return NoveltyResult(
            input_smiles=smiles,
            canonical_smiles="",
            target_id=target.get("id", ""),
            threshold=threshold,
            max_tanimoto=1.0,
            nearest_chembl_id=None,
            nearest_chembl_name=None,
            nearest_chembl_smiles=None,
            passed=False,
            offline=offline,
            error="RDKit could not parse SMILES",
        )

    canonical = Chem.MolToSmiles(mol)
    scaffold_hits = check_scaffold_blocklist(mol, blocklist)

    if offline:
        return NoveltyResult(
            input_smiles=smiles,
            canonical_smiles=canonical,
            target_id=target.get("id", ""),
            threshold=threshold,
            max_tanimoto=float("nan"),  # unknown in offline mode
            nearest_chembl_id=None,
            nearest_chembl_name=None,
            nearest_chembl_smiles=None,
            scaffold_hits=scaffold_hits,
            passed=not scaffold_hits,  # offline ⇒ scaffold-only gate
            offline=True,
        )

    try:
        hits = query_chembl_similarity(canonical)
        hits = _recompute_similarity_via_fingerprint(mol, hits)
    except Exception as exc:
        return NoveltyResult(
            input_smiles=smiles,
            canonical_smiles=canonical,
            target_id=target.get("id", ""),
            threshold=threshold,
            max_tanimoto=float("nan"),
            nearest_chembl_id=None,
            nearest_chembl_name=None,
            nearest_chembl_smiles=None,
            scaffold_hits=scaffold_hits,
            passed=False,
            offline=False,
            error=f"ChEMBL query failed: {exc}",
        )

    hits.sort(key=lambda h: h.similarity, reverse=True)
    max_tc = hits[0].similarity if hits else 0.0
    top = hits[0] if hits else None

    passed = (max_tc < threshold) and (not scaffold_hits)

    return NoveltyResult(
        input_smiles=smiles,
        canonical_smiles=canonical,
        target_id=target.get("id", ""),
        threshold=threshold,
        max_tanimoto=max_tc,
        nearest_chembl_id=top.chembl_id if top else None,
        nearest_chembl_name=top.pref_name if top else None,
        nearest_chembl_smiles=top.smiles if top else None,
        chembl_hits=[asdict(h) for h in hits[:10]],
        scaffold_hits=scaffold_hits,
        passed=passed,
    )


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------


def _main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--target", type=Path, required=True)
    ap.add_argument("--smiles", required=True)
    ap.add_argument(
        "--offline",
        action="store_true",
        help="skip ChEMBL; only run the SMARTS scaffold blocklist",
    )
    args = ap.parse_args()

    try:
        result = check_novelty(args.smiles, args.target, offline=args.offline)
    except Exception as exc:
        print(f"setup error: {exc}", file=sys.stderr)
        return 2

    print(json.dumps(asdict(result), indent=2))
    if result.error and "parse" in result.error.lower():
        return 2
    return 0 if result.passed else 1


if __name__ == "__main__":
    raise SystemExit(_main())
