#!/usr/bin/env python3
"""
Layers 1-4 of the verification chain: structural sanity for a candidate SMILES.

Layer 1: RDKit parse + valence / charge sanity.
Layer 2: PAINS (pan-assay interference compounds) + REOS filters.
Layer 3: Drug-likeness — Lipinski, Veber, Ghose.
Layer 4: Synthetic accessibility (Ertl SA Score, 1 = easy, 10 = hard).

Usage (library):
    from tools.validate_smiles import validate
    report = validate("CCOc1ccc2nc(S(=O)Cc3ncc(C)c(OC)c3C)[nH]c2c1")
    report["passed"]  # bool
    report["layers"]  # dict per layer with gate outcomes

Usage (CLI):
    python -m tools.validate_smiles "SMILES_STRING"
    python -m tools.validate_smiles --batch file.smi    # one SMILES per line

Exit codes: 0 = passed, 1 = failed, 2 = parse error.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, rdMolDescriptors
    from rdkit.Chem import FilterCatalog
    from rdkit.Chem import RDConfig
except ImportError as e:  # pragma: no cover - runtime environment issue
    sys.stderr.write(
        "[validate_smiles] RDKit not installed. `pip install rdkit` "
        "(or `conda install -c conda-forge rdkit`).\n"
    )
    raise

# Silence RDKit chatter when parsing many SMILES.
RDLogger.DisableLog("rdApp.*")

# Ertl SA Score ships as an RDKit contrib module.
_SA_SCORE_AVAILABLE = False
try:
    _contrib_path = os.path.join(RDConfig.RDContribDir, "SA_Score")
    if _contrib_path not in sys.path:
        sys.path.append(_contrib_path)
    import sascorer  # type: ignore[import-not-found]
    _SA_SCORE_AVAILABLE = True
except Exception:
    sascorer = None  # type: ignore[assignment]


# ----------------------------------------------------------------------
# Gates
# ----------------------------------------------------------------------


@dataclass
class GateResult:
    name: str
    passed: bool
    detail: dict[str, Any] = field(default_factory=dict)
    reason: str = ""


@dataclass
class LayerResult:
    layer: int
    name: str
    passed: bool
    gates: list[GateResult] = field(default_factory=list)


# REOS (rapid elimination of swill) — a small curated set of unwanted groups.
# See Walters et al., Advanced Drug Delivery Reviews, 2002.
_REOS_SMARTS = {
    "reactive_aldehyde": "[CX3H1](=O)[#6]",
    "reactive_michael_acceptor": "[CH]=[CH][C]=O",
    "nitro_aromatic": "c[N+](=O)[O-]",
    "thiol": "[SH]",
    "disulfide": "[S][S]",
    "polyhalide": "[Cl,Br,I][Cl,Br,I][Cl,Br,I]",
    "azo": "N=N",
    "epoxide_unactivated": "C1OC1",
    "aziridine": "C1NC1",
    "peroxide": "OO",
    "quaternary_nitrogen": "[N+](=O)(=O)",
}


def _parse(smiles: str) -> Chem.Mol | None:
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    try:
        Chem.SanitizeMol(m)
    except Exception:
        return None
    return m


# ----------------------------------------------------------------------
# Layer 1: parse + valence / charge
# ----------------------------------------------------------------------


def _layer1_structure(smiles: str, mol: Chem.Mol | None) -> LayerResult:
    gates: list[GateResult] = []
    passed = True

    if mol is None:
        gates.append(
            GateResult("parse", False, reason="RDKit could not parse the SMILES"),
        )
        return LayerResult(1, "structure", False, gates)

    gates.append(GateResult("parse", True, {"canonical_smiles": Chem.MolToSmiles(mol)}))

    # Heavy atom count sanity
    hac = mol.GetNumHeavyAtoms()
    gates.append(
        GateResult(
            "heavy_atom_count",
            10 <= hac <= 70,
            {"heavy_atoms": hac},
            ""
            if 10 <= hac <= 70
            else f"heavy atoms={hac}; outside 10-70 drug-like range",
        )
    )
    if not (10 <= hac <= 70):
        passed = False

    # Net formal charge (prefer neutral or single-charge)
    charge = Chem.rdmolops.GetFormalCharge(mol)
    charge_ok = abs(charge) <= 1
    gates.append(
        GateResult(
            "net_formal_charge",
            charge_ok,
            {"charge": charge},
            "" if charge_ok else f"net formal charge {charge} outside [-1, +1]",
        )
    )
    if not charge_ok:
        passed = False

    # Unusual valence via sanitize fallback — already covered by _parse.
    # Fragment count: reject if not a single connected molecule.
    frags = Chem.GetMolFrags(mol, asMols=False)
    single_frag = len(frags) == 1
    gates.append(
        GateResult(
            "single_fragment",
            single_frag,
            {"fragments": len(frags)},
            "" if single_frag else f"{len(frags)} disconnected fragments",
        )
    )
    if not single_frag:
        passed = False

    if not all(g.passed for g in gates):
        passed = False

    return LayerResult(1, "structure", passed, gates)


# ----------------------------------------------------------------------
# Layer 2: PAINS + REOS
# ----------------------------------------------------------------------


_PAINS_CATALOG: FilterCatalog.FilterCatalog | None = None


def _pains_catalog() -> FilterCatalog.FilterCatalog:
    global _PAINS_CATALOG
    if _PAINS_CATALOG is None:
        params = FilterCatalog.FilterCatalogParams()
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A)
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B)
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C)
        _PAINS_CATALOG = FilterCatalog.FilterCatalog(params)
    return _PAINS_CATALOG


def _layer2_hygiene(mol: Chem.Mol) -> LayerResult:
    gates: list[GateResult] = []

    # PAINS
    cat = _pains_catalog()
    pains_matches = cat.GetMatches(mol)
    pains_hit = len(pains_matches) > 0
    pains_names = [m.GetDescription() for m in pains_matches]
    gates.append(
        GateResult(
            "pains",
            not pains_hit,
            {"matches": pains_names},
            "" if not pains_hit else f"PAINS hits: {pains_names}",
        )
    )

    # REOS
    reos_hits: list[str] = []
    for label, smarts in _REOS_SMARTS.items():
        patt = Chem.MolFromSmarts(smarts)
        if patt is not None and mol.HasSubstructMatch(patt):
            reos_hits.append(label)
    gates.append(
        GateResult(
            "reos",
            not reos_hits,
            {"matches": reos_hits},
            "" if not reos_hits else f"REOS hits: {reos_hits}",
        )
    )

    passed = all(g.passed for g in gates)
    return LayerResult(2, "hygiene", passed, gates)


# ----------------------------------------------------------------------
# Layer 3: drug-likeness — Lipinski + Veber + Ghose
# ----------------------------------------------------------------------


def _descriptor_bundle(mol: Chem.Mol) -> dict[str, float]:
    return {
        "mw": Descriptors.MolWt(mol),
        "logp": Crippen.MolLogP(mol),
        "hbd": Lipinski.NumHDonors(mol),
        "hba": Lipinski.NumHAcceptors(mol),
        "tpsa": rdMolDescriptors.CalcTPSA(mol),
        "rotatable_bonds": Lipinski.NumRotatableBonds(mol),
        "heavy_atoms": mol.GetNumHeavyAtoms(),
        "molar_refractivity": Crippen.MolMR(mol),
    }


def _layer3_druglikeness(mol: Chem.Mol) -> LayerResult:
    d = _descriptor_bundle(mol)
    gates: list[GateResult] = []

    # Lipinski: MW ≤ 500, LogP ≤ 5, HBD ≤ 5, HBA ≤ 10 (up to 1 violation in
    # classical Pfizer Ro5 — here we require strict adherence).
    lip_violations: list[str] = []
    if d["mw"] > 500:
        lip_violations.append(f"MW {d['mw']:.1f} > 500")
    if d["logp"] > 5:
        lip_violations.append(f"LogP {d['logp']:.2f} > 5")
    if d["hbd"] > 5:
        lip_violations.append(f"HBD {d['hbd']} > 5")
    if d["hba"] > 10:
        lip_violations.append(f"HBA {d['hba']} > 10")
    gates.append(
        GateResult(
            "lipinski",
            not lip_violations,
            {k: d[k] for k in ("mw", "logp", "hbd", "hba")},
            "" if not lip_violations else "; ".join(lip_violations),
        )
    )

    # Veber: TPSA ≤ 140, rotatable bonds ≤ 10
    veber_violations: list[str] = []
    if d["tpsa"] > 140:
        veber_violations.append(f"TPSA {d['tpsa']:.1f} > 140")
    if d["rotatable_bonds"] > 10:
        veber_violations.append(f"rot bonds {d['rotatable_bonds']} > 10")
    gates.append(
        GateResult(
            "veber",
            not veber_violations,
            {k: d[k] for k in ("tpsa", "rotatable_bonds")},
            "" if not veber_violations else "; ".join(veber_violations),
        )
    )

    # Ghose: MW 160-480, LogP -0.4-5.6, MR 40-130, atoms 20-70
    ghose_violations: list[str] = []
    if not (160 <= d["mw"] <= 480):
        ghose_violations.append(f"MW {d['mw']:.1f} not in [160, 480]")
    if not (-0.4 <= d["logp"] <= 5.6):
        ghose_violations.append(f"LogP {d['logp']:.2f} not in [-0.4, 5.6]")
    if not (40 <= d["molar_refractivity"] <= 130):
        ghose_violations.append(
            f"MR {d['molar_refractivity']:.1f} not in [40, 130]"
        )
    if not (20 <= d["heavy_atoms"] <= 70):
        ghose_violations.append(
            f"heavy atoms {d['heavy_atoms']} not in [20, 70]"
        )
    gates.append(
        GateResult(
            "ghose",
            not ghose_violations,
            {k: d[k] for k in ("mw", "logp", "molar_refractivity", "heavy_atoms")},
            "" if not ghose_violations else "; ".join(ghose_violations),
        )
    )

    passed = all(g.passed for g in gates)
    return LayerResult(3, "druglikeness", passed, gates)


# ----------------------------------------------------------------------
# Layer 4: Synthetic accessibility (Ertl SA Score)
# ----------------------------------------------------------------------


def _layer4_sa_score(mol: Chem.Mol, max_sa: float = 6.0) -> LayerResult:
    gates: list[GateResult] = []
    if not _SA_SCORE_AVAILABLE:
        gates.append(
            GateResult(
                "sa_score_available",
                False,
                reason=(
                    "RDKit SA_Score contrib not available on this system. "
                    "`pip install rdkit` includes it by default; verify "
                    "RDConfig.RDContribDir/SA_Score/sascorer.py exists."
                ),
            )
        )
        return LayerResult(4, "synthesizability", False, gates)

    sa = float(sascorer.calculateScore(mol))
    passed = sa <= max_sa
    gates.append(
        GateResult(
            "sa_score",
            passed,
            {"sa_score": round(sa, 3), "max_sa": max_sa},
            "" if passed else f"SA {sa:.2f} > {max_sa} (too hard to synthesize)",
        )
    )
    return LayerResult(4, "synthesizability", passed, gates)


# ----------------------------------------------------------------------
# Public API
# ----------------------------------------------------------------------


def validate(smiles: str, *, max_sa: float = 6.0) -> dict[str, Any]:
    """Run Layers 1-4 on a single SMILES. Returns a structured report."""
    mol = _parse(smiles)

    layers: list[LayerResult] = [_layer1_structure(smiles, mol)]

    # If Layer 1 fails at the parse gate we cannot proceed.
    if mol is None:
        return {
            "input": smiles,
            "passed": False,
            "error": "parse_error",
            "layers": [asdict(layers[0])],
        }

    layers.append(_layer2_hygiene(mol))
    layers.append(_layer3_druglikeness(mol))
    layers.append(_layer4_sa_score(mol, max_sa=max_sa))

    all_passed = all(layer.passed for layer in layers)

    return {
        "input": smiles,
        "canonical_smiles": Chem.MolToSmiles(mol),
        "passed": all_passed,
        "layers": [asdict(layer) for layer in layers],
        "descriptors": _descriptor_bundle(mol),
    }


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------


def _main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("smiles", nargs="?", help="single SMILES string")
    ap.add_argument("--batch", type=Path, help="file with one SMILES per line")
    ap.add_argument(
        "--max-sa",
        type=float,
        default=6.0,
        help="SA Score threshold (default 6.0; 1=easy, 10=hard)",
    )
    ap.add_argument(
        "--json-only", action="store_true", help="emit only JSON to stdout"
    )
    args = ap.parse_args()

    if args.batch:
        results = []
        for line in args.batch.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            results.append(validate(line, max_sa=args.max_sa))
        print(json.dumps(results, indent=2 if not args.json_only else None))
        return 0 if all(r["passed"] for r in results) else 1

    if not args.smiles:
        ap.error("provide SMILES or --batch")

    report = validate(args.smiles, max_sa=args.max_sa)
    print(json.dumps(report, indent=2 if not args.json_only else None))

    if report.get("error") == "parse_error":
        return 2
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(_main())
