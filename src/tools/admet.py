#!/usr/bin/env python3
"""
Layer 8: ADMET / toxicity screen.

The goal of this layer is to catch candidates that are computationally plausible
binders (passed Layers 1-7) but carry red flags for absorption, distribution,
metabolism, excretion, or acute toxicity.

Strategy:
  - Local (offline) signals via RDKit: QED, TPSA, LogP-derived BBB + GI proxies,
    structural alerts via Brenk / NIH / ZINC FilterCatalogs, common toxicophore
    SMARTS that aren't covered by PAINS/REOS.
  - Optional online signals (live APIs): SwissADME and ProTox-II. These are
    `--online` only; the default offline mode is the CI-safe path.
  - Gates are loaded from `admet_gates` in the target YAML. Any violation is a
    red flag; the layer passes only when no red flags are raised.

CLI:
  python -m tools.admet --target targets/mrsa-gyrb.yaml --smiles "..."
  python -m tools.admet --target targets/mrsa-gyrb.yaml --smiles "..." --online

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
from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, QED, rdMolDescriptors
from rdkit.Chem import FilterCatalog

RDLogger.DisableLog("rdApp.*")

REPO_ROOT = Path(__file__).parent.parent.parent


# Toxicophores beyond PAINS/REOS. These SMARTS flag substructures that are
# widely accepted as red flags in medicinal chemistry but are not in RDKit's
# FilterCatalog PAINS A/B/C set.
TOXICOPHORE_SMARTS: dict[str, str] = {
    # Polyhalogenated aromatics: bioaccumulation, hepatotoxicity risk.
    "polyhalogenated_arene": "c1(:c:c:c:c:c:1)[F,Cl,Br,I]" ,
    # Nitro aromatic — bacterial mutagen (Ames-positive) risk.
    "nitro_aromatic": "c[N+](=O)[O-]",
    # Isocyanate / isothiocyanate — highly reactive.
    "isocyanate": "N=C=O",
    "isothiocyanate": "N=C=S",
    # Aliphatic aldehyde outside of glyco/pyranose context — electrophile.
    "aliphatic_aldehyde": "[CX3H1](=O)[#6;!a]",
    # Hydrazine / hydrazone — Ames-positive hotspot.
    "hydrazine": "[NH1,NH2][NH1,NH2]",
    "hydrazone": "C=N[NH1,NH2]",
    # Aromatic amine → potential Ames-positive metabolite.
    "aniline_primary": "c[NX3;H2]",
    # Michael acceptor (α,β-unsaturated carbonyl) in open-chain context.
    "michael_acceptor_openchain": "[CX3]=[CX3][CX3](=O)[!a]",
    # Acyl halide / sulfonyl halide — hyper-reactive.
    "acyl_halide": "[CX3](=O)[Cl,Br,I,F]",
    "sulfonyl_halide": "[SX4](=O)(=O)[Cl,Br,I,F]",
    # Methylene-dioxyphenyl — CYP3A4 mechanism-based inhibitor risk.
    "methylenedioxyphenyl": "c1ccc2OCOc2c1",
    # Quaternary ammonium salt — poor oral bioavailability.
    "quaternary_nitrogen": "[N+;X4]",
}


@dataclass
class AdmetResult:
    input_smiles: str
    canonical_smiles: str
    target_id: str
    qed: float | None
    mw: float | None
    logp: float | None
    tpsa: float | None
    hbd: int | None
    hba: int | None
    rotatable_bonds: int | None
    bbb_likely: bool | None
    oral_bioavailability_proxy: str | None
    structural_alerts: list[str] = field(default_factory=list)
    toxicophore_hits: list[str] = field(default_factory=list)
    red_flags: list[str] = field(default_factory=list)
    passed: bool = False
    online_used: bool = False
    online_results: dict[str, Any] = field(default_factory=dict)
    error: str | None = None


# ----------------------------------------------------------------------
# RDKit structural-alert catalogs (Brenk, NIH, ZINC)
# ----------------------------------------------------------------------


_ALERT_CATALOG: FilterCatalog.FilterCatalog | None = None


def _alert_catalog() -> FilterCatalog.FilterCatalog:
    global _ALERT_CATALOG
    if _ALERT_CATALOG is None:
        params = FilterCatalog.FilterCatalogParams()
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.NIH)
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.ZINC)
        _ALERT_CATALOG = FilterCatalog.FilterCatalog(params)
    return _ALERT_CATALOG


def check_structural_alerts(mol: Chem.Mol) -> list[str]:
    matches = _alert_catalog().GetMatches(mol)
    return [m.GetDescription() for m in matches]


def check_toxicophores(mol: Chem.Mol) -> list[str]:
    hits: list[str] = []
    for label, smarts in TOXICOPHORE_SMARTS.items():
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            continue
        if mol.HasSubstructMatch(patt):
            hits.append(label)
    return hits


# ----------------------------------------------------------------------
# Derived proxies
# ----------------------------------------------------------------------


def bbb_proxy(mw: float, logp: float, tpsa: float, hbd: int) -> bool:
    """Egan + variant heuristic for blood-brain barrier permeability.

    Returns True iff the molecule is likely to cross the BBB. For an
    anti-bacterial targeting a peripheral infection we generally *don't*
    want BBB permeation unless the target is CNS-resident, so this is a
    flag-raiser rather than a hard gate.
    """
    if mw >= 500:
        return False
    if logp > 5.0 or logp < -1.0:
        return False
    if tpsa > 90:
        return False
    if hbd > 5:
        return False
    return True


def oral_bioavailability_proxy(
    mw: float, logp: float, tpsa: float, rotb: int
) -> str:
    """Coarse Veber-style F(oral) bucket: high / medium / low."""
    if mw <= 500 and logp <= 5 and tpsa <= 140 and rotb <= 10:
        return "high"
    if mw <= 600 and tpsa <= 160:
        return "medium"
    return "low"


# ----------------------------------------------------------------------
# Top-level check
# ----------------------------------------------------------------------


def check_admet(
    smiles: str,
    target_yaml: Path,
    *,
    online: bool = False,
) -> AdmetResult:
    target = yaml.safe_load(Path(target_yaml).read_text())
    gates = target.get("admet_gates", {})

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return AdmetResult(
            input_smiles=smiles,
            canonical_smiles="",
            target_id=target.get("id", ""),
            qed=None, mw=None, logp=None, tpsa=None,
            hbd=None, hba=None, rotatable_bonds=None,
            bbb_likely=None, oral_bioavailability_proxy=None,
            passed=False,
            error="RDKit could not parse SMILES",
        )

    canonical = Chem.MolToSmiles(mol)

    qed = float(QED.qed(mol))
    mw = float(Descriptors.MolWt(mol))
    logp = float(Crippen.MolLogP(mol))
    tpsa = float(rdMolDescriptors.CalcTPSA(mol))
    hbd = int(Lipinski.NumHDonors(mol))
    hba = int(Lipinski.NumHAcceptors(mol))
    rotb = int(Lipinski.NumRotatableBonds(mol))

    bbb = bbb_proxy(mw, logp, tpsa, hbd)
    oral = oral_bioavailability_proxy(mw, logp, tpsa, rotb)

    alerts = check_structural_alerts(mol)
    toxies = check_toxicophores(mol)

    red_flags: list[str] = []

    # YAML gates
    lip = gates.get("lipinski", {})
    if lip:
        if mw > lip.get("max_mw", 500):
            red_flags.append(f"lipinski.mw({mw:.1f}>{lip['max_mw']})")
        if logp > lip.get("max_logp", 5.0):
            red_flags.append(f"lipinski.logp({logp:.2f}>{lip['max_logp']})")
        if hbd > lip.get("max_hbd", 5):
            red_flags.append(f"lipinski.hbd({hbd}>{lip['max_hbd']})")
        if hba > lip.get("max_hba", 10):
            red_flags.append(f"lipinski.hba({hba}>{lip['max_hba']})")
    veber = gates.get("veber", {})
    if veber:
        if tpsa > veber.get("max_tpsa", 140):
            red_flags.append(f"veber.tpsa({tpsa:.1f}>{veber['max_tpsa']})")
        if rotb > veber.get("max_rotatable_bonds", 10):
            red_flags.append(
                f"veber.rotb({rotb}>{veber['max_rotatable_bonds']})"
            )
    pk = gates.get("pharmacokinetics", {})
    min_bioavail = pk.get("oral_bioavailability_pct_min")
    if min_bioavail is not None and oral == "low":
        red_flags.append(
            f"pk.oral_bioavailability({oral} < needed {min_bioavail}%)"
        )

    # Structural-alert catalog (Brenk/NIH/ZINC)
    if alerts:
        red_flags.append(f"structural_alerts:{len(alerts)}")
    # Extra toxicophores
    if toxies:
        red_flags.append("toxicophore:" + ",".join(toxies))

    online_used = False
    online_results: dict[str, Any] = {}
    if online:
        try:
            online_results = _online_admet(canonical)
            online_used = True
            if online_results.get("ames_positive") is True:
                red_flags.append("online.ames_positive")
            if online_results.get("herg_positive") is True:
                red_flags.append("online.herg_positive")
            if online_results.get("hepatotox_positive") is True:
                red_flags.append("online.hepatotox_positive")
        except Exception as exc:
            online_results = {"error": str(exc)}

    passed = not red_flags

    return AdmetResult(
        input_smiles=smiles,
        canonical_smiles=canonical,
        target_id=target.get("id", ""),
        qed=round(qed, 3),
        mw=round(mw, 2),
        logp=round(logp, 3),
        tpsa=round(tpsa, 2),
        hbd=hbd,
        hba=hba,
        rotatable_bonds=rotb,
        bbb_likely=bbb,
        oral_bioavailability_proxy=oral,
        structural_alerts=alerts,
        toxicophore_hits=toxies,
        red_flags=red_flags,
        passed=passed,
        online_used=online_used,
        online_results=online_results,
    )


# ----------------------------------------------------------------------
# Online ADMET (stub — wires SwissADME / ProTox-II / pkCSM)
# ----------------------------------------------------------------------


def _online_admet(smiles: str) -> dict[str, Any]:  # pragma: no cover - IO
    """Placeholder for online ADMET predictions. The production agent will
    wire SwissADME (HTML form post) or an ADMET-AI hosted model here. We
    keep this abstract to avoid flaky network dependencies in CI; the
    expected return shape is documented below."""
    return {
        "note": "online ADMET not implemented; wire SwissADME/ProTox/pkCSM here",
        "ames_positive": None,
        "herg_positive": None,
        "hepatotox_positive": None,
    }


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------


def _main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--target", type=Path, required=True)
    ap.add_argument("--smiles", required=True)
    ap.add_argument("--online", action="store_true", help="use online ADMET APIs")
    args = ap.parse_args()

    try:
        r = check_admet(args.smiles, args.target, online=args.online)
    except Exception as exc:
        print(f"setup error: {exc}", file=sys.stderr)
        return 2

    print(json.dumps(asdict(r), indent=2))
    if r.error:
        return 2
    return 0 if r.passed else 1


if __name__ == "__main__":
    raise SystemExit(_main())
