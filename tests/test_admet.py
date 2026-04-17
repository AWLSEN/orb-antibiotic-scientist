"""Tests for src/tools/admet (Layer 8)."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT / "src"))

from tools.admet import (  # noqa: E402
    bbb_proxy,
    check_admet,
    check_structural_alerts,
    check_toxicophores,
    oral_bioavailability_proxy,
)

from rdkit import Chem  # noqa: E402

TARGET_YAML = ROOT / "targets" / "mrsa-gyrb.yaml"


def test_bbb_proxy_amphetamine_likely_bbb():
    # Amphetamine: MW 135, LogP ~1.8, TPSA 26 — clearly BBB-penetrant.
    assert bbb_proxy(mw=135, logp=1.8, tpsa=26, hbd=1) is True


def test_bbb_proxy_big_polar_unlikely_bbb():
    # Polar macromolecule-ish — TPSA too high.
    assert bbb_proxy(mw=480, logp=3.0, tpsa=120, hbd=3) is False


def test_oral_bioavailability_proxy_tiers():
    assert oral_bioavailability_proxy(mw=250, logp=2.0, tpsa=60, rotb=3) == "high"
    assert oral_bioavailability_proxy(mw=550, logp=3.0, tpsa=150, rotb=8) == "medium"
    assert oral_bioavailability_proxy(mw=700, logp=7.0, tpsa=200, rotb=15) == "low"


def test_toxicophore_nitro_aromatic_flagged():
    mol = Chem.MolFromSmiles("[O-][N+](=O)c1ccccc1")
    hits = check_toxicophores(mol)
    assert "nitro_aromatic" in hits


def test_toxicophore_isocyanate_flagged():
    mol = Chem.MolFromSmiles("O=C=NCCCC")
    hits = check_toxicophores(mol)
    assert "isocyanate" in hits


def test_toxicophore_clean_molecule_no_hits():
    mol = Chem.MolFromSmiles("CC(C)Cc1ccc(C(C)C(=O)O)cc1")  # ibuprofen
    hits = check_toxicophores(mol)
    # Ibuprofen is clean of our manual toxicophore list.
    assert hits == []


def test_structural_alerts_brenk_nitro_aromatic():
    # Brenk catalog flags nitro_aromatic; we expect at least 1 alert.
    mol = Chem.MolFromSmiles("[O-][N+](=O)c1ccccc1")
    alerts = check_structural_alerts(mol)
    assert len(alerts) >= 1


def test_check_admet_clean_candidate_passes():
    # Moxifloxacin — clinical antibiotic, should pass the Layer 8 gates
    # (it does trigger the fluoroquinolone scaffold alert in Layer 9 but
    # that is a separate layer; Layer 8 cares about ADMET, not novelty).
    smi = "COc1c(N2CC3CCCNC3C2)c(F)cc2c1N(C1CC1)C=C(C(=O)O)C2=O"
    r = check_admet(smi, TARGET_YAML, online=False)
    # Even clinical antibiotics sometimes hit Brenk alerts (heterocycles can
    # pattern-match structural-alert SMARTS). What we assert: the YAML
    # ADMET gates (Lipinski, Veber, PK) should not flag moxifloxacin.
    yaml_flags = [f for f in r.red_flags if f.startswith(("lipinski", "veber", "pk."))]
    assert yaml_flags == [], yaml_flags
    # And the basic descriptors are populated.
    assert r.qed is not None and 0.0 <= r.qed <= 1.0
    assert r.mw is not None and r.mw > 300


def test_check_admet_toxic_candidate_flagged():
    # Nitroaromatic alcohol — Brenk + toxicophore hit.
    smi = "[O-][N+](=O)c1ccc(CO)cc1"
    r = check_admet(smi, TARGET_YAML, online=False)
    assert r.passed is False
    assert any("nitro_aromatic" in f or "structural_alerts" in f for f in r.red_flags)


def test_check_admet_parse_error_returns_failure():
    r = check_admet("not a molecule", TARGET_YAML, online=False)
    assert r.passed is False
    assert r.error is not None


def test_check_admet_oversize_candidate_flagged_by_yaml_gates():
    # 11 MW-inflating carbons, LogP > 5.
    smi = "CCCCCCCCCCCCCCCCCCCCCCC(=O)O"
    r = check_admet(smi, TARGET_YAML, online=False)
    assert r.passed is False
    assert any(f.startswith("lipinski") for f in r.red_flags)
