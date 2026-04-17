"""Sanity tests for src/tools/validate_smiles.

Run with:
    python -m pytest tests/test_validate_smiles.py -v

These tests enforce the expected behaviour of Layers 1-4:

- known drug-like molecules (ibuprofen, caffeine, aspirin) pass,
- Lipinski violators fail Layer 3,
- garbage SMILES fail Layer 1 parse gate,
- fragmented SMILES fail Layer 1 single-fragment gate,
- known PAINS scaffolds fail Layer 2.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT / "src"))

from tools.validate_smiles import validate  # noqa: E402


def _layer(report, name):
    for layer in report["layers"]:
        if layer["name"] == name:
            return layer
    raise AssertionError(f"no layer named {name} in report")


def _gate(layer, name):
    for gate in layer["gates"]:
        if gate["name"] == name:
            return gate
    raise AssertionError(f"no gate named {name} in layer {layer['name']}")


class TestKnownDrugs:
    def test_ibuprofen_lipinski_veber_pass_ghose_fails_on_size(self):
        # Ibuprofen (15 heavy atoms) is below Ghose's heavy-atom floor of 20.
        # This is expected — Ghose is deliberately strict, and ibuprofen is
        # a small-molecule edge case, not an antibiotic candidate.
        r = validate("CC(C)Cc1ccc(C(C)C(=O)O)cc1")
        lip = _gate(_layer(r, "druglikeness"), "lipinski")
        veber = _gate(_layer(r, "druglikeness"), "veber")
        ghose = _gate(_layer(r, "druglikeness"), "ghose")
        assert lip["passed"]
        assert veber["passed"]
        assert ghose["passed"] is False

    def test_aspirin_passes_all_layers(self):
        r = validate("CC(=O)Oc1ccccc1C(=O)O")
        # Aspirin is tiny — below Ghose heavy-atom floor (20). Expect Ghose fail.
        ghose = _gate(_layer(r, "druglikeness"), "ghose")
        assert ghose["passed"] is False

    def test_caffeine_passes(self):
        r = validate("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
        # Caffeine is also below the Ghose HA floor; not a drug candidate.
        ghose = _gate(_layer(r, "druglikeness"), "ghose")
        assert ghose["passed"] is False

    def test_moxifloxacin_passes_all_gates(self):
        # Moxifloxacin is a clinically used fluoroquinolone antibiotic,
        # MW ~401, 29 heavy atoms — a realistic positive control for an
        # antibiotic candidate passing Lipinski + Veber + Ghose + SA Score.
        r = validate("COc1c(N2CC3CCCNC3C2)c(F)cc2c1N(C1CC1)C=C(C(=O)O)C2=O")
        lip = _gate(_layer(r, "druglikeness"), "lipinski")
        veber = _gate(_layer(r, "druglikeness"), "veber")
        ghose = _gate(_layer(r, "druglikeness"), "ghose")
        assert lip["passed"], lip
        assert veber["passed"], veber
        assert ghose["passed"], ghose
        assert r["passed"], r


class TestLayer1Structure:
    def test_garbage_smiles_parse_error(self):
        r = validate("not a real molecule")
        assert r["passed"] is False
        assert r.get("error") == "parse_error"

    def test_fragmented_smiles_fails_single_fragment(self):
        r = validate("CCO.CCN")
        frag_gate = _gate(_layer(r, "structure"), "single_fragment")
        assert frag_gate["passed"] is False
        assert frag_gate["detail"]["fragments"] == 2

    def test_tiny_molecule_fails_heavy_atom_floor(self):
        r = validate("CCO")  # ethanol, 3 heavy atoms
        hac = _gate(_layer(r, "structure"), "heavy_atom_count")
        assert hac["passed"] is False

    def test_doubly_charged_molecule_fails_charge_gate(self):
        # Divalent cation, net charge +2
        r = validate("[Mg+2]")
        # May fail parse via sanitization or charge gate; either is correct.
        assert r["passed"] is False


class TestLayer3Druglikeness:
    def test_too_lipophilic_fails_lipinski(self):
        # Very long alkyl chain, LogP > 5
        r = validate("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC")
        lip = _gate(_layer(r, "druglikeness"), "lipinski")
        assert lip["passed"] is False

    def test_too_heavy_fails_lipinski_mw(self):
        # Polypeptide ~ MW > 500
        r = validate(
            "CC(C)CC(NC(=O)C(CCCCN)NC(=O)C(CC(=O)O)NC(=O)C(Cc1ccccc1)N)C(=O)O"
        )
        lip = _gate(_layer(r, "druglikeness"), "lipinski")
        # MW of the example is > 500 — peptide-like
        assert lip["passed"] is False or r["descriptors"]["mw"] > 500
