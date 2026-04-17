"""Tests for src/tools/mechanism (Layer 7).

The gate logic (check_gate) is a pure function over the target YAML and a
residue set, so we test it exhaustively without needing a docking pose.
The distance backend is exercised against a synthetic mini-PDB.
"""

from __future__ import annotations

import sys
import textwrap
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT / "src"))

from tools.mechanism import (  # noqa: E402
    analyze_pose,
    check_gate,
    _distance_backend_contacts,
)

TARGET_YAML = ROOT / "targets" / "mrsa-gyrb.yaml"


@pytest.fixture(scope="module")
def target():
    return yaml.safe_load(TARGET_YAML.read_text())


# ----------------------------------------------------------------------
# check_gate — pure logic
# ----------------------------------------------------------------------


def test_gate_passes_with_two_catalytic_contacts(target):
    # Asp73 (catalytic, must-include) + Asn46 (catalytic, must-include)
    gate = check_gate([73, 46], target)
    assert gate["passed"], gate
    assert 73 in gate["catalytic_residues_hit"]
    assert 46 in gate["catalytic_residues_hit"]
    assert set(gate["must_include_hits"]) == {46, 73}


def test_gate_fails_when_only_one_contact(target):
    # Single contact — below minimum of 2.
    gate = check_gate([73], target)
    assert gate["passed"] is False
    assert any("too few contacts" in r for r in gate["reasons"])


def test_gate_fails_when_no_must_include_hit(target):
    # Two pocket contacts but neither in {46, 50, 73, 76}.
    # Val43 and Ile51 are pocket residues but not catalytic / must-include.
    gate = check_gate([43, 51], target)
    assert gate["passed"] is False
    assert any("must_include_any_of" in r for r in gate["reasons"])


def test_gate_fails_when_only_resistance_hotspots(target):
    # Gly85 and Arg136 are single-site resistance hotspots in this YAML.
    # Even though there are 2 contacts, both are evolvable → reject.
    gate = check_gate([85, 136], target)
    assert gate["passed"] is False
    assert gate["resistance_only"] is True
    assert any("resistance" in r for r in gate["reasons"])


def test_gate_passes_when_resistance_plus_non_hotspot(target):
    # Hitting a resistance residue is OK *if* the ligand also hits a
    # non-hotspot catalytic residue — resistance mutations at the hotspot
    # cannot single-handedly abolish binding.
    gate = check_gate([85, 73], target)
    # 73 is a catalytic must-include residue, 85 is resistance hotspot.
    # resistance_only should be False; must_include_hits contains 73.
    assert gate["resistance_only"] is False
    assert 73 in gate["must_include_hits"]
    assert gate["passed"] is True


def test_gate_ignores_residues_outside_pocket(target):
    # Random residue 200 is not in pocket; does not contribute to hits.
    gate = check_gate([73, 200], target)
    assert 200 not in gate["pocket_residues_hit"]
    # 73 alone isn't enough for minimum=2 contacts in pocket, but
    # check_gate counts total contacts (len of set), not pocket-only.
    # 73 is catalytic must-include and contacts_set has 2 entries, so passes.
    assert gate["passed"] is True


def test_gate_empty_contact_set_fails(target):
    gate = check_gate([], target)
    assert gate["passed"] is False
    assert "too few" in " ".join(gate["reasons"])


# ----------------------------------------------------------------------
# Distance backend — synthetic PDB
# ----------------------------------------------------------------------


SYNTHETIC_COMPLEX_PDB = textwrap.dedent("""\
ATOM      1  CA  ASP A  73       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  CB  ASP A  73       1.000   0.000   0.000  1.00  0.00           C
ATOM      3  CA  ASN A  46       2.000   0.000   0.000  1.00  0.00           C
ATOM      4  CA  VAL A  43      10.000   0.000   0.000  1.00  0.00           C
ATOM      5  CA  ILE A  51      20.000   0.000   0.000  1.00  0.00           C
ATOM      6  CA  VAL B  43       0.500   0.000   0.000  1.00  0.00           C
HETATM    7  C1  LIG X   1       0.000   0.000   2.000  1.00  0.00           C
HETATM    8  C2  LIG X   1       2.000   0.000   2.000  1.00  0.00           C
END
""")


def test_distance_backend_flags_close_residues(tmp_path: Path):
    p = tmp_path / "complex.pdb"
    p.write_text(SYNTHETIC_COMPLEX_PDB)

    # Ligand atoms at (0,0,2) and (2,0,2). Cutoff 3.0 Å.
    #   ASP 73 CA at (0,0,0) → 2.0 Å from first ligand atom → HIT
    #   ASN 46 CA at (2,0,0) → 2.0 Å from second ligand atom → HIT
    #   VAL 43 CA at (10,0,0) → 10 Å → MISS
    #   ILE 51 CA at (20,0,0) → 20 Å → MISS
    #   Chain-B VAL 43 at (0.5,0,0) → should be excluded by chain filter
    contacts = _distance_backend_contacts(
        p, chain="A", pocket_residues=[73, 46, 43, 51], cutoff_a=3.0,
    )
    res_hit = {c.residue_num for c in contacts}
    assert res_hit == {73, 46}


def test_distance_backend_respects_pocket_residue_list(tmp_path: Path):
    p = tmp_path / "complex.pdb"
    p.write_text(SYNTHETIC_COMPLEX_PDB)

    # Restrict pocket to just 46 — only ASN 46 should surface even though
    # ASP 73 is also within cutoff.
    contacts = _distance_backend_contacts(
        p, chain="A", pocket_residues=[46], cutoff_a=3.0,
    )
    assert [c.residue_num for c in contacts] == [46]


# ----------------------------------------------------------------------
# analyze_pose end-to-end with synthetic pose
# ----------------------------------------------------------------------


def test_analyze_pose_uses_distance_backend_and_passes_gate(tmp_path: Path):
    p = tmp_path / "complex.pdb"
    p.write_text(SYNTHETIC_COMPLEX_PDB)

    # Force fallback backend to avoid PLIP's openbabel dependency in CI.
    result = analyze_pose(
        p, TARGET_YAML, "cand-test", prefer_plip=False, cutoff_a=3.0,
    )
    assert result.backend == "distance"
    assert 73 in result.catalytic_residues_hit
    assert 46 in result.catalytic_residues_hit
    assert result.passed is True
    assert result.num_contacts == 2
