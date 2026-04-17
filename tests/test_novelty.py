"""Tests for src/tools/novelty (Layer 9).

Network-free (`--offline` path): exercises the scaffold SMARTS blocklist
from targets/mrsa-gyrb.yaml plus RDKit Morgan fingerprint math. The live
ChEMBL round-trip is tested manually (opt-in) and is not in CI.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT / "src"))

from tools.novelty import (  # noqa: E402
    check_novelty,
    check_scaffold_blocklist,
    morgan_fingerprint,
    tanimoto,
)

from rdkit import Chem  # noqa: E402

TARGET_YAML = ROOT / "targets" / "mrsa-gyrb.yaml"


def test_tanimoto_self_identity():
    mol = Chem.MolFromSmiles("CC(C)Cc1ccc(C(C)C(=O)O)cc1")
    fp = morgan_fingerprint(mol)
    assert tanimoto(fp, fp) == pytest.approx(1.0)


def test_tanimoto_distinct_molecules_below_threshold():
    a = morgan_fingerprint(Chem.MolFromSmiles("CCO"))
    b = morgan_fingerprint(Chem.MolFromSmiles("c1ccncc1"))
    assert tanimoto(a, b) < 0.4


def test_scaffold_blocklist_beta_lactam_detected():
    import yaml as _yaml
    cfg = _yaml.safe_load(TARGET_YAML.read_text())
    blocklist = cfg["novelty"]["scaffold_class_blocklist"]
    # Penicillin V — classic β-lactam + thiazolidine scaffold
    penicillin_v = "CC1(C)SC2C(NC(=O)COc3ccccc3)C(=O)N2C1C(=O)O"
    mol = Chem.MolFromSmiles(penicillin_v)
    hits = check_scaffold_blocklist(mol, blocklist)
    assert "beta_lactam_core" in hits


def test_scaffold_blocklist_fluoroquinolone_detected():
    import yaml as _yaml
    cfg = _yaml.safe_load(TARGET_YAML.read_text())
    blocklist = cfg["novelty"]["scaffold_class_blocklist"]
    # Ciprofloxacin
    cipro = "O=C(O)C1=CN(C2CC2)c3cc(N4CCNCC4)c(F)cc3C1=O"
    mol = Chem.MolFromSmiles(cipro)
    hits = check_scaffold_blocklist(mol, blocklist)
    assert "fluoroquinolone_core" in hits


def test_scaffold_blocklist_aminocoumarin_detected():
    import yaml as _yaml
    cfg = _yaml.safe_load(TARGET_YAML.read_text())
    blocklist = cfg["novelty"]["scaffold_class_blocklist"]
    # Novobiocin (full molecule including the 3-amido-4-hydroxycoumarin)
    novobiocin = (
        "CC1(C)OC2C(OC(=O)N)C(OC)C(Oc3ccc(C(=O)Nc4c(O)c5cc(CC=C(C)C)"
        "c(O)cc5oc4=O)cc3C)OC2O1"
    )
    mol = Chem.MolFromSmiles(novobiocin)
    hits = check_scaffold_blocklist(mol, blocklist)
    assert "aminocoumarin_3_amide" in hits


def test_scaffold_blocklist_aminoglycoside_detected():
    import yaml as _yaml
    cfg = _yaml.safe_load(TARGET_YAML.read_text())
    blocklist = cfg["novelty"]["scaffold_class_blocklist"]
    # Kanamycin
    kanamycin = "NC1CC(N)C(OC2OC(CO)C(O)C(N)C2O)C(OC2OC(CN)C(O)C(O)C2O)C1O"
    mol = Chem.MolFromSmiles(kanamycin)
    hits = check_scaffold_blocklist(mol, blocklist)
    assert "aminoglycoside_aminosugar" in hits


def test_scaffold_blocklist_non_antibiotic_clean():
    import yaml as _yaml
    cfg = _yaml.safe_load(TARGET_YAML.read_text())
    blocklist = cfg["novelty"]["scaffold_class_blocklist"]
    # Ibuprofen — aryl propionic acid, not an antibiotic scaffold
    mol = Chem.MolFromSmiles("CC(C)Cc1ccc(C(C)C(=O)O)cc1")
    hits = check_scaffold_blocklist(mol, blocklist)
    assert hits == []


def test_offline_check_rejects_beta_lactam():
    # Offline check short-circuits ChEMBL; the scaffold blocklist alone decides.
    penicillin_v = "CC1(C)SC2C(NC(=O)COc3ccccc3)C(=O)N2C1C(=O)O"
    r = check_novelty(penicillin_v, TARGET_YAML, offline=True)
    assert r.offline is True
    assert r.passed is False
    assert "beta_lactam_core" in r.scaffold_hits


def test_offline_check_accepts_clean_novel_scaffold():
    # A random imidazopyridine — not on our blocklist and not a known class.
    candidate = "Cc1cc2ncn(C3CCN(CC3)C(=O)c3ccco3)c2cc1"
    r = check_novelty(candidate, TARGET_YAML, offline=True)
    assert r.offline is True
    assert r.scaffold_hits == []
    assert r.passed is True


def test_offline_parse_error_returns_failure():
    r = check_novelty("not a molecule", TARGET_YAML, offline=True)
    assert r.passed is False
    assert r.error is not None
