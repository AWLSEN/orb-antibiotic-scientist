"""Tests for src/tools/dock_vina.

These tests exercise the parts of the docking pipeline that do NOT require
the `vina` binary or Meeko to be installed:

  - Vina stdout → mode-energy parsing (golden reference string).
  - PDB residue centroid parsing against a synthetic minimal PDB.
  - Grid box centroid computation.
  - Target YAML loading (mrsa-gyrb).

A separate integration test (skipped if binaries absent) will run an
end-to-end dock of a reference inhibitor.
"""

from __future__ import annotations

import math
import sys
import textwrap
from pathlib import Path

import pytest

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT / "src"))

from tools.dock_vina import (  # noqa: E402
    compute_grid_box,
    parse_pdb_residue_centroids,
    parse_vina_stdout,
    VinaDocker,
)


# ----------------------------------------------------------------------
# parse_vina_stdout
# ----------------------------------------------------------------------


SAMPLE_VINA_STDOUT = """
Computing Vina grid ... done.
Performing docking (random seed: 42) ...
0%   10   20   30   40   50   60   70   80   90   100%
|----|----|----|----|----|----|----|----|----|----|
***************************************************

mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
   1       -9.2      0.000      0.000
   2       -8.7      1.234      2.456
   3       -8.3      2.105      3.100
   4       -7.9      2.899      4.012

Writing output ... done.
"""


def test_parse_vina_stdout_extracts_energies_in_order():
    energies = parse_vina_stdout(SAMPLE_VINA_STDOUT)
    assert energies == [-9.2, -8.7, -8.3, -7.9]


def test_parse_vina_stdout_handles_empty():
    assert parse_vina_stdout("") == []


def test_parse_vina_stdout_ignores_non_table_lines():
    noise = "some pre-amble\nmore noise\n" + SAMPLE_VINA_STDOUT
    assert parse_vina_stdout(noise) == [-9.2, -8.7, -8.3, -7.9]


# ----------------------------------------------------------------------
# parse_pdb_residue_centroids + compute_grid_box
# ----------------------------------------------------------------------


SYNTHETIC_PDB = textwrap.dedent("""\
ATOM      1  N   ASP A  73      10.000  20.000  30.000  1.00  0.00           N
ATOM      2  CA  ASP A  73      11.000  21.000  31.000  1.00  0.00           C
ATOM      3  C   ASP A  73      12.000  22.000  32.000  1.00  0.00           C
ATOM      4  N   GLU A  50      20.000  10.000  40.000  1.00  0.00           N
ATOM      5  CA  GLU A  50      21.000  11.000  41.000  1.00  0.00           C
ATOM      6  C   GLU A  50      22.000  12.000  42.000  1.00  0.00           C
ATOM      7  N   VAL B  73      99.000  99.000  99.000  1.00  0.00           N
ATOM      8  CA  VAL B  73      99.000  99.000  99.000  1.00  0.00           C
TER
END
""")


@pytest.fixture
def synthetic_pdb(tmp_path: Path) -> Path:
    p = tmp_path / "synthetic.pdb"
    p.write_text(SYNTHETIC_PDB)
    return p


def test_parse_pdb_residue_centroids_chain_scoped(synthetic_pdb: Path):
    # Chain A only. Both ASP73 and GLU50 should appear; the chain-B VAL73
    # must be excluded even though its resnum matches.
    centroids = parse_pdb_residue_centroids(
        synthetic_pdb, chain="A", residues=[73, 50],
    )
    assert set(centroids.keys()) == {73, 50}
    # ASP73 centroid = mean of (10,20,30), (11,21,31), (12,22,32) = (11, 21, 31)
    assert centroids[73] == pytest.approx((11.0, 21.0, 31.0))
    # GLU50 centroid = (21, 11, 41)
    assert centroids[50] == pytest.approx((21.0, 11.0, 41.0))


def test_parse_pdb_residue_centroids_missing_residue_is_absent(synthetic_pdb: Path):
    centroids = parse_pdb_residue_centroids(
        synthetic_pdb, chain="A", residues=[73, 999],
    )
    assert 999 not in centroids
    assert 73 in centroids


def test_compute_grid_box_is_centroid_of_pocket(synthetic_pdb: Path):
    gb = compute_grid_box(
        synthetic_pdb, chain="A", pocket_residues=[73, 50],
        padding_a=8.0, side_length_a=22.0,
    )
    # Centroid of ASP73 (11,21,31) and GLU50 (21,11,41) = (16, 16, 36)
    assert gb.center_x == pytest.approx(16.0)
    assert gb.center_y == pytest.approx(16.0)
    assert gb.center_z == pytest.approx(36.0)
    # Cube side is fixed at side_length_a
    assert gb.size_x == pytest.approx(22.0)
    assert gb.size_y == pytest.approx(22.0)
    assert gb.size_z == pytest.approx(22.0)


def test_compute_grid_box_raises_when_no_residues_found(synthetic_pdb: Path):
    with pytest.raises(RuntimeError, match="No pocket residues"):
        compute_grid_box(
            synthetic_pdb, chain="A", pocket_residues=[999],
            padding_a=8.0, side_length_a=22.0,
        )


# ----------------------------------------------------------------------
# VinaDocker target YAML loading (no downloads)
# ----------------------------------------------------------------------


def test_vina_docker_loads_mrsa_yaml():
    # Use the real target file. We do NOT call receptor_pdb() which would
    # hit the network; this test just verifies config wiring.
    target_yaml = ROOT / "targets" / "mrsa-gyrb.yaml"
    assert target_yaml.exists(), "mrsa-gyrb.yaml should be committed"
    d = VinaDocker(target_yaml)
    assert d.target_id == "mrsa-gyrb"
    assert d.pdb_id == "4URL"
    assert d.chain == "A"
    assert 73 in d.pocket_residues          # Asp73, catalytic
    assert 46 in d.pocket_residues          # Asn46, catalytic
    assert d.threshold_kcalmol == -8.0
    assert d.exhaustiveness == 32
    assert d.num_modes == 20
    assert d.seed == 42
    assert d.grid_side_a == 22.0


def test_benchmark_sees_reference_inhibitors():
    # Without actually running vina, confirm the docker exposes the reference
    # list from the YAML for loops/positive_control.py to iterate over.
    target_yaml = ROOT / "targets" / "mrsa-gyrb.yaml"
    d = VinaDocker(target_yaml)
    refs = d.target.get("reference_inhibitors", [])
    with_smiles = [r for r in refs if r.get("smiles")]
    assert len(with_smiles) >= 3, "at least 3 references should have SMILES"
    names = {r["name"] for r in with_smiles}
    assert "novobiocin" in names
    assert "GSK299423" in names
