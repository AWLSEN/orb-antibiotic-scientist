"""Tests for src/scoring."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT / "src"))

from scoring import (  # noqa: E402
    score_admet,
    score_candidate,
    score_docking,
    score_druglikeness,
    score_mechanism,
    score_novelty,
    score_synthesizability,
    update_leaderboard,
)

TARGET_YAML = ROOT / "targets" / "mrsa-gyrb.yaml"


# ----------------------------------------------------------------------
# Per-layer scoring
# ----------------------------------------------------------------------


class TestScoreDocking:
    def test_weak_dg_zero(self):
        assert score_docking({"best_energy_kcalmol": -5.0}) == pytest.approx(0.0)

    def test_strong_dg_one(self):
        assert score_docking({"best_energy_kcalmol": -11.0}) == pytest.approx(1.0)

    def test_midrange_linear(self):
        # −8.0 → (8 − 5) / 6 = 0.5
        assert score_docking({"best_energy_kcalmol": -8.0}) == pytest.approx(0.5)

    def test_super_strong_capped_at_one(self):
        assert score_docking({"best_energy_kcalmol": -15.0}) == pytest.approx(1.0)

    def test_missing_energy_zero(self):
        assert score_docking({}) == 0.0
        assert score_docking(None) == 0.0


class TestScoreNovelty:
    def test_scaffold_hit_zeroes(self):
        assert score_novelty({"scaffold_hits": ["beta_lactam"], "max_tanimoto": 0.1}) == 0.0

    def test_low_tanimoto_high_score(self):
        assert score_novelty({"scaffold_hits": [], "max_tanimoto": 0.1}) == pytest.approx(0.9)

    def test_high_tanimoto_low_score(self):
        assert score_novelty({"scaffold_hits": [], "max_tanimoto": 0.8}) == pytest.approx(0.2)

    def test_nan_gives_partial_credit(self):
        assert score_novelty({"scaffold_hits": [], "max_tanimoto": float("nan")}) == 0.5


class TestScoreDruglikeness:
    def test_all_gates_pass_score_one(self):
        report = {
            "layers": [
                {"name": "druglikeness", "gates": [
                    {"name": "lipinski", "passed": True},
                    {"name": "veber",    "passed": True},
                    {"name": "ghose",    "passed": True},
                ]}
            ]
        }
        assert score_druglikeness(report) == pytest.approx(1.0)

    def test_two_of_three_pass_score_two_thirds(self):
        report = {
            "layers": [
                {"name": "druglikeness", "gates": [
                    {"name": "lipinski", "passed": True},
                    {"name": "veber",    "passed": True},
                    {"name": "ghose",    "passed": False},
                ]}
            ]
        }
        assert score_druglikeness(report) == pytest.approx(2.0 / 3.0)


class TestScoreAdmet:
    def test_no_red_flags_one(self):
        assert score_admet({"red_flags": []}) == 1.0

    def test_three_red_flags_drops_score(self):
        # 1 - 3*0.15 = 0.55
        assert score_admet({"red_flags": ["a", "b", "c"]}) == pytest.approx(0.55)

    def test_many_flags_clipped_at_zero(self):
        assert score_admet({"red_flags": ["a"] * 20}) == 0.0


class TestScoreSa:
    def test_sa_one_scores_one(self):
        assert score_synthesizability(1.0) == pytest.approx(1.0)

    def test_sa_six_scores_zero(self):
        assert score_synthesizability(6.0) == 0.0

    def test_sa_missing_half_credit(self):
        assert score_synthesizability(None) == 0.5


class TestScoreMechanism:
    def test_gate_failed_zero(self):
        assert score_mechanism({"passed": False}) == 0.0

    def test_gate_passed_base_half(self):
        assert score_mechanism(
            {"passed": True, "catalytic_residues_hit": []}
        ) == pytest.approx(0.5)

    def test_many_catalytic_hits_cap_one(self):
        assert score_mechanism(
            {"passed": True, "catalytic_residues_hit": [46, 50, 73, 76, 128]}
        ) == pytest.approx(1.0)


# ----------------------------------------------------------------------
# Composite score
# ----------------------------------------------------------------------


def _dummy_validate_all_pass():
    return {
        "passed": True,
        "layers": [
            {"name": "structure", "passed": True, "gates": [
                {"name": "parse", "passed": True},
                {"name": "heavy_atom_count", "passed": True},
                {"name": "net_formal_charge", "passed": True},
                {"name": "single_fragment", "passed": True},
            ]},
            {"name": "hygiene", "passed": True, "gates": [
                {"name": "pains", "passed": True},
                {"name": "reos", "passed": True},
            ]},
            {"name": "druglikeness", "passed": True, "gates": [
                {"name": "lipinski", "passed": True},
                {"name": "veber", "passed": True},
                {"name": "ghose", "passed": True},
            ]},
            {"name": "synthesizability", "passed": True, "gates": [
                {"name": "sa_score", "passed": True, "detail": {"sa_score": 3.0}},
            ]},
        ],
    }


def test_score_candidate_strong_all_gates(tmp_path: Path):
    s = score_candidate(
        target_yaml=TARGET_YAML,
        candidate_id="cand-good",
        smiles="c1ccncc1",
        docking={"best_energy_kcalmol": -11.0},
        novelty={"scaffold_hits": [], "max_tanimoto": 0.1},
        validate=_dummy_validate_all_pass(),
        admet={"red_flags": []},
        sa_score=3.0,
        mechanism={"passed": True, "catalytic_residues_hit": [46, 73]},
        redteam={"substantive_flaw": False},
    )
    # Expected per-layer scores:
    #   docking 1.0 * 0.30 = 0.30
    #   novelty 0.9 * 0.25 = 0.225
    #   druglikeness 1.0 * 0.15 = 0.15
    #   admet 1.0 * 0.15 = 0.15
    #   synth (6-3)/5 = 0.6 * 0.10 = 0.06
    #   mechanism 0.5 + 2*0.125 = 0.75 * 0.05 = 0.0375
    #   total = 0.9225  (above threshold 0.83)
    assert s.rigor_score == pytest.approx(0.9225, abs=1e-3)
    assert s.above_threshold is True
    assert s.veto_applied is False


def test_red_team_substantive_veto_zeroes_score():
    s = score_candidate(
        target_yaml=TARGET_YAML,
        candidate_id="cand-veto",
        smiles="c1ccncc1",
        docking={"best_energy_kcalmol": -11.0},
        novelty={"scaffold_hits": [], "max_tanimoto": 0.1},
        validate=_dummy_validate_all_pass(),
        admet={"red_flags": []},
        sa_score=3.0,
        mechanism={"passed": True, "catalytic_residues_hit": [46, 73]},
        redteam={"substantive_flaw": True},
    )
    assert s.rigor_score == 0.0
    assert s.veto_applied is True
    assert s.above_threshold is False


def test_score_candidate_below_threshold():
    # Strip everything to minimal scores; should fall below 0.83.
    s = score_candidate(
        target_yaml=TARGET_YAML,
        candidate_id="cand-weak",
        smiles="CCO",
        docking={"best_energy_kcalmol": -6.0},
        novelty={"scaffold_hits": [], "max_tanimoto": 0.5},
        validate=_dummy_validate_all_pass(),
        admet={"red_flags": ["flag1", "flag2"]},
        sa_score=5.5,
        mechanism={"passed": True, "catalytic_residues_hit": [73]},
        redteam={"substantive_flaw": False},
    )
    assert s.above_threshold is False
    assert s.rigor_score < 0.83


# ----------------------------------------------------------------------
# Leaderboard
# ----------------------------------------------------------------------


def test_update_leaderboard_inserts_and_sorts(tmp_path: Path):
    board_path = tmp_path / "leaderboard.json"

    def mk(cid: str, rigor: float, above: bool):
        return score_candidate(
            target_yaml=TARGET_YAML,
            candidate_id=cid,
            smiles="c1ccncc1",
            docking={"best_energy_kcalmol": -5.0 - 6.0 * rigor},  # linear to energy
            novelty={"scaffold_hits": [], "max_tanimoto": 1.0 - rigor},
            validate=_dummy_validate_all_pass(),
            admet={"red_flags": []},
            sa_score=3.0,
            mechanism={"passed": True, "catalytic_residues_hit": [46, 73]},
            redteam={"substantive_flaw": False},
        )

    s1 = mk("cand-1", 0.9, True)
    s2 = mk("cand-2", 0.95, True)
    s3 = mk("cand-3", 0.6, False)

    update_leaderboard(s1, leaderboard_path=board_path, top_n=10)
    update_leaderboard(s2, leaderboard_path=board_path, top_n=10)
    update_leaderboard(s3, leaderboard_path=board_path, top_n=10)

    data = json.loads(board_path.read_text())
    ids = [e["candidate_id"] for e in data["candidates"]]
    # cand-3 below threshold → not in leaderboard
    assert "cand-3" not in ids
    # sorted descending by rigor_score
    assert len(data["candidates"]) == 2
    scores = [e["rigor_score"] for e in data["candidates"]]
    assert scores == sorted(scores, reverse=True)


def test_update_leaderboard_replaces_same_id(tmp_path: Path):
    board_path = tmp_path / "leaderboard.json"

    first = score_candidate(
        target_yaml=TARGET_YAML,
        candidate_id="cand-x",
        smiles="c1ccncc1",
        docking={"best_energy_kcalmol": -11.0},
        novelty={"scaffold_hits": [], "max_tanimoto": 0.1},
        validate=_dummy_validate_all_pass(),
        admet={"red_flags": []},
        sa_score=3.0,
        mechanism={"passed": True, "catalytic_residues_hit": [46, 73]},
        redteam={"substantive_flaw": False},
    )
    update_leaderboard(first, leaderboard_path=board_path)
    data = json.loads(board_path.read_text())
    assert len(data["candidates"]) == 1
    first_score = data["candidates"][0]["rigor_score"]

    # Update same id with slightly worse mechanism (only 1 catalytic hit).
    second = score_candidate(
        target_yaml=TARGET_YAML,
        candidate_id="cand-x",
        smiles="c1ccncc1",
        docking={"best_energy_kcalmol": -11.0},
        novelty={"scaffold_hits": [], "max_tanimoto": 0.1},
        validate=_dummy_validate_all_pass(),
        admet={"red_flags": []},
        sa_score=3.0,
        mechanism={"passed": True, "catalytic_residues_hit": [73]},
        redteam={"substantive_flaw": False},
    )
    update_leaderboard(second, leaderboard_path=board_path)
    data = json.loads(board_path.read_text())
    assert len(data["candidates"]) == 1
    assert data["candidates"][0]["rigor_score"] != first_score


def test_score_candidate_is_deterministic():
    kwargs = dict(
        target_yaml=TARGET_YAML,
        candidate_id="cand-d",
        smiles="c1ccncc1",
        docking={"best_energy_kcalmol": -9.5},
        novelty={"scaffold_hits": [], "max_tanimoto": 0.15},
        validate=_dummy_validate_all_pass(),
        admet={"red_flags": ["one"]},
        sa_score=2.8,
        mechanism={"passed": True, "catalytic_residues_hit": [46, 50, 73]},
        redteam={"substantive_flaw": False},
    )
    s1 = score_candidate(**kwargs)
    s2 = score_candidate(**kwargs)
    assert s1.rigor_score == s2.rigor_score
    assert s1.sub_scores == s2.sub_scores
