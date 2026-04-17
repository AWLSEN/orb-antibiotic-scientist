"""Tests for src/tools/redteam (Layer 12).

Covers the dry-run path (deterministic, CI-safe) and the LLM-response
JSON parser. Live mode (`_call_claude_adversarial`) is not exercised in
CI to avoid network flake + API-key requirements.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT / "src"))

from tools.redteam import (  # noqa: E402
    parse_redteam_json,
    red_team,
)


# ----------------------------------------------------------------------
# parse_redteam_json
# ----------------------------------------------------------------------


SAMPLE_RESPONSE = '''Here is my adversarial review.
{
  "flaws": [
    {"category": "pharmacokinetics", "description": "LogP too high, poor oral bioavail",
     "severity": 8, "substantive": true},
    {"category": "synthesis", "description": "Retrosynthesis requires hazardous triflate",
     "severity": 6, "substantive": false}
  ],
  "reasoning": "Overall the candidate is marginal."
}
end of response'''


def test_parse_redteam_json_happy_path():
    flaws, reasoning = parse_redteam_json(SAMPLE_RESPONSE)
    assert len(flaws) == 2
    assert flaws[0].category == "pharmacokinetics"
    assert flaws[0].severity == 8
    assert flaws[0].substantive is True
    assert "marginal" in reasoning


def test_parse_redteam_json_raises_on_missing_json():
    with pytest.raises(ValueError):
        parse_redteam_json("no JSON at all here, just prose")


# ----------------------------------------------------------------------
# dry-run critique — rule outcomes
# ----------------------------------------------------------------------


def _minimal_good_dossier() -> dict:
    return {
        "candidate_id": "cand-000",
        "docking": {"best_energy_kcalmol": -9.1},
        "admet": {
            "logp": 3.0,
            "bbb_likely": False,
            "toxicophore_hits": [],
        },
        "novelty": {
            "scaffold_hits": [],
            "max_tanimoto": 0.22,
            "nearest_chembl_id": "CHEMBL999",
        },
        "mechanism": {
            "num_contacts": 4,
            "minimum_required": 2,
            "resistance_only": False,
        },
        "sa_score": 3.1,
    }


def test_dry_run_passes_clean_candidate(tmp_path: Path):
    d = _minimal_good_dossier()
    r = red_team(d, write_dir=tmp_path)
    # Must raise ≥3 flaws even on a clean dossier (explicit rule), but none
    # should be substantive.
    assert len(r.flaws) >= 3
    assert r.substantive_flaw is False
    assert r.passed is True
    assert r.mode == "dry-run"


def test_dry_run_flags_weak_docking(tmp_path: Path):
    d = _minimal_good_dossier()
    d["docking"]["best_energy_kcalmol"] = -6.0  # weak
    r = red_team(d, write_dir=tmp_path)
    assert r.substantive_flaw is True
    assert r.passed is False
    assert any("unlikely to show measurable activity" in f.description for f in r.flaws)


def test_dry_run_flags_lipophilic_candidate(tmp_path: Path):
    d = _minimal_good_dossier()
    d["admet"]["logp"] = 6.5
    r = red_team(d, write_dir=tmp_path)
    assert any(f.category == "pharmacokinetics" and f.substantive for f in r.flaws)
    assert r.passed is False


def test_dry_run_flags_scaffold_rediscovery(tmp_path: Path):
    d = _minimal_good_dossier()
    d["novelty"]["scaffold_hits"] = ["fluoroquinolone_core"]
    r = red_team(d, write_dir=tmp_path)
    assert any(f.category == "novelty" and f.substantive for f in r.flaws)
    assert r.passed is False


def test_dry_run_flags_resistance_only(tmp_path: Path):
    d = _minimal_good_dossier()
    d["mechanism"]["resistance_only"] = True
    r = red_team(d, write_dir=tmp_path)
    assert any(f.category == "resistance" and f.substantive for f in r.flaws)
    assert r.passed is False


def test_dry_run_flags_synthesis_hard(tmp_path: Path):
    d = _minimal_good_dossier()
    d["sa_score"] = 7.5
    r = red_team(d, write_dir=tmp_path)
    assert any(f.category == "synthesis" and f.substantive for f in r.flaws)


def test_dry_run_missing_docking_is_substantive(tmp_path: Path):
    d = _minimal_good_dossier()
    d["docking"] = {}  # no best_energy
    r = red_team(d, write_dir=tmp_path)
    assert any("docking" in f.description.lower() for f in r.flaws)
    assert r.substantive_flaw is True


def test_red_team_persists_markdown(tmp_path: Path):
    d = _minimal_good_dossier()
    r = red_team(d, write_dir=tmp_path)
    md_path = tmp_path / "cand-000.md"
    assert md_path.exists()
    body = md_path.read_text()
    assert "# Red-Team Critique" in body
    assert "Substantive flaw" in body
    json_path = tmp_path / "cand-000.json"
    assert json_path.exists()
    # Round-trip JSON
    parsed = json.loads(json_path.read_text())
    assert parsed["candidate_id"] == "cand-000"


def test_dossier_hash_is_deterministic(tmp_path: Path):
    d = _minimal_good_dossier()
    r1 = red_team(d, write_dir=tmp_path)
    r2 = red_team(d, write_dir=tmp_path)
    assert r1.dossier_hash == r2.dossier_hash
