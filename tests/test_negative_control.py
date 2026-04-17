"""Tests for loops/negative_control — pure EF1%, ROC-AUC, gate logic."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

from loops.negative_control import (  # noqa: E402
    ControlRun,
    EvaluationResult,
    RunRecord,
    build_alert,
    enrichment_factor,
    evaluate_record,
    load_decoys,
    load_history,
    roc_auc,
    append_record,
)


# ----------------------------------------------------------------------
# enrichment_factor
# ----------------------------------------------------------------------


def test_ef_perfect_enrichment():
    # 10 compounds, top 1% = top 1. If the single best is an active and
    # actives_total / N = 1/10 = 0.1, then EF = 1.0 / 0.1 = 10.0.
    ranked = [(-10.0, True)] + [(-5.0 - 0.1 * i, False) for i in range(9)]
    assert enrichment_factor(ranked, top_fraction=0.1) == pytest.approx(10.0)


def test_ef_zero_when_no_active_in_top():
    # Decoys sorted first (better scores), one active buried at the bottom
    ranked = [(-10.0 + 0.01 * i, False) for i in range(99)] + [(0.0, True)]
    assert enrichment_factor(ranked, top_fraction=0.01) == 0.0


def test_ef_handles_no_actives():
    ranked = [(-5.0 - i * 0.1, False) for i in range(100)]
    assert enrichment_factor(ranked, top_fraction=0.01) == 0.0


def test_ef_rounds_top_fraction_up():
    # 50 compounds, top_fraction = 0.01 → ceil(0.5) = 1 entry in top
    ranked = [(-10.0, True)] + [(-5.0 - 0.01 * i, False) for i in range(49)]
    assert enrichment_factor(ranked, top_fraction=0.01) == pytest.approx(50.0, abs=0.1)


# ----------------------------------------------------------------------
# roc_auc
# ----------------------------------------------------------------------


def test_roc_auc_perfect_separation():
    ranked = [(-10.0, True), (-9.5, True), (-5.0, False), (-4.5, False)]
    assert roc_auc(ranked) == pytest.approx(1.0)


def test_roc_auc_random_at_half():
    # Construct a list where one active is at the top and one at the
    # bottom: top active beats both decoys, bottom active beats none.
    # U = 2 + 0 = 2; denom = 2 * 2 = 4; AUC = 0.5.
    ranked = [(-10.0, True), (-8.0, False), (-5.0, False), (-1.0, True)]
    assert roc_auc(ranked) == pytest.approx(0.5)


def test_roc_auc_inverted_is_below_half():
    # Decoys scored better than actives ⇒ AUC < 0.5 (anti-correlation).
    ranked = [(-10.0, False), (-9.0, True), (-8.0, False), (-7.0, True)]
    auc = roc_auc(ranked)
    assert auc < 0.5
    assert auc == pytest.approx(0.25)


def test_roc_auc_all_ties():
    ranked = [(-5.0, True), (-5.0, False), (-5.0, True), (-5.0, False)]
    # All ties → 0.5
    assert roc_auc(ranked) == pytest.approx(0.5)


def test_roc_auc_handles_empty():
    assert roc_auc([]) == 0.0


def test_roc_auc_handles_single_class():
    assert roc_auc([(-1, True), (-2, True)]) == 0.0


# ----------------------------------------------------------------------
# Evaluate + alert
# ----------------------------------------------------------------------


def _good_record() -> RunRecord:
    # 5 actives clustered at the top (strong separation); 95 decoys at bottom
    runs = []
    for i in range(5):
        runs.append(ControlRun(
            name=f"active-{i}", delta_g_kcalmol=-10.0 - 0.1 * i, is_active=True,
        ))
    for i in range(95):
        runs.append(ControlRun(
            name=f"decoy-{i}", delta_g_kcalmol=-6.0 + 0.01 * i, is_active=False,
        ))
    return RunRecord(timestamp="t1", target_id="mrsa-gyrb", runs=runs)


def _broken_record() -> RunRecord:
    # Decoys and actives shuffled randomly — no discrimination
    runs = []
    for i in range(10):
        runs.append(ControlRun(
            name=f"active-{i}", delta_g_kcalmol=-5.0 - (i % 3), is_active=True,
        ))
    for i in range(90):
        runs.append(ControlRun(
            name=f"decoy-{i}", delta_g_kcalmol=-5.0 - (i % 5), is_active=False,
        ))
    return RunRecord(timestamp="t1", target_id="mrsa-gyrb", runs=runs)


def test_evaluate_good_record_passes():
    ev = evaluate_record(_good_record())
    assert ev.passed is True
    assert ev.roc_auc > 0.95
    assert ev.ef1 > 5.0


def test_evaluate_broken_record_fails_both_gates():
    ev = evaluate_record(_broken_record())
    assert ev.passed is False


def test_build_alert_clean_returns_none():
    ev = evaluate_record(_good_record())
    assert build_alert(ev) is None


def test_build_alert_broken_has_reasons():
    ev = evaluate_record(_broken_record())
    alert = build_alert(ev)
    assert alert is not None
    assert alert["severity"] == "halt"
    assert len(alert["reasons"]) >= 1


# ----------------------------------------------------------------------
# Decoys + history roundtrip
# ----------------------------------------------------------------------


def test_load_decoys_from_explicit_path(tmp_path: Path):
    p = tmp_path / "decoys.smi"
    p.write_text("CCO\nCCN CCN-name\n# comment\n\nc1ccccc1\n")
    decoys = load_decoys("mrsa-gyrb", explicit_path=p)
    # Expect 3 valid SMILES; comments and blank lines skipped; names trimmed.
    assert decoys == ["CCO", "CCN", "c1ccccc1"]


def test_history_roundtrip(tmp_path: Path):
    rec = _good_record()
    ev = evaluate_record(rec)
    rec.ef1 = ev.ef1
    rec.roc_auc = ev.roc_auc
    path = tmp_path / "negative-control.jsonl"
    append_record(rec, path=path)
    loaded = load_history(path)
    assert len(loaded) == 1
    assert loaded[0].ef1 == pytest.approx(rec.ef1)
    assert loaded[0].roc_auc == pytest.approx(rec.roc_auc)
    assert loaded[0].runs[0].is_active is True
