"""Tests for loops/positive_control. Focus on the pure drift/alert logic."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

from loops.positive_control import (  # noqa: E402
    DriftEvaluation,
    ReferenceRun,
    RunRecord,
    append_record,
    build_alert,
    evaluate_drift,
    load_history,
    rolling_mean_by_reference,
)


def _mk(t: str, refs: dict[str, float | None]) -> RunRecord:
    runs = [ReferenceRun(name=n, delta_g_kcalmol=v) for n, v in refs.items()]
    return RunRecord(timestamp=t, target_id="mrsa-gyrb", runs=runs)


def test_rolling_mean_handles_window():
    history = [
        _mk("t1", {"novobiocin": -10.0}),
        _mk("t2", {"novobiocin": -10.2}),
        _mk("t3", {"novobiocin": -10.1}),
        _mk("t4", {"novobiocin": -10.4}),
    ]
    means = rolling_mean_by_reference(history, window=3)
    # Last 3 only: (-10.2, -10.1, -10.4) → mean = -10.2333
    assert means["novobiocin"][0] == pytest.approx(-10.2333, abs=1e-3)
    assert means["novobiocin"][1] == 3


def test_rolling_mean_skips_missing():
    history = [
        _mk("t1", {"novobiocin": -10.0}),
        _mk("t2", {"novobiocin": None}),
        _mk("t3", {"novobiocin": -10.2}),
    ]
    means = rolling_mean_by_reference(history, window=10)
    assert means["novobiocin"][1] == 2


def test_evaluate_drift_under_threshold():
    history = [
        _mk("t1", {"novobiocin": -10.0}),
        _mk("t2", {"novobiocin": -10.2}),
        _mk("t3", {"novobiocin": -10.1}),
    ]
    current = _mk("t4", {"novobiocin": -9.5})
    evals = evaluate_drift(current, history, max_drift_kcalmol=1.5)
    assert len(evals) == 1
    e = evals[0]
    assert e.drift_abs == pytest.approx(abs(-9.5 - ((-10.0 - 10.2 - 10.1) / 3)), abs=1e-3)
    assert e.over_threshold is False


def test_evaluate_drift_over_threshold_raises_flag():
    history = [
        _mk("t1", {"novobiocin": -10.0}),
        _mk("t2", {"novobiocin": -10.2}),
        _mk("t3", {"novobiocin": -10.1}),
    ]
    # Drift > 1.5 kcal/mol
    current = _mk("t4", {"novobiocin": -8.0})
    evals = evaluate_drift(current, history, max_drift_kcalmol=1.5)
    assert evals[0].over_threshold is True


def test_evaluate_drift_bootstrap_tolerant():
    # With < 3 history points, drift should NOT be flagged even if large.
    history = [_mk("t1", {"novobiocin": -10.0})]
    current = _mk("t2", {"novobiocin": -5.0})
    evals = evaluate_drift(current, history, max_drift_kcalmol=1.5)
    assert evals[0].over_threshold is False
    assert evals[0].rolling_mean is None


def test_evaluate_drift_failed_dock_is_over_threshold():
    history = [_mk(f"t{i}", {"novobiocin": -10.0}) for i in range(5)]
    current = _mk("tnow", {"novobiocin": None})
    evals = evaluate_drift(current, history)
    assert evals[0].over_threshold is True
    assert evals[0].current is None


def test_build_alert_from_evaluations():
    evals = [
        DriftEvaluation(
            name="novobiocin", current=-7.0, rolling_mean=-10.0,
            rolling_n=20, drift_abs=3.0, over_threshold=True,
            threshold=1.5,
        ),
        DriftEvaluation(
            name="clorobiocin", current=-10.1, rolling_mean=-10.0,
            rolling_n=20, drift_abs=0.1, over_threshold=False,
            threshold=1.5,
        ),
    ]
    alert = build_alert(evals, target_id="mrsa-gyrb")
    assert alert is not None
    assert alert["severity"] == "halt"
    assert len(alert["flagged"]) == 1
    assert alert["flagged"][0]["name"] == "novobiocin"


def test_build_alert_empty_when_clean():
    evals = [
        DriftEvaluation(
            name="novobiocin", current=-10.1, rolling_mean=-10.0,
            rolling_n=20, drift_abs=0.1, over_threshold=False,
            threshold=1.5,
        ),
    ]
    assert build_alert(evals, target_id="mrsa-gyrb") is None


def test_jsonl_roundtrip(tmp_path: Path):
    path = tmp_path / "positive-control.jsonl"
    rec1 = _mk("t1", {"novobiocin": -10.0, "clorobiocin": -10.5})
    rec2 = _mk("t2", {"novobiocin": -10.2, "clorobiocin": -10.4})
    append_record(rec1, path=path)
    append_record(rec2, path=path)

    loaded = load_history(path)
    assert len(loaded) == 2
    assert loaded[0].timestamp == "t1"
    assert loaded[1].runs[0].name == "novobiocin"
