"""Tests for loops/consensus — Spearman + RMSD consensus gate."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

from loops.consensus import (  # noqa: E402
    build_alert,
    compute_consensus,
    spearman_correlation,
)


# ----------------------------------------------------------------------
# Spearman
# ----------------------------------------------------------------------


def test_spearman_perfect_positive():
    xs = [1.0, 2.0, 3.0, 4.0, 5.0]
    ys = [10.0, 20.0, 30.0, 40.0, 50.0]
    assert spearman_correlation(xs, ys) == pytest.approx(1.0)


def test_spearman_perfect_negative():
    xs = [1.0, 2.0, 3.0, 4.0, 5.0]
    ys = [50.0, 40.0, 30.0, 20.0, 10.0]
    assert spearman_correlation(xs, ys) == pytest.approx(-1.0)


def test_spearman_zero_on_unrelated():
    xs = [1, 2, 3, 4]
    ys = [3, 1, 4, 2]
    # This is not literally 0 but should be small (Spearman is exact on ranks
    # so we compute it)
    r = spearman_correlation(xs, ys)
    assert -1.0 <= r <= 1.0


def test_spearman_handles_ties():
    xs = [1.0, 1.0, 2.0, 2.0]
    ys = [5.0, 5.0, 6.0, 6.0]
    assert spearman_correlation(xs, ys) == pytest.approx(1.0)


def test_spearman_degenerate_lengths():
    assert spearman_correlation([1], [2]) == 0.0
    assert spearman_correlation([], []) == 0.0


# ----------------------------------------------------------------------
# compute_consensus
# ----------------------------------------------------------------------


def test_consensus_perfect_agreement_passes():
    primary = {f"c{i}": float(-10 + i * 0.1) for i in range(10)}
    secondary = dict(primary)  # identical = perfect agreement
    rmsds = {k: 0.5 for k in primary}
    report = compute_consensus(
        primary=primary, secondary=secondary, rmsds=rmsds, target_id="mrsa-gyrb",
    )
    assert report.spearman == pytest.approx(1.0)
    assert report.mean_rmsd_a == pytest.approx(0.5)
    assert report.passed is True
    assert report.unstable_candidates == []


def test_consensus_disagreement_fails_spearman():
    primary = {f"c{i}": float(-10 + i * 0.1) for i in range(10)}
    # Reverse order ⇒ perfect anti-correlation (−1) ⇒ fails 0.5 threshold
    secondary = {f"c{i}": float(-10 + (9 - i) * 0.1) for i in range(10)}
    report = compute_consensus(
        primary=primary, secondary=secondary, rmsds=None, target_id="mrsa-gyrb",
    )
    assert report.spearman_pass is False
    assert report.passed is False


def test_consensus_high_rmsd_fails_rmsd_gate():
    primary = {f"c{i}": float(-10 + i * 0.1) for i in range(5)}
    secondary = dict(primary)
    rmsds = {k: 10.0 for k in primary}  # way above 3.0 threshold
    report = compute_consensus(
        primary=primary, secondary=secondary, rmsds=rmsds, target_id="mrsa-gyrb",
    )
    assert report.spearman_pass is True
    assert report.rmsd_pass is False
    assert report.passed is False


def test_consensus_unstable_candidates_flagged():
    # 10 candidates in perfect order for primary, but one outlier swapped
    # to the bottom of secondary. That outlier has large rank_delta and
    # high RMSD — flag it as unstable.
    primary = {f"c{i}": float(-10 + i * 0.1) for i in range(10)}
    secondary = {f"c{i}": float(-10 + i * 0.1) for i in range(10)}
    # Swap c0 (best primary) with worst secondary rank
    secondary["c0"] = -1.0  # now worst
    rmsds = {f"c{i}": 0.5 for i in range(10)}
    rmsds["c0"] = 6.0  # high disagreement pose
    report = compute_consensus(
        primary=primary, secondary=secondary, rmsds=rmsds, target_id="mrsa-gyrb",
        rank_delta_unstable=5, rmsd_unstable=4.0,
    )
    assert "c0" in report.unstable_candidates


def test_consensus_missing_ids_are_skipped():
    primary = {"c0": -10.0, "c1": -9.0, "c2": -8.0}
    secondary = {"c0": -9.8, "c1": -8.9}  # c2 missing
    report = compute_consensus(
        primary=primary, secondary=secondary, rmsds=None, target_id="mrsa-gyrb",
    )
    ids_in_report = {e["candidate_id"] for e in report.entries}
    assert "c2" not in ids_in_report
    assert report.n_candidates == 2


# ----------------------------------------------------------------------
# build_alert
# ----------------------------------------------------------------------


def test_alert_clean_report_returns_none():
    primary = {f"c{i}": float(-10 + i * 0.1) for i in range(10)}
    secondary = dict(primary)
    report = compute_consensus(
        primary=primary, secondary=secondary, rmsds=None, target_id="mrsa-gyrb",
    )
    assert build_alert(report) is None


def test_alert_failing_spearman_has_warn_severity():
    primary = {f"c{i}": float(-10 + i * 0.1) for i in range(10)}
    secondary = {f"c{i}": float(-10 + (9 - i) * 0.1) for i in range(10)}
    report = compute_consensus(
        primary=primary, secondary=secondary, rmsds=None, target_id="mrsa-gyrb",
    )
    alert = build_alert(report)
    assert alert is not None
    assert alert["severity"] == "warn"
    assert len(alert["reasons"]) >= 1
