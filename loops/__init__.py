"""Continuous pipeline-health loops for orb-antibiotic-scientist.

Each module in this package is a periodic job that verifies the
end-to-end docking / scoring pipeline is still producing trustworthy
results while the main agent keeps generating candidates.

  positive_control   hourly — re-dock known inhibitors, track rank/ΔG drift
  negative_control   hourly — actives vs decoys EF1% + ROC-AUC
  consensus          weekly — cross-method docking agreement
  retrosynthesis     per-top-candidate — synthesis feasibility check
"""
