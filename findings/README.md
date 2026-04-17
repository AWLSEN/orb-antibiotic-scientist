# findings/

Everything the agent produces goes here. This directory is the agent's
only persistent memory across auto-compaction / restart boundaries —
treat it as the source of truth.

## Layout

```
findings/
├── leaderboard.json            # current top-N candidates, updated by src/scoring.py
├── candidates/                 # one directory per validated candidate
│   └── cand-YYYYMMDD-xxxxxx/
│       ├── candidate.json      # id + SMILES + design rationale
│       ├── validate.json       # Layers 1-4 (structure, hygiene, druglike, SA)
│       ├── docking.json        # Layer 5 primary Vina
│       ├── docking-secondary.json  # Layer 6 secondary (Vina-seed2 / DiffDock)
│       ├── mechanism.json      # Layer 7 pocket-contact analysis
│       ├── admet.json          # Layer 8 QED + ADMET + toxicophores
│       ├── novelty.json        # Layer 9 ChEMBL + scaffold blocklist
│       ├── retrosynthesis.json # Layer 11 BRICS/AiZynth
│       ├── redteam.json        # Layer 12 adversarial critique
│       └── scored.json         # composite rigor score + sub-scores
├── docking/                    # legacy docking outputs (.pdbqt poses — gitignored)
├── admet/                      # freestanding ADMET reports outside a candidate
├── novelty/                    # standalone novelty reports
├── mechanism/                  # standalone mechanism reports
├── red-team/                   # <cid>.json + <cid>.md red-team critiques
├── retrosynthesis/             # standalone retrosynthesis outputs
├── weekly-reports/             # weekly agent-authored report Markdown
└── loop-health/                # time series for the pipeline-health loops
    ├── positive-control.jsonl
    ├── positive-control-alert.json   # present iff the loop failed
    ├── negative-control.jsonl
    ├── negative-control-alert.json
    └── consensus.jsonl
```

## What's NOT here

- Raw PDB / PDBQT pose files (large; `.gitignored`) live under each
  `findings/candidates/<id>/` during runtime but are excluded from the
  commit stream. The `docking.json` summary is what's versioned.
- The ChEMBL similarity cache lives in `data/chembl-cache/`.

## First run

The agent bootstraps the structure on start — see the `FINDING_DIRS`
list in `src/agent.py`. On a fresh deploy, `leaderboard.json` is empty
and the `.jsonl` files don't exist yet. The dashboard handles this state
gracefully (shows "no data yet" placeholders).

## Candidate ID format

`cand-YYYYMMDD-xxxxxx` where `xxxxxx` is the first 6 hex chars of
`sha256(canonical_smiles_of_candidate)`. Deterministic: the same
molecule always gets the same ID, preventing duplicate work across
session boundaries.
