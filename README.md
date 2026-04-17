# orb-antibiotic-scientist

An autonomous AI scientist that designs novel antibiotic candidates for drug-resistant bacteria. Runs 24/7 on Orb Cloud with Claude Opus via the Agent SDK. Sleeps between LLM calls. Wakes when results arrive.

> Antibiotic resistance kills ~1.27 million people every year. Big pharma left the space in the 2010s — the economics don't work at their scale. We need R&D at 100× lower cost.

## What it does

The agent runs a continuous design-and-validate loop, 24 hours a day:

1. **Designs** a novel candidate — small molecule (RDKit scaffold hopping) or antimicrobial peptide (ESM-2).
2. **Validates** the candidate through 12 orthogonal computational filters (see [Verification](#verification)).
3. **Scores** rigour + binding + ADMET + novelty into a composite leaderboard rank.
4. **Documents** the candidate with SMILES / FASTA, 3D pose, mechanism hypothesis, and literature citations.
5. **Critiques** itself — a separate adversarial Claude call tries to find reasons the candidate will fail.
6. **Sleeps** until the next docking job or LLM call returns.

Initial target: **MRSA DNA gyrase B (GyrB)**, PDB `4URL`.
Target #2 (planned): *M. tuberculosis* **InhA**, PDB `4TZK`.

## What this is NOT

**These are computational candidates, not validated therapies.** They have not been synthesized, tested in any biological assay, or reviewed by any regulatory body. See [DISCLAIMER.md](./DISCLAIMER.md) before drawing any conclusion from the output.

## Verification

Every candidate passes through a 12-layer filter chain before it appears on the leaderboard:

| # | Layer | Tool | Gate |
|---|---|---|---|
| 1 | SMILES validity | RDKit parse | valid atoms + valence |
| 2 | Structural hygiene | PAINS + REOS | no pan-assay-interference scaffolds |
| 3 | Drug-likeness | Lipinski + Veber + Ghose | all three pass |
| 4 | Synthesizability | SA Score (RDKit) | SA ≤ 6 |
| 5 | Primary docking | AutoDock Vina | ΔG < −8.0 kcal/mol |
| 6 | Secondary docking | DiffDock / GNINA | pose consensus with Vina |
| 7 | Mechanism sanity | PLIP interaction fingerprint | contacts ≥ 2 catalytic residues |
| 8 | ADMET | SwissADME + ProTox-II + pkCSM | no critical red flags |
| 9 | Novelty | ChEMBL Tanimoto (Morgan r=2) | max similarity < 0.4 |
| 10 | Resistance-proof | target hotspot lookup | not hitting known resistance site alone |
| 11 | Retrosynthesis | AiZynthFinder | ≥ 1 feasible route in ≤ 4 steps |
| 12 | Red team | adversarial Claude call | no substantive flaw found |

Beyond the per-candidate chain, **continuous pipeline-health loops** run in the background:

- **Positive-control loop** (hourly) — re-docks known GyrB inhibitors (novobiocin, GSK299423, clorobiocin); their rank must stay in the top 5%.
- **Negative-control loop** (hourly) — docks DUD-E property-matched decoys; the pipeline must maintain EF1% ≥ 5 and ROC-AUC ≥ 0.75 on the actives/decoys split.
- **Cross-method consensus** (weekly) — top 20 re-docked with Vina + DiffDock + GNINA; disagreements published.
- **Novelty-drift loop** (weekly) — newly published antibiotics pulled from ChEMBL; existing candidates re-scored for novelty.
- **Self-consistency loop** — top candidates re-scored 3× with scrambled prompts; high variance demotes.
- **Reproducibility loop** — CI re-runs scoring on every candidate; same SMILES + target must produce same scores.

All loop metrics are surfaced on the public dashboard in real time. The whole point is that *you can watch the pipeline verify itself*.

## Architecture

```
orb-antibiotic-scientist/
├── src/
│   ├── agent.py               # Claude Agent SDK runner (Opus, auto-compact, resume)
│   ├── agent-prompt.md        # System prompt: mission, rubric, output contract
│   ├── watchdog.sh            # Process supervisor (restart on crash / 5h stall)
│   ├── control.sh             # start / stop / status / logs / findings
│   └── tools/
│       ├── validate_smiles.py
│       ├── dock_vina.py
│       ├── dock_diffdock.py
│       ├── admet.py
│       ├── novelty.py
│       ├── mechanism.py
│       ├── retrosynthesis.py
│       └── redteam.py
├── loops/
│   ├── positive_control.py
│   ├── negative_control.py
│   ├── consensus.py
│   └── novelty_drift.py
├── targets/
│   ├── mrsa-gyrb.yaml
│   └── tb-inha.yaml
├── findings/
│   ├── candidates/            # SMILES + metadata
│   ├── docking/               # Vina poses + scores
│   ├── admet/                 # SwissADME / ProTox-II / pkCSM reports
│   ├── novelty/               # ChEMBL similarity reports
│   ├── mechanism/             # Binding-pocket interaction analysis
│   ├── red-team/              # Adversarial critiques
│   ├── leaderboard.json       # Top candidates, updated hourly
│   └── weekly-reports/        # Agent-written weekly summaries
├── dashboard/                 # Next.js public site
├── orb.toml                   # Orb Cloud deploy config
└── DISCLAIMER.md
```

## Running the agent

```bash
# Start the agent (runs in background, survives terminal close)
./src/control.sh start

# Check if it's running and see stats
./src/control.sh status

# Tail logs
./src/control.sh logs

# See what it's found
./src/control.sh findings

# Stop the agent
./src/control.sh stop
```

The agent picks up where it left off by inspecting `findings/` on restart.

## Deployment

Runs on [Orb Cloud](https://orbcloud.dev). `orb.toml` configures checkpoint-to-NVMe during LLM waits and docking jobs — the agent hibernates between bursts of activity, so continuous operation costs fractions of a cent per hour.

## License

MIT — see [LICENSE](./LICENSE).

## Contributing

Academic and industry wet labs interested in testing top-ranked candidates: please open an issue or email the maintainer. The top 5 candidates (by composite rigor score) each week are prepared as SMILES + mechanism brief + suggested assay for wet-lab validation.

---

*Built as part of the Orb Cloud open-source repo strategy. The Orb-specific moat here is that an autonomous scientist with continuous pipeline-health loops is economically viable only when compute sleeps during the ~95% of time the agent waits on docking jobs and LLM responses.*
