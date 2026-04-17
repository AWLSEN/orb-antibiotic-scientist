# orb-antibiotic-scientist — agent context

## What This Is

An autonomous AI research agent that designs and computationally validates novel antibiotic candidates for drug-resistant bacterial targets. Built on the Claude Agent SDK (Opus), the same harness pattern as SPOQ-Food (`null-bytes/ai-nutrition-meat-pipeline`), adapted for medicinal chemistry against antibacterial targets.

Initial target: MRSA DNA gyrase B (PDB 4URL).
Planned target #2: *M. tuberculosis* InhA (PDB 4TZK).

## Architecture

- `src/agent.py` — Claude Agent SDK runner (Opus, auto-compaction, session resume)
- `src/agent-prompt.md` — system prompt: mission, scoring rubric, output contract
- `src/watchdog.sh` — process supervisor (stall detection, crash-restart)
- `src/control.sh` — start / stop / status / logs / findings
- `src/tools/` — per-layer validation tools (SMILES validate, Vina, DiffDock, SwissADME, ChEMBL novelty, PLIP mechanism, AiZynthFinder retrosynthesis, red-team)
- `loops/` — continuous pipeline-health loops (positive/negative controls, consensus, novelty-drift)
- `targets/` — YAML specs for each bacterial target
- `findings/` — agent output (candidates, docking, admet, novelty, mechanism, red-team, leaderboard)

## Running

```bash
./src/control.sh start    # Start agent in background
./src/control.sh stop     # Stop agent
./src/control.sh status   # Status + stats
./src/control.sh logs     # Tail watchdog logs
./src/control.sh findings # Show findings summary
```

## Key Rules

1. **Every candidate must pass all 12 verification layers** before appearing on the leaderboard. See README.md §Verification.
2. **Positive and negative controls run continuously.** If known inhibitors stop ranking in the top 5%, or actives/decoys separation (EF1%, ROC-AUC) drops below threshold, halt and alert — the pipeline is broken.
3. **Novelty floor is Tanimoto < 0.4** vs. ChEMBL (Morgan fingerprint, radius 2). Rediscovering known antibiotics is a failure mode.
4. **Every candidate carries a red-team critique** — a separate Claude call with an adversarial prompt. No red team = no leaderboard entry.
5. **Save everything to disk immediately.** Each candidate is a directory with SMILES, poses, ADMET, novelty, mechanism, red team — fully reproducible.
6. **Do not repeat work.** Check `findings/candidates/` before proposing. Same SMILES + same target = same scores.
7. **Use structured formats** — JSON, JSONL, CSV, SMILES, FASTA, PDB.
8. **Nothing therapeutic is claimed.** Outputs are computational hypotheses. See DISCLAIMER.md.

## Verification Loops (summary)

See README.md for the full 12-layer filter chain and continuous-loop description. Short version:

- Per-candidate filter chain (synchronous): SMILES validity → PAINS/REOS → druglike → SA Score → Vina → DiffDock/GNINA → PLIP mechanism → ADMET → novelty → resistance-proof → retrosynthesis → red team.
- Hourly: positive-control inhibitor ranks, negative-control DUD-E decoys EF1% + ROC-AUC.
- Weekly: cross-method consensus, novelty drift against ChEMBL updates, weekly report.

## External Services

- **AutoDock Vina** (local binary) — primary docking
- **DiffDock / GNINA** (local or via API) — secondary docking for consensus
- **ChEMBL API** — novelty + known-antibiotic lookups
- **PubChem / DrugBank / ZINC** — additional similarity searches
- **PubMed / bioRxiv** (Tavily-backed) — literature grounding
- **SwissADME / ProTox-II / pkCSM** — ADMET predictions
- **AiZynthFinder** — retrosynthesis feasibility
- **PLIP** — interaction fingerprinting

## Cadence

Commit after every small step. Push after every 5 commits. PRs on `nextbysam` fork, not upstream org.

## Related

- SPOQ-Food: `/Users/sam/workspaces/programming/null-bytes/ai-nutrition-meat-pipeline` — reference harness.
- Orb Cloud docs: https://docs.orbcloud.dev — deployment target.
