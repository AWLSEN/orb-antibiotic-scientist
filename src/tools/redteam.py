#!/usr/bin/env python3
"""
Layer 12: adversarial red-team critique.

A dedicated Claude call with a deliberately critical system prompt reviews
the full candidate dossier (SMILES + docking + ADMET + novelty + mechanism
+ retrosynthesis) and is asked to find reasons the candidate will fail in
a wet lab. The returned critique is stored under
`findings/red-team/<candidate_id>.json` alongside a human-readable
Markdown version in `findings/red-team/<candidate_id>.md`.

Gate: candidates with `substantive_flaw == true` are demoted off the
leaderboard regardless of how well they did in earlier layers. The
`severity_score` (0–10) is surfaced on the leaderboard for transparency.

This tool deliberately uses a *separate* prompt + separate API call from
the main agent harness so its reasoning is independent. The main agent is
optimising for candidates that pass filters; the red-team model is
optimising for finding reasons candidates fail.

Two execution modes:
  - live: calls the Anthropic API with a claude-opus model.
  - dry-run (default in tests): returns a pre-canned critique based on
    deterministic rules applied to the dossier (e.g. "LogP > 5 → flag
    membrane issue"). Used in CI to exercise the parsing + storage path.

CLI:
  python -m tools.redteam --dossier findings/candidates/<id>/dossier.json
  python -m tools.redteam --dossier <path> --live

Exit codes: 0 pass (no substantive flaw), 1 fail, 2 setup error.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from textwrap import dedent
from typing import Any

REPO_ROOT = Path(__file__).parent.parent.parent


# ----------------------------------------------------------------------
# Data classes
# ----------------------------------------------------------------------


@dataclass
class RedTeamFlaw:
    category: str            # e.g. "pharmacokinetics", "chemistry", "mechanism"
    description: str
    severity: int            # 1-10
    substantive: bool


@dataclass
class RedTeamResult:
    candidate_id: str
    dossier_hash: str
    flaws: list[RedTeamFlaw] = field(default_factory=list)
    substantive_flaw: bool = False
    severity_score: int = 0
    passed: bool = False
    mode: str = "dry-run"
    model: str | None = None
    reasoning: str = ""
    error: str | None = None
    created_at: str = ""


# ----------------------------------------------------------------------
# Prompts
# ----------------------------------------------------------------------


ADVERSARIAL_SYSTEM = dedent("""
You are a senior medicinal-chemistry red-team reviewer. Your ONE job is to
find reasons the candidate in the dossier will FAIL in a wet-lab experiment.
You are adversarial by design: the other reviewers are optimising for
candidates passing filters — you are optimising for catching the failures
they will miss. You review thousands of candidates every year. Be specific,
concrete, and mechanistic.

Output format — MUST be valid JSON only, no prose before or after:

{
  "flaws": [
    {
      "category": "pharmacokinetics" | "chemistry" | "mechanism" | "synthesis"
                  | "resistance" | "safety" | "novelty" | "other",
      "description": "<1-3 sentence precise description>",
      "severity": <integer 1-10>,
      "substantive": <true if this alone would likely sink the candidate>
    }
    // up to 5 flaws
  ],
  "reasoning": "<1 paragraph explaining your overall assessment>"
}

Rules:
- Find AT LEAST 3 flaws. It is almost impossible for a computational
  candidate to be flawless.
- `substantive: true` means this flaw alone is a high-probability killer
  (e.g. hERG liability, resistance would evolve in < 10 generations,
  synthesis requires a known hazardous reagent, the docking pose violates
  the mechanism of the target).
- Do not invent facts not in the dossier.
- If the dossier lacks key data (e.g. no ADMET section), that absence is
  itself a flaw of severity 6–7.
""").strip()


USER_TEMPLATE = dedent("""
Candidate dossier:

```json
{dossier_json}
```

Review this candidate adversarially and return JSON per the system prompt.
""").strip()


# ----------------------------------------------------------------------
# Dry-run: deterministic rule-based critique
# ----------------------------------------------------------------------


def _dry_run_flaws(dossier: dict[str, Any]) -> tuple[list[RedTeamFlaw], str]:
    """Deterministic red-team critique used in CI/tests. Applies a handful
    of concrete rules so tests can assert on stable output."""
    flaws: list[RedTeamFlaw] = []
    notes: list[str] = []

    # Docking-energy sanity
    docking = dossier.get("docking") or {}
    best_energy = docking.get("best_energy_kcalmol")
    if best_energy is None:
        flaws.append(RedTeamFlaw(
            category="other",
            description="No docking result provided; the candidate cannot be evaluated for binding.",
            severity=9, substantive=True,
        ))
    elif best_energy > -6.5:
        flaws.append(RedTeamFlaw(
            category="mechanism",
            description=(
                f"Vina ΔG {best_energy:.2f} kcal/mol is weak; "
                "candidate unlikely to show measurable activity in wet assay."
            ),
            severity=8, substantive=True,
        ))
    elif best_energy > -8.0:
        flaws.append(RedTeamFlaw(
            category="mechanism",
            description=(
                f"Vina ΔG {best_energy:.2f} kcal/mol is marginal — "
                "known inhibitors typically score ≤ −8.5."
            ),
            severity=5, substantive=False,
        ))

    # ADMET red flags
    admet = dossier.get("admet") or {}
    logp = admet.get("logp")
    if logp is not None and logp > 5.0:
        flaws.append(RedTeamFlaw(
            category="pharmacokinetics",
            description=(
                f"LogP {logp:.2f} is too lipophilic; expect poor solubility, "
                "high plasma-protein binding, and likely hERG liability."
            ),
            severity=7, substantive=True,
        ))
    if admet.get("bbb_likely") is True:
        notes.append("BBB-permeant candidate; for a peripheral MRSA target this is a distribution liability.")
    tox_hits = admet.get("toxicophore_hits") or []
    for t in tox_hits:
        flaws.append(RedTeamFlaw(
            category="safety",
            description=f"Toxicophore present: {t}. Known to correlate with wet-lab cytotoxicity.",
            severity=7, substantive=True,
        ))

    # Novelty: flag scaffold blocklist hits as high severity (rediscovery).
    novelty = dossier.get("novelty") or {}
    scaffold_hits = novelty.get("scaffold_hits") or []
    for s in scaffold_hits:
        flaws.append(RedTeamFlaw(
            category="novelty",
            description=(
                f"Candidate contains {s}; this is a re-derivative of a known "
                "antibiotic class — resistance mechanisms already exist in the wild."
            ),
            severity=8, substantive=True,
        ))
    max_tc = novelty.get("max_tanimoto")
    if isinstance(max_tc, (int, float)) and max_tc >= 0.7:
        flaws.append(RedTeamFlaw(
            category="novelty",
            description=(
                f"Max Tanimoto {max_tc:.2f} vs. ChEMBL nearest neighbour "
                f"({novelty.get('nearest_chembl_name') or novelty.get('nearest_chembl_id')}) "
                "is too close — effectively a known compound."
            ),
            severity=8, substantive=True,
        ))

    # Mechanism
    mech = dossier.get("mechanism") or {}
    if mech.get("resistance_only") is True:
        flaws.append(RedTeamFlaw(
            category="resistance",
            description=(
                "Pose contacts only residues listed as single-site resistance "
                "hotspots; a single point mutation would abolish binding."
            ),
            severity=9, substantive=True,
        ))
    if mech.get("num_contacts", 0) < mech.get("minimum_required", 2):
        flaws.append(RedTeamFlaw(
            category="mechanism",
            description="Too few pocket contacts; binding is unlikely to be specific.",
            severity=7, substantive=True,
        ))

    # Synthesisability (via the Layer-4 SA score passed through the dossier)
    sa = dossier.get("sa_score")
    if isinstance(sa, (int, float)) and sa > 6.0:
        flaws.append(RedTeamFlaw(
            category="synthesis",
            description=f"SA score {sa:.2f} is too hard to synthesise (>6).",
            severity=7, substantive=True,
        ))

    # Baseline pessimism: honour the "find at least 3 flaws" rule by adding
    # generic (non-substantive) wet-lab-reality reminders until we hit 3.
    baseline_flaws = [
        RedTeamFlaw(
            category="other",
            description=(
                "Computational scores correlate weakly with wet-lab activity; "
                "at least one orthogonal assay (MIC vs. S. aureus clinical "
                "isolates) is needed before calling this a hit."
            ),
            severity=4, substantive=False,
        ),
        RedTeamFlaw(
            category="other",
            description=(
                "Absence of wet-lab selectivity data (off-target hERG, CYP, "
                "mammalian cell tox) means safety profile is assumed, not known."
            ),
            severity=4, substantive=False,
        ),
        RedTeamFlaw(
            category="other",
            description=(
                "Docking scoring functions systematically over-estimate binding "
                "for larger flexible ligands; a consensus run with DiffDock/GNINA "
                "is required before trusting the Vina rank."
            ),
            severity=4, substantive=False,
        ),
    ]
    for bf in baseline_flaws:
        if len(flaws) < 3:
            flaws.append(bf)
        else:
            break

    reasoning = (
        "Deterministic dry-run critique. Live-mode critique is provided by "
        "the adversarial Claude call and replaces this output."
    )
    if notes:
        reasoning += " Notes: " + " | ".join(notes)
    return flaws, reasoning


# ----------------------------------------------------------------------
# Parsing live LLM JSON output
# ----------------------------------------------------------------------


def parse_redteam_json(text: str) -> tuple[list[RedTeamFlaw], str]:
    """Parse a JSON blob from an LLM response. Tolerates leading/trailing prose."""
    start = text.find("{")
    end = text.rfind("}")
    if start < 0 or end < 0 or end < start:
        raise ValueError("no JSON object found in red-team response")
    payload = json.loads(text[start : end + 1])
    raw_flaws = payload.get("flaws") or []
    reasoning = str(payload.get("reasoning") or "")
    flaws = [
        RedTeamFlaw(
            category=str(f.get("category", "other")),
            description=str(f.get("description", "")),
            severity=int(f.get("severity", 5)),
            substantive=bool(f.get("substantive", False)),
        )
        for f in raw_flaws
    ]
    return flaws, reasoning


# ----------------------------------------------------------------------
# Live Claude call (optional)
# ----------------------------------------------------------------------


def _call_claude_adversarial(
    dossier: dict[str, Any],
    model: str = "claude-opus-4-7",
) -> tuple[list[RedTeamFlaw], str, str]:  # pragma: no cover - network
    try:
        import anthropic  # type: ignore
    except ImportError as exc:
        raise RuntimeError(
            "`pip install anthropic` required for live red-team mode"
        ) from exc

    api_key = os.environ.get("ANTHROPIC_API_KEY")
    if not api_key:
        raise RuntimeError("ANTHROPIC_API_KEY not set")

    client = anthropic.Anthropic(api_key=api_key)
    msg = client.messages.create(
        model=model,
        max_tokens=2000,
        system=ADVERSARIAL_SYSTEM,
        messages=[{
            "role": "user",
            "content": USER_TEMPLATE.format(dossier_json=json.dumps(dossier, indent=2)),
        }],
    )
    text = "".join(block.text for block in msg.content if hasattr(block, "text"))
    flaws, reasoning = parse_redteam_json(text)
    return flaws, reasoning, text


# ----------------------------------------------------------------------
# Top-level critique
# ----------------------------------------------------------------------


def _dossier_hash(dossier: dict[str, Any]) -> str:
    import hashlib
    return hashlib.sha256(
        json.dumps(dossier, sort_keys=True).encode()
    ).hexdigest()[:16]


def red_team(
    dossier: dict[str, Any],
    *,
    live: bool = False,
    model: str = "claude-opus-4-7",
    write_dir: Path | None = None,
) -> RedTeamResult:
    candidate_id = str(dossier.get("candidate_id", "unknown"))
    dh = _dossier_hash(dossier)

    result = RedTeamResult(
        candidate_id=candidate_id,
        dossier_hash=dh,
        mode="live" if live else "dry-run",
        model=model if live else None,
        created_at=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    )

    try:
        if live:
            flaws, reasoning, _raw = _call_claude_adversarial(dossier, model=model)
        else:
            flaws, reasoning = _dry_run_flaws(dossier)
    except Exception as exc:
        result.error = str(exc)
        return result

    result.flaws = flaws
    result.reasoning = reasoning
    result.substantive_flaw = any(f.substantive for f in flaws)
    result.severity_score = max((f.severity for f in flaws), default=0)
    result.passed = not result.substantive_flaw

    # Persist
    if write_dir is None:
        write_dir = REPO_ROOT / "findings" / "red-team"
    write_dir.mkdir(parents=True, exist_ok=True)
    json_path = write_dir / f"{candidate_id}.json"
    md_path = write_dir / f"{candidate_id}.md"
    json_path.write_text(json.dumps(asdict(result), indent=2))
    md_path.write_text(_format_markdown(result))

    return result


def _format_markdown(r: RedTeamResult) -> str:
    header = (
        f"# Red-Team Critique — {r.candidate_id}\n\n"
        f"- **Mode:** {r.mode}\n"
        f"- **Model:** {r.model or '—'}\n"
        f"- **Created:** {r.created_at}\n"
        f"- **Dossier hash:** `{r.dossier_hash}`\n"
        f"- **Substantive flaw:** {'YES' if r.substantive_flaw else 'no'}\n"
        f"- **Max severity:** {r.severity_score}/10\n"
        f"- **Verdict:** {'PASS' if r.passed else 'FAIL'}\n\n"
    )
    lines = [header, "## Flaws\n"]
    for i, f in enumerate(r.flaws, 1):
        marker = "⚠️" if f.substantive else "·"
        lines.append(
            f"{i}. {marker} **[{f.category}]** (severity {f.severity}) {f.description}\n"
        )
    lines.append("\n## Reasoning\n")
    lines.append(r.reasoning)
    return "".join(lines) + "\n"


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------


def _main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dossier", type=Path, required=True)
    ap.add_argument("--live", action="store_true", help="call Claude API")
    ap.add_argument("--model", default="claude-opus-4-7")
    ap.add_argument("--out-dir", type=Path, default=None)
    args = ap.parse_args()

    try:
        dossier = json.loads(args.dossier.read_text())
    except Exception as exc:
        print(f"dossier load error: {exc}", file=sys.stderr)
        return 2

    result = red_team(
        dossier, live=args.live, model=args.model, write_dir=args.out_dir
    )
    print(json.dumps(asdict(result), indent=2))
    if result.error:
        return 2
    return 0 if result.passed else 1


if __name__ == "__main__":
    raise SystemExit(_main())
