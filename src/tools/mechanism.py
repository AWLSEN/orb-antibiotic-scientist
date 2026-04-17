#!/usr/bin/env python3
"""
Layer 7: mechanism sanity via protein-ligand interaction analysis.

A candidate passes Layer 7 only when its top Vina pose makes the right kind
of contacts with the binding pocket:

  1. At least `binding_pocket.required_contacts.minimum` contacts with
     pocket residues listed in the target YAML.
  2. At least one of those contacts is with a residue in
     `binding_pocket.required_contacts.must_include_any_of`.
  3. The candidate does NOT contact only residues listed as
     `resistance_hotspots.single_site`. Candidates binding exclusively
     to evolvable residues are demoted (resistance will evolve quickly).

Contact extraction backend:

  - Preferred: PLIP (Protein-Ligand Interaction Profiler). When `plip` is
    importable, we run the full interaction profile (H-bonds, hydrophobic,
    salt bridges, π-stacking, halogen bonds).
  - Fallback: a simple residue-centroid distance analyzer that flags any
    pocket residue whose centroid is within `contact_cutoff_a` Å of any
    ligand heavy atom. Coarser but PLIP-independent, which keeps CI green
    on hosts without PLIP / openbabel installed.

The PLIP backend runs on a combined protein + ligand PDB file produced
by dock_vina.py (receptor + best pose). The fallback reads the Vina
output pose directly.

CLI:
  python -m tools.mechanism \
      --target targets/mrsa-gyrb.yaml \
      --complex-pdb findings/docking/<cid>/complex.pdb \
      --candidate-id cand-0001
  python -m tools.mechanism --check-gate --contacts 46,73,76 --target ...

Exit codes: 0 pass, 1 fail, 2 setup error.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Iterable

import yaml


# ----------------------------------------------------------------------
# Data classes
# ----------------------------------------------------------------------


@dataclass
class Contact:
    residue_num: int
    residue_name: str
    chain: str
    interaction_type: str
    distance_a: float | None = None


@dataclass
class MechanismResult:
    candidate_id: str
    target_id: str
    pocket_residues_hit: list[int] = field(default_factory=list)
    catalytic_residues_hit: list[int] = field(default_factory=list)
    must_include_hits: list[int] = field(default_factory=list)
    resistance_hotspots_hit: list[int] = field(default_factory=list)
    contacts: list[dict[str, Any]] = field(default_factory=list)
    num_contacts: int = 0
    minimum_required: int = 2
    resistance_only: bool = False
    passed: bool = False
    backend: str = "unknown"
    reasons: list[str] = field(default_factory=list)
    error: str | None = None


# ----------------------------------------------------------------------
# Gate logic — depends only on the target YAML and a contact set.
# ----------------------------------------------------------------------


def check_gate(
    contact_residues: Iterable[int],
    target: dict[str, Any],
) -> dict[str, Any]:
    """Pure function: given a set of residue numbers the ligand contacts
    and a loaded target YAML, return a dict with pocket / catalytic /
    must-include hits, resistance-hotspot hits, resistance_only boolean,
    and passed boolean with reasons.
    """
    contact_set = set(int(r) for r in contact_residues)

    pocket_res = {int(r["resnum"]) for r in target["binding_pocket"]["pocket_residues"]}
    catalytic = set(target["binding_pocket"].get("catalytic_residues", []))
    req_cfg = target["binding_pocket"]["required_contacts"]
    minimum = int(req_cfg.get("minimum", 2))
    must_include_any_of = set(req_cfg.get("must_include_any_of", []))
    resistance_single = {
        int(r["resnum"])
        for r in target.get("resistance_hotspots", {}).get("single_site", [])
    }

    pocket_hits = sorted(contact_set & pocket_res)
    catalytic_hits = sorted(contact_set & catalytic)
    must_hits = sorted(contact_set & must_include_any_of)
    resistance_hits = sorted(contact_set & resistance_single)

    reasons: list[str] = []

    if len(contact_set) < minimum:
        reasons.append(
            f"too few contacts: {len(contact_set)} < required {minimum}"
        )

    if must_include_any_of and not must_hits:
        reasons.append(
            f"no contact with any of must_include_any_of={sorted(must_include_any_of)}"
        )

    # "Resistance-only" = the only residues the ligand hits are in
    # resistance_hotspots (evolvable). A candidate that also hits at least
    # one non-hotspot is fine.
    non_hotspot_pocket_hits = [r for r in pocket_hits if r not in resistance_hits]
    resistance_only = bool(resistance_hits) and not non_hotspot_pocket_hits
    if resistance_only:
        reasons.append(
            f"candidate binds only resistance hotspots {resistance_hits} — "
            "resistance would evolve quickly"
        )

    passed = not reasons
    return {
        "pocket_residues_hit": pocket_hits,
        "catalytic_residues_hit": catalytic_hits,
        "must_include_hits": must_hits,
        "resistance_hotspots_hit": resistance_hits,
        "resistance_only": resistance_only,
        "minimum_required": minimum,
        "reasons": reasons,
        "passed": passed,
    }


# ----------------------------------------------------------------------
# Backend: PLIP (preferred)
# ----------------------------------------------------------------------


def _plip_contact_residues(complex_pdb: Path, chain: str) -> list[Contact]:
    """Run PLIP and return a list of Contact objects."""
    from plip.structure.preparation import PDBComplex  # type: ignore
    comp = PDBComplex()
    comp.load_pdb(str(complex_pdb))
    comp.analyze()

    contacts: list[Contact] = []
    for _bsid, interaction in comp.interaction_sets.items():
        for h in getattr(interaction, "hbonds_ldon", []) + \
                 getattr(interaction, "hbonds_pdon", []):
            contacts.append(Contact(
                residue_num=h.resnr, residue_name=h.restype,
                chain=h.reschain, interaction_type="hbond",
                distance_a=round(float(h.distance_ad), 2),
            ))
        for hc in getattr(interaction, "hydrophobic_contacts", []):
            contacts.append(Contact(
                residue_num=hc.resnr, residue_name=hc.restype,
                chain=hc.reschain, interaction_type="hydrophobic",
                distance_a=round(float(hc.distance), 2),
            ))
        for sb in getattr(interaction, "saltbridge_lneg", []) + \
                  getattr(interaction, "saltbridge_pneg", []):
            contacts.append(Contact(
                residue_num=sb.resnr, residue_name=sb.restype,
                chain=sb.reschain, interaction_type="salt_bridge",
                distance_a=round(float(sb.distance), 2),
            ))
        for ps in getattr(interaction, "pistacking", []):
            contacts.append(Contact(
                residue_num=ps.resnr, residue_name=ps.restype,
                chain=ps.reschain, interaction_type="pi_stacking",
                distance_a=round(float(ps.distance), 2),
            ))
        for hx in getattr(interaction, "halogen_bonds", []):
            contacts.append(Contact(
                residue_num=hx.resnr, residue_name=hx.restype,
                chain=hx.reschain, interaction_type="halogen",
                distance_a=round(float(hx.distance_az), 2),
            ))

    if chain:
        contacts = [c for c in contacts if not c.chain or c.chain == chain]
    return contacts


# ----------------------------------------------------------------------
# Backend: simple distance (fallback)
# ----------------------------------------------------------------------


def _parse_pdb_atoms(pdb_path: Path):
    """Yield (record_type, chain, resnum, resname, atom_name, x, y, z)."""
    with pdb_path.open() as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec not in ("ATOM", "HETATM"):
                continue
            chain = line[21].strip() or " "
            try:
                resnum = int(line[22:26])
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            except ValueError:
                continue
            atom_name = line[12:16].strip()
            resname = line[17:20].strip()
            yield rec, chain, resnum, resname, atom_name, x, y, z


def _distance_backend_contacts(
    complex_pdb: Path,
    chain: str,
    pocket_residues: list[int],
    cutoff_a: float = 4.5,
) -> list[Contact]:
    """Simple backend: flag any pocket residue whose any atom is within
    `cutoff_a` Å of any HETATM ligand atom. Interaction type is reported
    as `proximity`. Coarse but always available."""
    pocket_set = set(pocket_residues)
    protein_atoms: list[tuple[int, str, str, float, float, float]] = []
    ligand_atoms: list[tuple[float, float, float]] = []

    for rec, ch, resnum, resname, _aname, x, y, z in _parse_pdb_atoms(complex_pdb):
        if rec == "ATOM" and (not chain or ch == chain) and resnum in pocket_set:
            protein_atoms.append((resnum, resname, ch, x, y, z))
        elif rec == "HETATM":
            ligand_atoms.append((x, y, z))

    contacts: list[Contact] = []
    seen: dict[int, float] = {}
    for resnum, resname, ch, px, py, pz in protein_atoms:
        best = None
        for lx, ly, lz in ligand_atoms:
            d2 = (px - lx) ** 2 + (py - ly) ** 2 + (pz - lz) ** 2
            if d2 <= cutoff_a * cutoff_a:
                d = d2 ** 0.5
                if best is None or d < best:
                    best = d
        if best is not None:
            prev = seen.get(resnum)
            if prev is None or best < prev:
                seen[resnum] = best
                contacts.append(Contact(
                    residue_num=resnum, residue_name=resname, chain=ch,
                    interaction_type="proximity", distance_a=round(best, 2),
                ))
    return contacts


# ----------------------------------------------------------------------
# Public API
# ----------------------------------------------------------------------


def analyze_pose(
    complex_pdb: Path,
    target_yaml: Path,
    candidate_id: str,
    *,
    prefer_plip: bool = True,
    cutoff_a: float = 4.5,
) -> MechanismResult:
    target = yaml.safe_load(Path(target_yaml).read_text())
    chain = target["structure"].get("chain_for_docking", "A")
    pocket_residues = [
        int(r["resnum"]) for r in target["binding_pocket"]["pocket_residues"]
    ]

    backend = "distance"
    contacts: list[Contact] = []
    error: str | None = None

    if prefer_plip:
        try:
            contacts = _plip_contact_residues(complex_pdb, chain)
            backend = "plip"
        except ImportError:
            contacts = _distance_backend_contacts(
                complex_pdb, chain, pocket_residues, cutoff_a=cutoff_a,
            )
        except Exception as exc:
            error = f"PLIP failed, falling back to distance backend: {exc}"
            contacts = _distance_backend_contacts(
                complex_pdb, chain, pocket_residues, cutoff_a=cutoff_a,
            )
    else:
        contacts = _distance_backend_contacts(
            complex_pdb, chain, pocket_residues, cutoff_a=cutoff_a,
        )

    contact_residues = {c.residue_num for c in contacts}
    gate = check_gate(contact_residues, target)

    return MechanismResult(
        candidate_id=candidate_id,
        target_id=target.get("id", ""),
        pocket_residues_hit=gate["pocket_residues_hit"],
        catalytic_residues_hit=gate["catalytic_residues_hit"],
        must_include_hits=gate["must_include_hits"],
        resistance_hotspots_hit=gate["resistance_hotspots_hit"],
        contacts=[asdict(c) for c in contacts],
        num_contacts=len(contacts),
        minimum_required=gate["minimum_required"],
        resistance_only=gate["resistance_only"],
        passed=gate["passed"],
        backend=backend,
        reasons=gate["reasons"],
        error=error,
    )


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------


def _main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--target", type=Path, required=True)
    ap.add_argument("--complex-pdb", type=Path, help="combined protein+ligand PDB")
    ap.add_argument("--candidate-id", default="cand-xxxx")
    ap.add_argument(
        "--check-gate",
        action="store_true",
        help="bypass PDB parsing; --contacts supplies the residue list",
    )
    ap.add_argument(
        "--contacts",
        type=str,
        default="",
        help="comma-separated residue numbers for --check-gate mode",
    )
    ap.add_argument("--cutoff", type=float, default=4.5)
    ap.add_argument("--no-plip", action="store_true")
    args = ap.parse_args()

    try:
        target = yaml.safe_load(args.target.read_text())
    except Exception as exc:
        print(f"setup error: {exc}", file=sys.stderr)
        return 2

    if args.check_gate:
        if not args.contacts:
            print("--check-gate needs --contacts=73,50,46", file=sys.stderr)
            return 2
        residues = [int(x) for x in args.contacts.split(",") if x.strip()]
        gate = check_gate(residues, target)
        print(json.dumps(gate, indent=2))
        return 0 if gate["passed"] else 1

    if not args.complex_pdb:
        print("--complex-pdb required unless --check-gate", file=sys.stderr)
        return 2

    result = analyze_pose(
        args.complex_pdb,
        args.target,
        args.candidate_id,
        prefer_plip=not args.no_plip,
        cutoff_a=args.cutoff,
    )
    print(json.dumps(asdict(result), indent=2))
    return 0 if result.passed else 1


if __name__ == "__main__":
    raise SystemExit(_main())
