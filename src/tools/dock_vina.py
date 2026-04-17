#!/usr/bin/env python3
"""
Layer 5 of the verification chain: primary docking via AutoDock Vina.

Responsibilities:
  - Load a target YAML (e.g. targets/mrsa-gyrb.yaml).
  - Download + cache the preferred receptor PDB from RCSB.
  - Compute the docking grid box from the pocket-residue centroid listed
    in the YAML (`binding_pocket.grid_box.mode = centroid_of_residues`).
  - Prepare the candidate SMILES → 3D ligand PDBQT (via RDKit + Meeko when
    available; falls back to obabel CLI if installed).
  - Invoke the `vina` binary with the computed grid box.
  - Parse the top pose energies and return a structured result.
  - Benchmark mode: re-dock all reference inhibitors listed in the target
    YAML; used by loops/positive_control.py hourly.

External binaries expected at runtime on the Orb deploy:
  - vina                          (AutoDock Vina ≥ 1.2)
  - prepare_receptor or meeko     (PDB → PDBQT)
  - obabel (optional fallback)

CLI:
  python -m tools.dock_vina --target targets/mrsa-gyrb.yaml \
      --smiles "CCOc1..." --candidate-id cand-0001
  python -m tools.dock_vina --target targets/mrsa-gyrb.yaml --benchmark

Exit codes: 0 on pass (ΔG ≤ threshold), 1 on fail, 2 on setup error.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import yaml


REPO_ROOT = Path(__file__).parent.parent.parent
CACHE_DIR = REPO_ROOT / "data" / "pdb-cache"
CACHE_DIR.mkdir(parents=True, exist_ok=True)

RCSB_PDB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"


# ----------------------------------------------------------------------
# Data classes
# ----------------------------------------------------------------------


@dataclass
class GridBox:
    center_x: float
    center_y: float
    center_z: float
    size_x: float
    size_y: float
    size_z: float

    def as_vina_args(self) -> list[str]:
        return [
            "--center_x", f"{self.center_x:.3f}",
            "--center_y", f"{self.center_y:.3f}",
            "--center_z", f"{self.center_z:.3f}",
            "--size_x",   f"{self.size_x:.3f}",
            "--size_y",   f"{self.size_y:.3f}",
            "--size_z",   f"{self.size_z:.3f}",
        ]


@dataclass
class DockingResult:
    candidate_id: str
    smiles: str
    target_id: str
    receptor_pdb: str
    engine: str
    best_energy_kcalmol: float | None
    mode_energies: list[float] = field(default_factory=list)
    grid_box: dict[str, float] = field(default_factory=dict)
    threshold_kcalmol: float | None = None
    passed: bool = False
    pose_pdb_path: str | None = None
    duration_s: float = 0.0
    error: str | None = None


# ----------------------------------------------------------------------
# PDB download + parsing
# ----------------------------------------------------------------------


def fetch_pdb(pdb_id: str, cache_dir: Path = CACHE_DIR) -> Path:
    """Download a PDB file from RCSB if not cached. Returns local path."""
    pdb_id = pdb_id.upper()
    local = cache_dir / f"{pdb_id}.pdb"
    if local.exists() and local.stat().st_size > 1000:
        return local
    try:
        import urllib.request
        url = RCSB_PDB_URL.format(pdb_id=pdb_id)
        with urllib.request.urlopen(url, timeout=60) as resp:
            data = resp.read()
        if len(data) < 1000 or not data.startswith((b"HEADER", b"ATOM", b"REMARK")):
            raise RuntimeError(f"RCSB returned unexpected content for {pdb_id}")
        local.write_bytes(data)
        return local
    except Exception as exc:
        raise RuntimeError(f"failed to download PDB {pdb_id}: {exc}") from exc


def parse_pdb_residue_centroids(
    pdb_path: Path,
    chain: str,
    residues: list[int],
) -> dict[int, tuple[float, float, float]]:
    """Return a map resnum → (x, y, z) centroid of the residue's atoms.

    Reads standard ATOM records only.
    """
    centroids: dict[int, list[tuple[float, float, float]]] = {r: [] for r in residues}
    target_residues = set(residues)

    with pdb_path.open() as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            rec_chain = line[21].strip() or " "
            try:
                resnum = int(line[22:26].strip())
            except ValueError:
                continue
            if chain and rec_chain != chain:
                continue
            if resnum not in target_residues:
                continue
            try:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            except ValueError:
                continue
            centroids[resnum].append((x, y, z))

    result: dict[int, tuple[float, float, float]] = {}
    for r, coords in centroids.items():
        if not coords:
            continue
        n = len(coords)
        cx = sum(c[0] for c in coords) / n
        cy = sum(c[1] for c in coords) / n
        cz = sum(c[2] for c in coords) / n
        result[r] = (cx, cy, cz)
    return result


def compute_grid_box(
    pdb_path: Path,
    chain: str,
    pocket_residues: list[int],
    padding_a: float,
    side_length_a: float,
) -> GridBox:
    """Compute Vina grid box centered on the centroid of the pocket residues.

    Size is fixed (side_length_a) regardless of padding; padding_a is used
    when falling back to pocket extent mode. Keeping a fixed cube is the
    safest default for a single binding pocket.
    """
    centroids = parse_pdb_residue_centroids(pdb_path, chain, pocket_residues)
    if not centroids:
        raise RuntimeError(
            f"No pocket residues {pocket_residues} found in chain {chain} "
            f"of {pdb_path.name}; check YAML residue numbering vs PDB."
        )
    n = len(centroids)
    cx = sum(c[0] for c in centroids.values()) / n
    cy = sum(c[1] for c in centroids.values()) / n
    cz = sum(c[2] for c in centroids.values()) / n
    return GridBox(
        center_x=cx, center_y=cy, center_z=cz,
        size_x=side_length_a, size_y=side_length_a, size_z=side_length_a,
    )


# ----------------------------------------------------------------------
# Ligand preparation
# ----------------------------------------------------------------------


def smiles_to_pdbqt(smiles: str, out_dir: Path, name: str = "ligand") -> Path:
    """Prepare a SMILES string as an AutoDock PDBQT ligand.

    Preference order:
      1. meeko (Python) — the canonical modern path.
      2. obabel CLI — fallback.
    Raises RuntimeError if neither is available.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")
    mol = Chem.AddHs(mol)
    embed_code = AllChem.EmbedMolecule(mol, randomSeed=42)
    if embed_code != 0:
        # Retry with a different random seed before giving up
        embed_code = AllChem.EmbedMolecule(mol, randomSeed=1729)
        if embed_code != 0:
            raise RuntimeError(f"3D embedding failed for SMILES: {smiles!r}")
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:
        # Fall through; partial charges still OK for Vina.
        pass

    out_dir.mkdir(parents=True, exist_ok=True)
    pdb_path = out_dir / f"{name}.pdb"
    pdbqt_path = out_dir / f"{name}.pdbqt"
    Chem.MolToPDBFile(mol, str(pdb_path))

    # Path 1: Meeko
    try:
        from meeko import MoleculePreparation, PDBQTWriterLegacy  # type: ignore
        prep = MoleculePreparation()
        prep.prepare(mol)
        pdbqt_str, _, _ = PDBQTWriterLegacy.write_string(prep.setup)
        pdbqt_path.write_text(pdbqt_str)
        return pdbqt_path
    except Exception:
        pass

    # Path 2: obabel CLI
    obabel = shutil.which("obabel")
    if obabel is not None:
        subprocess.run(
            [obabel, str(pdb_path), "-O", str(pdbqt_path), "--gen3d"],
            check=True, capture_output=True,
        )
        return pdbqt_path

    raise RuntimeError(
        "No ligand-prep tool available. `pip install meeko` or install "
        "OpenBabel (obabel) on the system PATH."
    )


def prepare_receptor_pdbqt(pdb_path: Path, chain: str, out_dir: Path) -> Path:
    """Prepare a receptor PDBQT. Uses ADFR `prepare_receptor` if available,
    otherwise obabel. Removes waters and cofactors (HETATM) by default."""
    out_dir.mkdir(parents=True, exist_ok=True)
    clean_pdb = out_dir / f"{pdb_path.stem}_clean.pdb"

    # Keep only ATOM records of the requested chain (strip HETATM: waters,
    # cofactors, co-crystal ligands). Keeps TER for chain termination.
    with pdb_path.open() as src, clean_pdb.open("w") as dst:
        for line in src:
            if line.startswith("ATOM"):
                rec_chain = line[21].strip() or " "
                if chain and rec_chain != chain:
                    continue
                dst.write(line)
            elif line.startswith("TER"):
                dst.write(line)
            elif line.startswith("END"):
                dst.write(line)

    pdbqt_path = out_dir / f"{pdb_path.stem}_{chain}.pdbqt"

    prepare_receptor = shutil.which("prepare_receptor")
    if prepare_receptor is not None:
        subprocess.run(
            [prepare_receptor, "-r", str(clean_pdb), "-o", str(pdbqt_path)],
            check=True, capture_output=True,
        )
        return pdbqt_path

    obabel = shutil.which("obabel")
    if obabel is not None:
        subprocess.run(
            [obabel, str(clean_pdb), "-O", str(pdbqt_path), "-xr"],
            check=True, capture_output=True,
        )
        return pdbqt_path

    raise RuntimeError(
        "No receptor-prep tool available. Install AutoDockFR "
        "(`prepare_receptor`) or OpenBabel (`obabel`)."
    )


# ----------------------------------------------------------------------
# Vina invocation + output parsing
# ----------------------------------------------------------------------


VINA_MODE_LINE = re.compile(
    r"^\s*(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$"
)


def parse_vina_stdout(stdout: str) -> list[float]:
    """Extract the binding-energy column (kcal/mol) from a Vina stdout."""
    energies: list[float] = []
    in_table = False
    for raw in stdout.splitlines():
        line = raw.rstrip()
        if "mode |   affinity" in line or line.strip().startswith("mode"):
            in_table = True
            continue
        if not in_table:
            continue
        if line.strip().startswith("-----"):
            continue
        m = VINA_MODE_LINE.match(line)
        if not m:
            if energies:
                break
            continue
        energies.append(float(m.group(2)))
    return energies


def run_vina(
    *,
    receptor_pdbqt: Path,
    ligand_pdbqt: Path,
    grid: GridBox,
    out_pdbqt: Path,
    exhaustiveness: int,
    num_modes: int,
    seed: int,
    cpu: int | None = None,
) -> tuple[list[float], str]:
    vina = shutil.which("vina")
    if vina is None:
        raise RuntimeError(
            "AutoDock Vina binary not on PATH. Install from "
            "https://vina.scripps.edu/ or `conda install -c bioconda autodock-vina`."
        )
    args = [
        vina,
        "--receptor", str(receptor_pdbqt),
        "--ligand",   str(ligand_pdbqt),
        "--out",      str(out_pdbqt),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
        "--seed", str(seed),
        *grid.as_vina_args(),
    ]
    if cpu is not None:
        args += ["--cpu", str(cpu)]
    proc = subprocess.run(args, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"vina failed (rc={proc.returncode}):\nstderr={proc.stderr}\n"
            f"stdout={proc.stdout}"
        )
    return parse_vina_stdout(proc.stdout), proc.stdout


# ----------------------------------------------------------------------
# Top-level docker
# ----------------------------------------------------------------------


class VinaDocker:
    def __init__(self, target_yaml: Path, *, work_dir: Path | None = None):
        self.target_path = Path(target_yaml)
        self.target = yaml.safe_load(self.target_path.read_text())
        self.target_id: str = self.target["id"]
        self.work_dir = Path(work_dir) if work_dir else REPO_ROOT / "findings" / "docking"
        self.work_dir.mkdir(parents=True, exist_ok=True)

        # Pick preferred receptor
        preferred = next(
            (p for p in self.target["structure"]["pdb_ids"]
             if p.get("preferred_docking_receptor")),
            self.target["structure"]["pdb_ids"][0],
        )
        self.pdb_id: str = preferred["id"]
        self.chain: str = self.target["structure"].get("chain_for_docking", "A")

        self.pocket_residues: list[int] = [
            r["resnum"] for r in self.target["binding_pocket"]["pocket_residues"]
        ]
        gb_cfg = self.target["binding_pocket"]["grid_box"]
        self.grid_padding_a: float = float(gb_cfg.get("padding_a", 8.0))
        self.grid_side_a: float = float(gb_cfg.get("side_length_a", 22.0))

        dock_cfg = self.target["docking"]["primary"]
        self.exhaustiveness: int = int(dock_cfg.get("exhaustiveness", 32))
        self.num_modes: int = int(dock_cfg.get("num_modes", 20))
        self.seed: int = int(dock_cfg.get("seed", 42))
        self.threshold_kcalmol: float = float(
            dock_cfg.get("delta_g_threshold_kcalmol", -8.0)
        )

        self._receptor_pdb: Path | None = None
        self._receptor_pdbqt: Path | None = None
        self._grid_box: GridBox | None = None

    # --- Preparation -------------------------------------------------

    def receptor_pdb(self) -> Path:
        if self._receptor_pdb is None:
            self._receptor_pdb = fetch_pdb(self.pdb_id)
        return self._receptor_pdb

    def grid_box(self) -> GridBox:
        if self._grid_box is None:
            self._grid_box = compute_grid_box(
                self.receptor_pdb(),
                self.chain,
                self.pocket_residues,
                self.grid_padding_a,
                self.grid_side_a,
            )
        return self._grid_box

    def receptor_pdbqt(self) -> Path:
        if self._receptor_pdbqt is None:
            self._receptor_pdbqt = prepare_receptor_pdbqt(
                self.receptor_pdb(), self.chain, self.work_dir / "receptors"
            )
        return self._receptor_pdbqt

    # --- Docking -----------------------------------------------------

    def dock(self, smiles: str, candidate_id: str) -> DockingResult:
        t0 = time.time()
        grid = self.grid_box()
        result = DockingResult(
            candidate_id=candidate_id,
            smiles=smiles,
            target_id=self.target_id,
            receptor_pdb=self.pdb_id,
            engine="autodock_vina",
            best_energy_kcalmol=None,
            grid_box={
                "center_x": grid.center_x, "center_y": grid.center_y,
                "center_z": grid.center_z, "size_x": grid.size_x,
                "size_y": grid.size_y, "size_z": grid.size_z,
            },
            threshold_kcalmol=self.threshold_kcalmol,
        )
        cand_dir = self.work_dir / candidate_id
        cand_dir.mkdir(parents=True, exist_ok=True)
        try:
            ligand_pdbqt = smiles_to_pdbqt(smiles, cand_dir, name="ligand")
            out_pdbqt = cand_dir / "poses.pdbqt"
            energies, stdout = run_vina(
                receptor_pdbqt=self.receptor_pdbqt(),
                ligand_pdbqt=ligand_pdbqt,
                grid=grid,
                out_pdbqt=out_pdbqt,
                exhaustiveness=self.exhaustiveness,
                num_modes=self.num_modes,
                seed=self.seed,
            )
            result.mode_energies = energies
            result.best_energy_kcalmol = min(energies) if energies else None
            result.pose_pdb_path = str(out_pdbqt)
            (cand_dir / "vina.stdout.txt").write_text(stdout)
            if result.best_energy_kcalmol is not None:
                result.passed = result.best_energy_kcalmol <= self.threshold_kcalmol
        except Exception as exc:
            result.error = str(exc)
        result.duration_s = round(time.time() - t0, 2)
        (cand_dir / "docking.json").write_text(json.dumps(asdict(result), indent=2))
        return result

    def benchmark(self) -> list[DockingResult]:
        """Dock every reference inhibitor in the target YAML (positive controls)."""
        refs = self.target.get("reference_inhibitors", [])
        results: list[DockingResult] = []
        for ref in refs:
            smiles = ref.get("smiles", "").strip()
            if not smiles:
                continue
            rid = f"bench-{ref['name'].replace(' ', '_').lower()}"
            results.append(self.dock(smiles, rid))
        return results


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------


def _main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--target", type=Path, required=True, help="target YAML")
    ap.add_argument("--smiles", help="candidate SMILES to dock")
    ap.add_argument("--candidate-id", default=None, help="candidate ID (default: derived)")
    ap.add_argument("--benchmark", action="store_true",
                    help="dock every reference inhibitor in the target YAML")
    ap.add_argument("--work-dir", type=Path, default=None,
                    help="output dir (default: findings/docking)")
    args = ap.parse_args()

    try:
        docker = VinaDocker(args.target, work_dir=args.work_dir)
    except Exception as exc:
        print(f"setup error: {exc}", file=sys.stderr)
        return 2

    if args.benchmark:
        results = docker.benchmark()
        print(json.dumps([asdict(r) for r in results], indent=2))
        return 0 if all(r.passed for r in results if r.best_energy_kcalmol is not None) else 1

    if not args.smiles:
        ap.error("--smiles required unless --benchmark")

    cid = args.candidate_id or f"cand-{int(time.time())}"
    result = docker.dock(args.smiles, cid)
    print(json.dumps(asdict(result), indent=2))
    if result.error:
        return 2
    return 0 if result.passed else 1


if __name__ == "__main__":
    raise SystemExit(_main())
