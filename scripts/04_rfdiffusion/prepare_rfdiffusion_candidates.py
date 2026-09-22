#!/usr/bin/env python3
"""
Convert screened C3 HA-scaffold candidate CIFs to RFdiffusion-ready PDBs.

The input candidates are combined HA+scaffold mmCIF files from
align/filter_linker_candidates.py. This script writes PDB copies under a
repo-level folder and, by default, canonicalizes the fitted C3 axis to the
RFdiffusion cyclic symmetry Z axis.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable

import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
RIGID_DIR = SCRIPT_DIR.parent
REPO_DIR = RIGID_DIR.parent
ALIGN_DIR = RIGID_DIR / "align"


@dataclass(frozen=True)
class Atom:
    group: str
    atom_id: int
    type_symbol: str
    atom_name: str
    alt_id: str
    res_name: str
    label_asym_id: str
    label_entity_id: str
    label_seq_id: int
    ins_code: str
    x: float
    y: float
    z: float
    occupancy: str
    b_iso: str
    formal_charge: str
    chain_id: str
    model_num: str

    @property
    def coord(self) -> np.ndarray:
        return np.array([self.x, self.y, self.z], dtype=float)


def parse_mmcif_atoms(path: Path) -> list[Atom]:
    atoms: list[Atom] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            parts = line.strip().split()
            if len(parts) < 18:
                raise ValueError(f"Cannot parse atom line in {path}: {line.rstrip()}")
            atoms.append(
                Atom(
                    group=parts[0],
                    atom_id=int(parts[1]),
                    type_symbol=parts[2],
                    atom_name=parts[3],
                    alt_id=parts[4],
                    res_name=parts[5],
                    label_asym_id=parts[6],
                    label_entity_id=parts[7],
                    label_seq_id=int(parts[8]),
                    ins_code=parts[9],
                    x=float(parts[10]),
                    y=float(parts[11]),
                    z=float(parts[12]),
                    occupancy=parts[13],
                    b_iso=parts[14],
                    formal_charge=parts[15],
                    chain_id=parts[16],
                    model_num=parts[17],
                )
            )
    if not atoms:
        raise ValueError(f"No ATOM/HETATM records found in {path}")
    return atoms


def normalize(vector: np.ndarray, label: str) -> np.ndarray:
    norm = float(np.linalg.norm(vector))
    if norm < 1e-8:
        raise ValueError(f"Cannot normalize near-zero {label}")
    return vector / norm


def rotation_about_axis(axis: np.ndarray, angle_rad: float) -> np.ndarray:
    axis = normalize(axis, "rotation axis")
    x, y, z = axis
    c = math.cos(angle_rad)
    s = math.sin(angle_rad)
    one_c = 1.0 - c
    return np.array(
        [
            [c + x * x * one_c, x * y * one_c - z * s, x * z * one_c + y * s],
            [y * x * one_c + z * s, c + y * y * one_c, y * z * one_c - x * s],
            [z * x * one_c - y * s, z * y * one_c + x * s, c + z * z * one_c],
        ]
    )


def rotation_between_vectors(source: np.ndarray, target: np.ndarray) -> np.ndarray:
    source = normalize(source, "source axis")
    target = normalize(target, "target axis")
    cross = np.cross(source, target)
    dot = float(np.dot(source, target))

    if dot > 1.0 - 1e-10:
        return np.eye(3)
    if dot < -1.0 + 1e-10:
        basis = np.array([1.0, 0.0, 0.0])
        if abs(float(np.dot(source, basis))) > 0.8:
            basis = np.array([0.0, 1.0, 0.0])
        axis = normalize(np.cross(source, basis), "180-degree rotation axis")
        return rotation_about_axis(axis, math.pi)

    skew = np.array(
        [
            [0.0, -cross[2], cross[1]],
            [cross[2], 0.0, -cross[0]],
            [-cross[1], cross[0], 0.0],
        ]
    )
    return np.eye(3) + skew + skew @ skew * ((1.0 - dot) / (np.linalg.norm(cross) ** 2))


def ca_by_chain_residue(
    atoms: list[Atom],
    chains: list[str],
    min_seq: int,
    max_seq: int,
) -> dict[str, dict[int, np.ndarray]]:
    ca: dict[str, dict[int, np.ndarray]] = {chain: {} for chain in chains}
    for atom in atoms:
        if (
            atom.chain_id in ca
            and min_seq <= atom.label_seq_id <= max_seq
            and atom.atom_name == "CA"
            and atom.alt_id in (".", "A", "?")
        ):
            ca[atom.chain_id][atom.label_seq_id] = atom.coord
    return ca


def fit_c3_axis(
    atoms: list[Atom],
    chains: list[str],
    min_seq: int,
    max_seq: int,
) -> tuple[np.ndarray, np.ndarray]:
    ca = ca_by_chain_residue(atoms, chains, min_seq, max_seq)
    common_residues = set(ca[chains[0]])
    for chain in chains[1:]:
        common_residues &= set(ca[chain])
    if len(common_residues) < 2:
        raise ValueError(f"Need at least two shared CA residues in {min_seq}-{max_seq}")

    centers = np.array(
        [np.mean([ca[chain][seq] for chain in chains], axis=0) for seq in sorted(common_residues)]
    )
    axis_point = np.mean(centers, axis=0)
    _, _, vh = np.linalg.svd(centers - axis_point, full_matrices=False)
    axis = normalize(vh[0], "C3 axis")
    return axis_point, axis


def atom_coord(atoms: list[Atom], chain: str, seq: int, atom_name: str) -> np.ndarray:
    for atom in atoms:
        if (
            atom.chain_id == chain
            and atom.label_seq_id == seq
            and atom.atom_name == atom_name
            and atom.alt_id in (".", "A", "?")
        ):
            return atom.coord
    raise ValueError(f"Missing {chain}:{seq}:{atom_name}")


def canonicalize_c3_axis(
    atoms: list[Atom],
    chains: list[str],
    ha_min_residue: int,
    ha_end_residue: int,
) -> list[Atom]:
    axis_point, axis = fit_c3_axis(atoms, chains, ha_min_residue, ha_end_residue)
    ha_c_centroid = np.mean(
        [atom_coord(atoms, chain, ha_end_residue, "C") for chain in chains],
        axis=0,
    )
    body_ca = np.array([atom.coord for atom in atoms if atom.atom_name == "CA"])
    body_centroid = np.mean(body_ca, axis=0)
    if float(np.dot(axis, ha_c_centroid - body_centroid)) < 0.0:
        axis = -axis

    rotation = rotation_between_vectors(axis, np.array([0.0, 0.0, 1.0]))
    transformed: list[Atom] = []
    for atom in atoms:
        coord = rotation @ (atom.coord - axis_point)
        transformed.append(replace(atom, x=float(coord[0]), y=float(coord[1]), z=float(coord[2])))
    return transformed


def pdb_atom_name(atom_name: str, element: str) -> str:
    if len(atom_name) >= 4:
        return atom_name[:4]
    if len(element) == 1:
        return f" {atom_name:<3}"
    return f"{atom_name:<4}"


def format_pdb_atom(atom: Atom, serial: int) -> str:
    record = "HETATM" if atom.group == "HETATM" else "ATOM  "
    element = (atom.type_symbol or atom.atom_name[0]).strip().upper()[:2]
    atom_name = pdb_atom_name(atom.atom_name, element)
    alt_id = " " if atom.alt_id in (".", "?") else atom.alt_id[:1]
    ins_code = " " if atom.ins_code in (".", "?") else atom.ins_code[:1]
    try:
        occupancy = float(atom.occupancy)
    except ValueError:
        occupancy = 1.0
    try:
        b_iso = float(atom.b_iso)
    except ValueError:
        b_iso = 0.0
    return (
        f"{record}{serial:5d} {atom_name}{alt_id}{atom.res_name:>3} {atom.chain_id[:1]}"
        f"{atom.label_seq_id:4d}{ins_code}   "
        f"{atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}"
        f"{occupancy:6.2f}{b_iso:6.2f}          {element:>2}\n"
    )


def write_pdb(path: Path, atoms: Iterable[Atom], remarks: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for remark in remarks:
            handle.write(f"REMARK {remark}\n")
        last_atom = None
        serial = 0
        for atom in atoms:
            if last_atom is not None and atom.chain_id != last_atom.chain_id:
                serial += 1
                handle.write(
                    f"TER   {serial:5d}      {last_atom.res_name:>3} "
                    f"{last_atom.chain_id[:1]}{last_atom.label_seq_id:4d}\n"
                )
            serial += 1
            handle.write(format_pdb_atom(atom, serial))
            last_atom = atom
        if last_atom is not None:
            serial += 1
            handle.write(
                f"TER   {serial:5d}      {last_atom.res_name:>3} "
                f"{last_atom.chain_id[:1]}{last_atom.label_seq_id:4d}\n"
            )
        handle.write("END\n")


def parse_rank(path: Path) -> int:
    name = path.stem
    if name.startswith("rank_"):
        try:
            return int(name.split("_", 2)[1])
        except (IndexError, ValueError):
            pass
    return 999999


def contig_for_chains(
    chains: list[str],
    ha_min_residue: int,
    ha_end_residue: int,
    scaffold_start_residue: int,
    scaffold_max_residue: int,
    linker_spec: str,
) -> str:
    parts = [
        f"{chain}{ha_min_residue}-{ha_end_residue}/{linker_spec}/"
        f"{chain}{scaffold_start_residue}-{scaffold_max_residue}"
        for chain in chains
    ]
    return "/0 ".join(parts) + "/0"


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Convert screened C3 candidate CIFs to canonical PDBs for RFdiffusion."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=ALIGN_DIR / "ck_xt_24_linker_candidates" / "cif",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=RIGID_DIR / "rfdiffusion_candidate_pdb_ck_xt_24_20aa",
    )
    parser.add_argument(
        "--rf-inputs-dir",
        type=Path,
        default=REPO_DIR / "inputs" / "ck_xt_24_20aa",
        help="Additional RFdiffusion input directory that receives the generated PDBs.",
    )
    parser.add_argument("--chains", default="A,B,C")
    parser.add_argument(
        "--ha-min-residue",
        type=int,
        default=1,
        help="First HA residue used for C3 axis fitting during canonicalization",
    )
    parser.add_argument(
        "--contig-ha-start",
        type=int,
        default=480,
        help="First HA residue included in the fixed RFdiffusion contig segment",
    )
    parser.add_argument("--ha-end-residue", type=int, default=525)
    parser.add_argument("--scaffold-start-residue", type=int, default=526)
    parser.add_argument("--scaffold-max-residue", type=int, default=652)
    parser.add_argument("--linker-spec", default="20-20", help="RFdiffusion contig length range for each linker")
    parser.add_argument("--num-designs", type=int, default=100, help="RFdiffusion designs per candidate")
    parser.add_argument(
        "--command-input-prefix",
        default="inputs/ck_xt_24_20aa",
        help="Optional POSIX prefix for command drafts, e.g. /inputs/rfdiffusion_candidate_pdb",
    )
    parser.add_argument("--no-canonicalize", action="store_true", help="Keep the current coordinate frame")
    return parser


def main() -> None:
    args = make_parser().parse_args()
    chains = [chain.strip() for chain in args.chains.split(",") if chain.strip()]
    if len(chains) != 3:
        raise ValueError("--chains must contain exactly three chain IDs")
    if args.contig_ha_start > args.ha_end_residue:
        raise ValueError("--contig-ha-start must be <= --ha-end-residue")

    cif_paths = sorted(args.input_dir.glob("rank_*_*.cif"), key=parse_rank)
    if not cif_paths:
        raise ValueError(f"No rank_*_*.cif files found in {args.input_dir}")

    args.outdir.mkdir(parents=True, exist_ok=True)
    args.rf_inputs_dir.mkdir(parents=True, exist_ok=True)
    contig = contig_for_chains(
        chains,
        args.contig_ha_start,
        args.ha_end_residue,
        args.scaffold_start_residue,
        args.scaffold_max_residue,
        args.linker_spec,
    )

    rows: list[dict[str, str]] = []
    commands: list[str] = []
    for cif_path in cif_paths:
        atoms = parse_mmcif_atoms(cif_path)
        if not args.no_canonicalize:
            atoms = canonicalize_c3_axis(atoms, chains, args.ha_min_residue, args.ha_end_residue)

        pdb_path = args.outdir / f"{cif_path.stem}.pdb"
        rf_input_path = args.rf_inputs_dir / pdb_path.name
        remarks = [
            f"source_cif {cif_path}",
            f"canonicalized_c3_to_z {not args.no_canonicalize}",
            f"rf_contig {contig}",
        ]
        write_pdb(pdb_path, atoms, remarks)
        write_pdb(rf_input_path, atoms, remarks)

        output_prefix = f"outputs/{pdb_path.stem}/{pdb_path.stem}"
        if args.command_input_prefix:
            command_input_pdb = f"{args.command_input_prefix.rstrip('/')}/{pdb_path.name}"
        else:
            command_input_pdb = pdb_path.relative_to(REPO_DIR).as_posix()
        command = (
            "./scripts/run_inference.py --config-name=symmetry "
            'inference.symmetry="C3" '
            f"inference.input_pdb={command_input_pdb} "
            f"inference.output_prefix={output_prefix} "
            f"inference.num_designs={args.num_designs} "
            f"'contigmap.contigs=[{contig}]'"
        )
        commands.append(command)
        rows.append(
            {
                "source_cif": str(cif_path),
                "pdb": str(pdb_path),
                "rf_input_pdb": str(rf_input_path),
                "canonicalized_c3_to_z": str(not args.no_canonicalize),
                "linker_spec": args.linker_spec,
                "num_designs": str(args.num_designs),
                "rf_contig": contig,
                "example_output_prefix": output_prefix,
            }
        )

    index_path = args.outdir / "index.csv"
    with index_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    commands_path = args.outdir / "rfdiffusion_commands.sh"
    with commands_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("#!/bin/bash\n")
        handle.write("set -euo pipefail\n\n")
        for command in commands:
            handle.write(command + "\n")

    print(f"Wrote {len(rows)} PDB candidates to {args.outdir}")
    print(f"Wrote RFdiffusion input PDB copies to {args.rf_inputs_dir}")
    print(f"Index: {index_path}")
    print(f"RFdiffusion command drafts: {commands_path}")
    print(f"Contig: {contig}")


if __name__ == "__main__":
    main()
