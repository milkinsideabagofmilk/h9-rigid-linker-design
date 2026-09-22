#!/usr/bin/env python3
"""
Align a C3 HA trimer C terminus to a C3 scaffold N terminus and sample poses.

The script keeps the HA fixed, moves the scaffold, and writes one mmCIF per
sample plus a CSV summary. It only needs numpy and is intended for simple
PyMOL-exported mmCIF files with an _atom_site loop.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable

import numpy as np


ATOM_SITE_FIELDS = [
    "_atom_site.group_PDB",
    "_atom_site.id",
    "_atom_site.type_symbol",
    "_atom_site.label_atom_id",
    "_atom_site.label_alt_id",
    "_atom_site.label_comp_id",
    "_atom_site.label_asym_id",
    "_atom_site.label_entity_id",
    "_atom_site.label_seq_id",
    "_atom_site.pdbx_PDB_ins_code",
    "_atom_site.Cartn_x",
    "_atom_site.Cartn_y",
    "_atom_site.Cartn_z",
    "_atom_site.occupancy",
    "_atom_site.B_iso_or_equiv",
    "_atom_site.auth_seq_id",
    "_atom_site.auth_asym_id",
    "_atom_site.pdbx_PDB_model_num",
]


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
    label_seq_id: str
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

    @property
    def seq_num(self) -> int:
        return int(self.label_seq_id)


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
                    label_seq_id=parts[8],
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


def write_mmcif(path: Path, atoms: Iterable[Atom], data_name: str, comments: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for comment in comments:
            handle.write(f"# {comment}\n")
        handle.write(f"data_{data_name}\n")
        handle.write(f"_entry.id {data_name}\n")
        handle.write("#\nloop_\n")
        for field in ATOM_SITE_FIELDS:
            handle.write(f"{field}\n")
        for idx, atom in enumerate(atoms, start=1):
            handle.write(format_atom_line(replace(atom, atom_id=idx)))
        handle.write("#\n")


def format_atom_line(atom: Atom) -> str:
    auth_seq_id = atom.formal_charge
    if auth_seq_id in {"?", ".", "0"}:
        auth_seq_id = atom.label_seq_id
    return (
        f"{atom.group:<6} {atom.atom_id:>6d} {atom.type_symbol:<2} "
        f"{atom.atom_name:<4} {atom.alt_id:<1} {atom.res_name:<3} "
        f"{atom.label_asym_id:<2} {atom.label_entity_id:>2} "
        f"{atom.label_seq_id:>4} {atom.ins_code:<1} "
        f"{atom.x:>9.3f} {atom.y:>9.3f} {atom.z:>9.3f} "
        f"{atom.occupancy:>4} {atom.b_iso:>7} {auth_seq_id:>4} "
        f"{atom.chain_id:<2} {atom.model_num:>2}\n"
    )


def split_prealigned_atoms(
    atoms: list[Atom],
    chains: list[str],
    ha_end_residue: int,
    scaffold_start_residue: int,
) -> tuple[list[Atom], list[Atom]]:
    chain_set = set(chains)
    ha_atoms = [
        atom for atom in atoms if atom.chain_id in chain_set and atom.seq_num <= ha_end_residue
    ]
    scaffold_atoms = [
        atom for atom in atoms if atom.chain_id in chain_set and atom.seq_num >= scaffold_start_residue
    ]
    if not ha_atoms:
        raise ValueError(
            f"No HA atoms found at residues <= {ha_end_residue} in chains {','.join(chains)}"
        )
    if not scaffold_atoms:
        raise ValueError(
            "No scaffold atoms found at residues >= "
            f"{scaffold_start_residue} in chains {','.join(chains)}"
        )
    return ha_atoms, scaffold_atoms


def normalize(vector: np.ndarray, label: str = "vector") -> np.ndarray:
    norm = float(np.linalg.norm(vector))
    if norm < 1e-8:
        raise ValueError(f"Cannot normalize near-zero {label}")
    return vector / norm


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


def signed_angle_around_axis(source: np.ndarray, target: np.ndarray, axis: np.ndarray) -> float:
    axis = normalize(axis, "signed-angle axis")
    source = source - np.dot(source, axis) * axis
    target = target - np.dot(target, axis) * axis
    source = normalize(source, "source radial vector")
    target = normalize(target, "target radial vector")
    sin_value = float(np.dot(axis, np.cross(source, target)))
    cos_value = float(np.dot(source, target))
    return math.atan2(sin_value, cos_value)


def collect_ca_by_chain_residue(atoms: list[Atom], chains: list[str]) -> dict[str, dict[int, np.ndarray]]:
    ca: dict[str, dict[int, np.ndarray]] = {chain: {} for chain in chains}
    for atom in atoms:
        if atom.chain_id in ca and atom.atom_name == "CA" and atom.alt_id in (".", "A", "?"):
            ca[atom.chain_id][atom.seq_num] = atom.coord
    return ca


def fit_c3_axis(atoms: list[Atom], chains: list[str]) -> tuple[np.ndarray, np.ndarray]:
    ca = collect_ca_by_chain_residue(atoms, chains)
    missing = [chain for chain, residues in ca.items() if not residues]
    if missing:
        raise ValueError(f"Chains missing CA atoms: {', '.join(missing)}")

    common_residues = set(ca[chains[0]])
    for chain in chains[1:]:
        common_residues &= set(ca[chain])
    if len(common_residues) < 2:
        raise ValueError("Need at least two shared CA residue positions to fit a C3 axis")

    centers = np.array(
        [np.mean([ca[chain][seq] for chain in chains], axis=0) for seq in sorted(common_residues)]
    )
    axis_point = np.mean(centers, axis=0)
    _, _, vh = np.linalg.svd(centers - axis_point, full_matrices=False)
    axis = normalize(vh[0], "fitted C3 axis")
    return axis_point, axis


def terminal_points(
    atoms: list[Atom],
    chains: list[str],
    terminus: str,
    atom_name: str,
) -> dict[str, np.ndarray]:
    if terminus not in {"n", "c"}:
        raise ValueError("terminus must be 'n' or 'c'")

    points: dict[str, np.ndarray] = {}
    for chain in chains:
        chain_atoms = [atom for atom in atoms if atom.chain_id == chain]
        if not chain_atoms:
            raise ValueError(f"Chain {chain} not found")
        seqs = [atom.seq_num for atom in chain_atoms]
        terminal_seq = min(seqs) if terminus == "n" else max(seqs)
        candidates = [
            atom
            for atom in chain_atoms
            if atom.seq_num == terminal_seq and atom.atom_name == atom_name and atom.alt_id in (".", "A", "?")
        ]
        if not candidates:
            raise ValueError(
                f"Could not find atom {atom_name} at {terminus.upper()} terminus "
                f"residue {terminal_seq} in chain {chain}"
            )
        points[chain] = candidates[0].coord
    return points


def orient_axis_toward_terminal(
    axis: np.ndarray,
    atoms: list[Atom],
    terminal_centroid: np.ndarray,
) -> np.ndarray:
    ca_coords = np.array([atom.coord for atom in atoms if atom.atom_name == "CA"])
    if len(ca_coords) == 0:
        ca_coords = np.array([atom.coord for atom in atoms])
    body_centroid = np.mean(ca_coords, axis=0)
    if float(np.dot(axis, terminal_centroid - body_centroid)) < 0.0:
        return -axis
    return axis


def radial_basis(axis: np.ndarray, terminal_points_by_chain: dict[str, np.ndarray], chains: list[str]) -> tuple[np.ndarray, np.ndarray]:
    centroid = np.mean([terminal_points_by_chain[chain] for chain in chains], axis=0)
    e1 = terminal_points_by_chain[chains[0]] - centroid
    e1 = e1 - np.dot(e1, axis) * axis
    e1 = normalize(e1, f"radial vector for chain {chains[0]}")
    e2 = normalize(np.cross(axis, e1), "second radial basis vector")
    return e1, e2


def transform_atoms(
    atoms: list[Atom],
    origin: np.ndarray,
    base_rotation: np.ndarray,
    target_origin: np.ndarray,
    target_axis: np.ndarray,
    twist_deg: float,
) -> list[Atom]:
    twist = rotation_about_axis(target_axis, math.radians(twist_deg))
    transformed: list[Atom] = []
    for atom in atoms:
        aligned = origin + base_rotation @ (atom.coord - origin)
        shifted = aligned + (target_origin - origin)
        twisted = target_origin + twist @ (shifted - target_origin)
        transformed.append(replace(atom, x=float(twisted[0]), y=float(twisted[1]), z=float(twisted[2])))
    return transformed


def transform_points(
    points_by_chain: dict[str, np.ndarray],
    origin: np.ndarray,
    base_rotation: np.ndarray,
    target_origin: np.ndarray,
    target_axis: np.ndarray,
    twist_deg: float,
) -> dict[str, np.ndarray]:
    twist = rotation_about_axis(target_axis, math.radians(twist_deg))
    out = {}
    for chain, point in points_by_chain.items():
        aligned = origin + base_rotation @ (point - origin)
        shifted = aligned + (target_origin - origin)
        out[chain] = target_origin + twist @ (shifted - target_origin)
    return out


def parse_samples(spec: str) -> list[float]:
    values: list[float] = []
    for chunk in spec.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        if ":" not in chunk:
            values.append(float(chunk))
            continue
        parts = [float(part) for part in chunk.split(":")]
        if len(parts) != 3:
            raise ValueError(f"Range '{chunk}' must be start:end:step")
        start, end, step = parts
        if abs(step) < 1e-12:
            raise ValueError(f"Range '{chunk}' has zero step")
        current = start
        if step > 0:
            while current <= end + 1e-9:
                values.append(round(current, 10))
                current += step
        else:
            while current >= end - 1e-9:
                values.append(round(current, 10))
                current += step
    if not values:
        raise ValueError(f"No sample values parsed from '{spec}'")
    return values


def sample_count(*sample_lists: list[float]) -> int:
    count = 1
    for samples in sample_lists:
        count *= len(samples)
    return count


SCRIPT_DIR = Path(__file__).resolve().parent
RIGID_DIR = SCRIPT_DIR.parent
ALIGN_DIR = RIGID_DIR / "align"


def align_path(path: str) -> Path:
    return ALIGN_DIR / path


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Align HA C termini to scaffold N termini along their C3 axes, "
            "then sample axial translations and rotations around the common axis."
        )
    )
    parser.add_argument(
        "--prealigned",
        type=Path,
        default=align_path("CK_XT_24_scaffold_prealigned.cif"),
        help=(
            "Combined HA+scaffold mmCIF to split into HA and scaffold before sampling. "
            "Use this for CK_XT_24_scaffold_prealigned.cif."
        ),
    )
    parser.add_argument(
        "--ha-end-residue",
        type=int,
        default=525,
        help="Last HA residue when --prealigned is used",
    )
    parser.add_argument(
        "--scaffold-start-residue",
        type=int,
        default=526,
        help="First scaffold residue when --prealigned is used",
    )
    parser.add_argument("--ha", default=align_path("CK_XT_24_1_525.cif"), type=Path)
    parser.add_argument("--scaffold", default=align_path("ChainABC_519_644.cif"), type=Path)
    parser.add_argument("--outdir", default=align_path("ck_xt_24_scaffold_samples"), type=Path)
    parser.add_argument("--chains", default="A,B,C", help="Comma-separated chain IDs in C3 order")
    parser.add_argument("--ha-terminal-atom", default="C", help="Atom at the HA C terminus")
    parser.add_argument("--scaffold-terminal-atom", default="N", help="Atom at the scaffold N terminus")
    parser.add_argument(
        "--z-shifts",
        default="20,24,28,32",
        help="Nonnegative axial translations in Angstroms, as values or start:end:step",
    )
    parser.add_argument(
        "--rotations",
        default="0:350:10",
        help=(
            "Twist rotations in degrees around the common C3 axis. "
            "The default covers the full 360-degree phase because fixed chain "
            "identity makes 120-degree rotations non-equivalent."
        ),
    )
    parser.add_argument("--x-shifts", default="0", help="Optional radial translations along chain-A direction")
    parser.add_argument("--y-shifts", default="0", help="Optional radial translations orthogonal to x-shifts")
    parser.add_argument(
        "--orientation",
        choices=["face-to-face", "same-direction"],
        default="face-to-face",
        help="Use face-to-face for C-to-N fusion geometry",
    )
    parser.add_argument(
        "--no-phase-align",
        action="store_true",
        help="Do not pre-align scaffold chain A radial vector to HA chain A at zero twist",
    )
    parser.add_argument(
        "--max-samples",
        type=int,
        default=1000,
        help="Abort if the Cartesian product produces more samples than this",
    )
    return parser


def main() -> None:
    args = make_parser().parse_args()
    chains = [chain.strip() for chain in args.chains.split(",") if chain.strip()]
    if len(chains) != 3:
        raise ValueError("--chains must contain exactly three chain IDs")

    z_shifts = parse_samples(args.z_shifts)
    rotations = parse_samples(args.rotations)
    x_shifts = parse_samples(args.x_shifts)
    y_shifts = parse_samples(args.y_shifts)
    if any(z_shift < -1e-9 for z_shift in z_shifts):
        raise ValueError("--z-shifts must be >= 0 so the termini are only moved farther apart")
    n_samples = sample_count(z_shifts, rotations, x_shifts, y_shifts)
    if n_samples > args.max_samples:
        raise ValueError(
            f"Requested {n_samples} samples, above --max-samples {args.max_samples}. "
            "Raise --max-samples if this is intentional."
        )

    if args.prealigned:
        prealigned_atoms = parse_mmcif_atoms(args.prealigned)
        ha_atoms, scaffold_atoms = split_prealigned_atoms(
            prealigned_atoms,
            chains,
            args.ha_end_residue,
            args.scaffold_start_residue,
        )
    else:
        ha_atoms = parse_mmcif_atoms(args.ha)
        scaffold_atoms = parse_mmcif_atoms(args.scaffold)

    ha_axis_point, ha_axis = fit_c3_axis(ha_atoms, chains)
    scaffold_axis_point, scaffold_axis = fit_c3_axis(scaffold_atoms, chains)

    ha_c_points = terminal_points(ha_atoms, chains, "c", args.ha_terminal_atom)
    scaffold_n_points = terminal_points(scaffold_atoms, chains, "n", args.scaffold_terminal_atom)
    ha_c_centroid = np.mean([ha_c_points[chain] for chain in chains], axis=0)
    scaffold_n_centroid = np.mean([scaffold_n_points[chain] for chain in chains], axis=0)

    ha_axis = orient_axis_toward_terminal(ha_axis, ha_atoms, ha_c_centroid)
    scaffold_axis = orient_axis_toward_terminal(scaffold_axis, scaffold_atoms, scaffold_n_centroid)
    desired_scaffold_axis = -ha_axis if args.orientation == "face-to-face" else ha_axis

    axis_rotation = rotation_between_vectors(scaffold_axis, desired_scaffold_axis)
    base_rotation = axis_rotation
    if not args.no_phase_align:
        ha_e1, _ = radial_basis(ha_axis, ha_c_points, chains)
        scaffold_e1, _ = radial_basis(scaffold_axis, scaffold_n_points, chains)
        scaffold_e1_aligned = axis_rotation @ scaffold_e1
        phase_angle = signed_angle_around_axis(scaffold_e1_aligned, ha_e1, ha_axis)
        base_rotation = rotation_about_axis(ha_axis, phase_angle) @ axis_rotation

    x_axis, y_axis = radial_basis(ha_axis, ha_c_points, chains)
    args.outdir.mkdir(parents=True, exist_ok=True)
    summary_path = args.outdir / "summary.csv"

    rows: list[dict[str, object]] = []
    sample_id = 0
    for z_shift in z_shifts:
        for x_shift in x_shifts:
            for y_shift in y_shifts:
                for rotation in rotations:
                    sample_id += 1
                    target_origin = (
                        ha_c_centroid
                        + z_shift * ha_axis
                        + x_shift * x_axis
                        + y_shift * y_axis
                    )
                    moved_scaffold = transform_atoms(
                        scaffold_atoms,
                        scaffold_n_centroid,
                        base_rotation,
                        target_origin,
                        ha_axis,
                        rotation,
                    )
                    moved_n_points = transform_points(
                        scaffold_n_points,
                        scaffold_n_centroid,
                        base_rotation,
                        target_origin,
                        ha_axis,
                        rotation,
                    )
                    distances = [
                        float(np.linalg.norm(ha_c_points[chain] - moved_n_points[chain]))
                        for chain in chains
                    ]
                    data_name = f"sample_{sample_id:04d}"
                    out_path = args.outdir / f"{data_name}.cif"
                    comments = [
                        "Generated by align/sample_c3_axis.py",
                        (
                            f"Prealigned source split at HA <= {args.ha_end_residue} "
                            f"and scaffold >= {args.scaffold_start_residue}: {args.prealigned}"
                            if args.prealigned
                            else f"HA fixed: {args.ha}"
                        ),
                        (
                            f"Moved scaffold from prealigned source: {args.prealigned}"
                            if args.prealigned
                            else f"Moved scaffold: {args.scaffold}"
                        ),
                        f"z_shift_A={z_shift:.3f} x_shift_A={x_shift:.3f} "
                        f"y_shift_A={y_shift:.3f} rotation_deg={rotation:.3f}",
                        f"orientation={args.orientation}",
                    ]
                    write_mmcif(out_path, [*ha_atoms, *moved_scaffold], data_name, comments)
                    rows.append(
                        {
                            "sample_id": data_name,
                            "file": str(out_path),
                            "z_shift_A": z_shift,
                            "x_shift_A": x_shift,
                            "y_shift_A": y_shift,
                            "rotation_deg": rotation,
                            "endpoint_rmsd_A": math.sqrt(sum(d * d for d in distances) / len(distances)),
                            "endpoint_mean_A": sum(distances) / len(distances),
                            "endpoint_min_A": min(distances),
                            "endpoint_max_A": max(distances),
                            **{f"endpoint_{chain}_A": dist for chain, dist in zip(chains, distances)},
                        }
                    )

    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {sample_id} sampled poses to {args.outdir}")
    print(f"Summary: {summary_path}")
    print(f"HA C-terminal centroid: {ha_c_centroid[0]:.3f} {ha_c_centroid[1]:.3f} {ha_c_centroid[2]:.3f}")
    print(
        "HA axis toward C termini: "
        f"{ha_axis[0]:.5f} {ha_axis[1]:.5f} {ha_axis[2]:.5f}"
    )
    print(
        "Scaffold axis point: "
        f"{scaffold_axis_point[0]:.3f} {scaffold_axis_point[1]:.3f} {scaffold_axis_point[2]:.3f}"
    )
    print(
        "HA axis point: "
        f"{ha_axis_point[0]:.3f} {ha_axis_point[1]:.3f} {ha_axis_point[2]:.3f}"
    )


if __name__ == "__main__":
    main()
