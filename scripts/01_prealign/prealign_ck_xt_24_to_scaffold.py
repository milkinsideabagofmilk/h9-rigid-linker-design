#!/usr/bin/env python3
"""
Pre-align CK_XT_24 HA(1-525) to ChainABC_518_644 scaffold on a shared C3 axis.

The HA coordinates are kept fixed. The scaffold is rigid-body transformed so
its C3 axis is collinear with the HA C3 axis in face-to-face C-to-N fusion
geometry, then translated so scaffold N termini sit within the requested
distance from HA residue 525 C termini.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_DIR = SCRIPT_DIR.parent
ALIGN_DIR = REPO_DIR / "align"
sys.path.insert(0, str(SCRIPT_DIR))

from sample_c3_axis import (  # noqa: E402
    Atom,
    fit_c3_axis,
    orient_axis_toward_terminal,
    parse_mmcif_atoms,
    rotation_about_axis,
    rotation_between_vectors,
    signed_angle_around_axis,
    terminal_points,
)


AF_ATOM_SITE_FIELDS = [
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


def truncate_atoms(atoms: list[Atom], max_residue: int) -> list[Atom]:
    return [atom for atom in atoms if atom.seq_num <= max_residue]


def renumber_atoms(
    atoms: list[Atom],
    new_start: int,
    entity_by_chain: dict[str, str],
) -> list[Atom]:
    min_residue = min(atom.seq_num for atom in atoms)
    offset = new_start - min_residue
    return [
        replace(
            atom,
            label_seq_id=str(atom.seq_num + offset),
            label_entity_id=entity_by_chain.get(atom.chain_id, atom.label_entity_id),
            formal_charge=str(atom.seq_num + offset),
        )
        for atom in atoms
    ]


def find_atom_site_fields(lines: list[str]) -> tuple[int, list[str]]:
    for idx, line in enumerate(lines):
        if line.startswith("_atom_site."):
            fields: list[str] = []
            cursor = idx
            while cursor < len(lines) and lines[cursor].startswith("_atom_site."):
                fields.append(lines[cursor].strip())
                cursor += 1
            return cursor, fields
    raise ValueError("Could not find _atom_site fields")


def filter_entity_poly_seq_prefix(lines: list[str], max_residue: int) -> list[str]:
    out: list[str] = []
    i = 0
    while i < len(lines):
        if (
            lines[i].strip() == "loop_"
            and i + 1 < len(lines)
            and lines[i + 1].startswith("_entity_poly_seq.")
        ):
            out.append(lines[i])
            i += 1
            fields: list[str] = []
            while i < len(lines) and lines[i].startswith("_entity_poly_seq."):
                fields.append(lines[i].strip())
                out.append(lines[i])
                i += 1
            try:
                num_idx = fields.index("_entity_poly_seq.num")
            except ValueError as exc:
                raise ValueError("Could not find _entity_poly_seq.num") from exc
            while i < len(lines) and lines[i].strip() != "#":
                parts = lines[i].split()
                if parts and int(parts[num_idx]) <= max_residue:
                    out.append(lines[i])
                i += 1
            if i < len(lines):
                out.append(lines[i])
                i += 1
            continue
        out.append(lines[i])
        i += 1
    return out


def write_truncated_original_cif(source: Path, out_path: Path, max_residue: int) -> None:
    lines = source.read_text(encoding="utf-8").splitlines()
    atom_start, fields = find_atom_site_fields(lines)
    field_index = {field: idx for idx, field in enumerate(fields)}
    seq_idx = field_index["_atom_site.label_seq_id"]
    atom_id_idx = field_index["_atom_site.id"]

    prefix = filter_entity_poly_seq_prefix(lines[:atom_start], max_residue)
    filtered_atom_lines: list[str] = []
    atom_id = 0
    for line in lines[atom_start:]:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        parts = line.split()
        if int(parts[seq_idx]) > max_residue:
            continue
        atom_id += 1
        parts[atom_id_idx] = str(atom_id)
        filtered_atom_lines.append(" ".join(parts))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f"# Truncated from {source}\n")
        handle.write(f"# Kept residues 1-{max_residue} and preserved original atom_site schema\n")
        for line in prefix:
            handle.write(f"{line}\n")
        for line in filtered_atom_lines:
            handle.write(f"{line}\n")
        handle.write("#\n")


def format_af_atom_line(atom: Atom, atom_id: int) -> str:
    auth_seq_id = atom.label_seq_id
    return (
        f"{atom.group:<6} {atom_id:<5d} {atom.type_symbol:<2} "
        f"{atom.atom_name:<4} {atom.alt_id:<1} {atom.res_name:<3} "
        f"{atom.label_asym_id:<2} {atom.label_entity_id:>2} "
        f"{atom.label_seq_id:>4} {atom.ins_code:<1} "
        f"{atom.x:>9.3f} {atom.y:>9.3f} {atom.z:>9.3f} "
        f"{atom.occupancy:>4} {atom.b_iso:>7} {auth_seq_id:>4} "
        f"{atom.chain_id:<2} {atom.model_num:>2}\n"
    )


def write_af_atom_site_mmcif(
    path: Path,
    atoms: list[Atom],
    data_name: str,
    comments: list[str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for comment in comments:
            handle.write(f"# {comment}\n")
        handle.write(f"data_{data_name}\n")
        handle.write(f"_entry.id {data_name}\n")
        handle.write("#\nloop_\n")
        handle.write("_struct_asym.entity_id\n")
        handle.write("_struct_asym.id\n")
        for entity_id, chain in (("1", "A"), ("2", "B"), ("3", "C")):
            handle.write(f"{entity_id} {chain}\n")
        handle.write("#\nloop_\n")
        for field in AF_ATOM_SITE_FIELDS:
            handle.write(f"{field}\n")
        for atom_id, atom in enumerate(atoms, start=1):
            handle.write(format_af_atom_line(atom, atom_id))
        handle.write("#\n")



def chain_atom(atoms: list[Atom], chain: str, seq: int, atom_name: str) -> np.ndarray:
    for atom in atoms:
        if (
            atom.chain_id == chain
            and atom.seq_num == seq
            and atom.atom_name == atom_name
            and atom.alt_id in (".", "A", "?")
        ):
            return atom.coord
    raise ValueError(f"Missing atom {chain}:{seq}:{atom_name}")


def point_projection_on_axis(point: np.ndarray, axis_point: np.ndarray, axis: np.ndarray) -> np.ndarray:
    return axis_point + np.dot(point - axis_point, axis) * axis


def radial_from_axis(point: np.ndarray, axis_point: np.ndarray, axis: np.ndarray) -> np.ndarray:
    radial = point - point_projection_on_axis(point, axis_point, axis)
    norm = float(np.linalg.norm(radial))
    if norm < 1e-8:
        raise ValueError("Cannot define radial vector for point on axis")
    return radial / norm


def transform_atoms_on_axis(
    atoms: list[Atom],
    source_axis_point: np.ndarray,
    total_rotation: np.ndarray,
    target_axis_point: np.ndarray,
) -> list[Atom]:
    transformed: list[Atom] = []
    for atom in atoms:
        coord = target_axis_point + total_rotation @ (atom.coord - source_axis_point)
        transformed.append(replace(atom, x=float(coord[0]), y=float(coord[1]), z=float(coord[2])))
    return transformed


def transform_points_on_axis(
    points_by_chain: dict[str, np.ndarray],
    source_axis_point: np.ndarray,
    total_rotation: np.ndarray,
    target_axis_point: np.ndarray,
) -> dict[str, np.ndarray]:
    return {
        chain: target_axis_point + total_rotation @ (point - source_axis_point)
        for chain, point in points_by_chain.items()
    }


def best_twist(
    scaffold_n_points: dict[str, np.ndarray],
    scaffold_axis_point: np.ndarray,
    base_rotation: np.ndarray,
    target_axis_point: np.ndarray,
    target_axis: np.ndarray,
    ha_c_points: dict[str, np.ndarray],
    chains: list[str],
    step_deg: float,
) -> tuple[float, dict[str, np.ndarray], list[float]]:
    best_rotation = 0.0
    best_points: dict[str, np.ndarray] | None = None
    best_distances: list[float] | None = None
    best_key: tuple[float, float] | None = None

    n_steps = max(1, int(round(360.0 / step_deg)))
    for i in range(n_steps):
        rotation = i * 360.0 / n_steps
        twist = rotation_about_axis(target_axis, math.radians(rotation))
        moved_n_points = transform_points_on_axis(
            scaffold_n_points,
            scaffold_axis_point,
            twist @ base_rotation,
            target_axis_point,
        )
        distances = [
            float(np.linalg.norm(ha_c_points[chain] - moved_n_points[chain]))
            for chain in chains
        ]
        key = (max(distances), math.sqrt(sum(d * d for d in distances) / len(distances)))
        if best_key is None or key < best_key:
            best_key = key
            best_rotation = rotation
            best_points = moved_n_points
            best_distances = distances

    if best_points is None or best_distances is None:
        raise RuntimeError("Failed to search twist rotations")
    return best_rotation, best_points, best_distances


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ha", default=REPO_DIR / "CK_XT_24.cif", type=Path)
    parser.add_argument("--scaffold", default=REPO_DIR / "ChainABC_518_644.cif", type=Path)
    parser.add_argument("--chains", default="A,B,C")
    parser.add_argument("--ha-end-residue", default=525, type=int)
    parser.add_argument("--target-distance", default=19.0, type=float)
    parser.add_argument("--twist-step-deg", default=1.0, type=float)
    parser.add_argument("--renumber-scaffold-start", default=526, type=int)
    parser.add_argument(
        "--truncated-out",
        default=ALIGN_DIR / "CK_XT_24_1_525.cif",
        type=Path,
    )
    parser.add_argument(
        "--combined-out",
        default=ALIGN_DIR / "CK_XT_24_scaffold_prealigned.cif",
        type=Path,
    )
    parser.add_argument(
        "--metrics-out",
        default=ALIGN_DIR / "CK_XT_24_scaffold_prealigned_metrics.csv",
        type=Path,
    )
    return parser


def main() -> None:
    args = make_parser().parse_args()
    chains = [chain.strip() for chain in args.chains.split(",") if chain.strip()]
    if len(chains) != 3:
        raise ValueError("--chains must contain exactly three chain IDs")
    if args.target_distance >= 20.0:
        raise ValueError("--target-distance should be below 20 A to keep all termini within 20 A")

    ha_atoms = truncate_atoms(parse_mmcif_atoms(args.ha), args.ha_end_residue)
    scaffold_atoms = parse_mmcif_atoms(args.scaffold)

    ha_axis_point, ha_axis = fit_c3_axis(ha_atoms, chains)
    scaffold_axis_point, scaffold_axis = fit_c3_axis(scaffold_atoms, chains)

    ha_c_points = {
        chain: chain_atom(ha_atoms, chain, args.ha_end_residue, "C") for chain in chains
    }
    scaffold_n_points = terminal_points(scaffold_atoms, chains, "n", "N")
    ha_c_centroid = np.mean([ha_c_points[chain] for chain in chains], axis=0)
    scaffold_n_centroid = np.mean([scaffold_n_points[chain] for chain in chains], axis=0)

    ha_axis = orient_axis_toward_terminal(ha_axis, ha_atoms, ha_c_centroid)
    scaffold_axis = orient_axis_toward_terminal(scaffold_axis, scaffold_atoms, scaffold_n_centroid)
    desired_scaffold_axis = -ha_axis

    axis_rotation = rotation_between_vectors(scaffold_axis, desired_scaffold_axis)
    ha_e1 = radial_from_axis(ha_c_points[chains[0]], ha_axis_point, ha_axis)
    scaffold_e1 = radial_from_axis(
        scaffold_n_points[chains[0]], scaffold_axis_point, scaffold_axis
    )
    phase_angle = signed_angle_around_axis(axis_rotation @ scaffold_e1, ha_e1, ha_axis)
    base_rotation = rotation_about_axis(ha_axis, phase_angle) @ axis_rotation

    ha_c_projection = point_projection_on_axis(ha_c_centroid, ha_axis_point, ha_axis)
    rotated_scaffold_n_centroid_offset = base_rotation @ (scaffold_n_centroid - scaffold_axis_point)
    desired_scaffold_n_projection = ha_c_projection + args.target_distance * ha_axis
    target_axis_point = (
        desired_scaffold_n_projection
        - np.dot(rotated_scaffold_n_centroid_offset, ha_axis) * ha_axis
    )
    best_rotation, moved_n_points, distances = best_twist(
        scaffold_n_points,
        scaffold_axis_point,
        base_rotation,
        target_axis_point,
        ha_axis,
        ha_c_points,
        chains,
        args.twist_step_deg,
    )
    if max(distances) >= 20.0:
        raise ValueError(
            f"Best HA C to scaffold N distance is {max(distances):.3f} A; "
            "lower --target-distance and rerun"
        )

    moved_scaffold = transform_atoms_on_axis(
        scaffold_atoms,
        scaffold_axis_point,
        rotation_about_axis(ha_axis, math.radians(best_rotation)) @ base_rotation,
        target_axis_point,
    )
    entity_by_chain = {chain: str(idx + 1) for idx, chain in enumerate(chains)}
    moved_scaffold = renumber_atoms(moved_scaffold, args.renumber_scaffold_start, entity_by_chain)

    moved_axis_point = target_axis_point
    moved_axis = rotation_about_axis(ha_axis, math.radians(best_rotation)) @ (base_rotation @ scaffold_axis)
    axis_alignment_angle_deg = math.degrees(
        math.acos(min(1.0, max(-1.0, abs(float(np.dot(moved_axis, ha_axis))))))
    )
    axis_offset = float(
        np.linalg.norm(np.cross(moved_axis, moved_axis_point - ha_axis_point))
    )

    write_truncated_original_cif(
        args.ha,
        args.truncated_out,
        args.ha_end_residue,
    )
    write_af_atom_site_mmcif(
        args.combined_out,
        [*ha_atoms, *moved_scaffold],
        "CK_XT_24_1_525_ChainABC_526_652_prealigned",
        [
            "Generated by align/prealign_ck_xt_24_to_scaffold.py",
            f"HA fixed: {args.ha}; kept residues 1-{args.ha_end_residue}",
            f"Moved scaffold: {args.scaffold}",
            f"Scaffold source residues renumbered to start at {args.renumber_scaffold_start}",
            "Face-to-face C3 axis alignment for HA C terminus to scaffold N terminus",
            f"target_centroid_distance_A={args.target_distance:.3f}",
            f"twist_deg={best_rotation:.3f}",
            f"axis_alignment_angle_deg={axis_alignment_angle_deg:.8f}",
            f"axis_line_offset_A={axis_offset:.8f}",
        ],
    )

    args.metrics_out.parent.mkdir(parents=True, exist_ok=True)
    rows = [
        {"metric": "axis_alignment_angle_deg", "value": axis_alignment_angle_deg},
        {"metric": "axis_line_offset_A", "value": axis_offset},
        {"metric": "target_centroid_distance_A", "value": args.target_distance},
        {"metric": "twist_deg", "value": best_rotation},
        {"metric": "endpoint_mean_A", "value": sum(distances) / len(distances)},
        {
            "metric": "endpoint_rmsd_A",
            "value": math.sqrt(sum(d * d for d in distances) / len(distances)),
        },
        {"metric": "endpoint_min_A", "value": min(distances)},
        {"metric": "endpoint_max_A", "value": max(distances)},
    ]
    rows.extend(
        {"metric": f"endpoint_{chain}_A", "value": distance}
        for chain, distance in zip(chains, distances)
    )
    with args.metrics_out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "value"])
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote truncated HA: {args.truncated_out}")
    print(f"Wrote prealigned complex: {args.combined_out}")
    print(f"Wrote metrics: {args.metrics_out}")
    print(f"Axis alignment angle: {axis_alignment_angle_deg:.8f} deg")
    print(f"Axis line offset: {axis_offset:.8f} A")
    print(
        "HA C to scaffold N distances: "
        + ", ".join(f"{chain}={distance:.3f} A" for chain, distance in zip(chains, distances))
    )


if __name__ == "__main__":
    main()
