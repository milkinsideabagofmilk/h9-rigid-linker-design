#!/usr/bin/env python3
"""Merge AF2 metrics, structural checks, and relaxed-Rosetta energies into the
final candidate table (08_af2_validation/output/final_candidates.csv).

Final candidate = ipTM>=0.65 AND struct_pass. Sorted by ipTM desc.
"""
import argparse
import csv
from pathlib import Path

STEP = Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--af2", type=Path,
                    default=STEP / "output" / "results_conservative_summary.csv")
    ap.add_argument("--struct", type=Path, default=STEP / "output" / "structural_checks.csv")
    ap.add_argument("--rosetta", type=Path, default=STEP / "output" / "rosetta_rescore"
                    / "ck_xt_24_rosetta_linker_scores_relaxed.csv")
    ap.add_argument("--out", type=Path, default=STEP / "output" / "final_candidates.csv")
    args = ap.parse_args()

    struct = {r["design_id"]: r for r in csv.DictReader(args.struct.open())}
    rosetta = {}
    if args.rosetta.exists():
        for r in csv.DictReader(args.rosetta.open()):
            if r["status"] == "PASS":
                design = r["name"].split("_unrelaxed_rank_001_")[0]
                rosetta[design] = r

    rows = []
    for r in csv.DictReader(args.af2.open()):
        d = r["design_id"]
        s = struct.get(d, {})
        ro = rosetta.get(d, {})
        rows.append({
            "design_id": d,
            "backbone": d.rsplit("_sample", 1)[0],
            "sample": r["sample"],
            "iptm": r["iptm"], "ptm": r["ptm"],
            "plddt_mean": r["plddt_mean"], "linker_plddt": r["linker_plddt_mean"],
            "struct_pass": s.get("struct_pass", ""),
            "bend": s.get("bend_max", ""), "sym_rmsd": s.get("sym_rmsd_max", ""),
            "core_rmsd": s.get("full_rmsd_max", ""),
            "linker_rmsd_info": s.get("linker_rmsd_max", ""),
            "rosetta_flankE_relaxed": ro.get("post_linker_flank_total_energy", ""),
            "rosetta_flankE_prerelax": ro.get("pre_linker_flank_total_energy", ""),
            "final_candidate": (float(r["iptm"]) >= 0.65 and s.get("struct_pass") == "True"),
        })
    rows.sort(key=lambda r: -float(r["iptm"]))

    cols = list(rows[0])
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)

    finals = [r for r in rows if r["final_candidate"]]
    print(f"{len(rows)} designs -> {len(finals)} final candidates  ({args.out})")
    print(f"{'design_id':<40} {'ipTM':>5} {'coreR':>5} {'symR':>5} {'flankE':>7}")
    for r in finals:
        print(f"{r['design_id']:<40} {r['iptm']:>5} {r['core_rmsd']:>5} "
              f"{r['sym_rmsd']:>5} {r['rosetta_flankE_relaxed']:>7}")


if __name__ == "__main__":
    main()
