#!/usr/bin/env python3
"""Summarize ColabFold/AF2-multimer results into results_conservative_summary.csv.

Columns match the pre-relax version:
design_id,source_target,sample,plddt_mean,linker_plddt_mean,
linker_flank_plddt_mean,ptm,iptm,pae_mean,inter_chain_pae_mean,
pdb_path,score_json_path

Residue layout (0-based, flat over A/B/C chains of 193 aa each):
linker       = 46-65, 239-258, 432-451   (chain-internal 47-66)
linker_flank = 43-68, 236-261, 429-454   (linker +/- 3)
Sorted by iptm descending.
"""
import argparse
import csv
import json
from pathlib import Path

STEP = Path(__file__).resolve().parents[1]

CHAIN_LEN = 193
LINKER = [i for c in range(3) for i in range(c * CHAIN_LEN + 46, c * CHAIN_LEN + 66)]
FLANK = [i for c in range(3) for i in range(c * CHAIN_LEN + 43, c * CHAIN_LEN + 69)]


def mean(xs):
    return sum(xs) / len(xs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", type=Path, default=STEP / "output" / "results_conservative")
    ap.add_argument("--out", type=Path,
                    default=STEP / "output" / "results_conservative_summary.csv")
    args = ap.parse_args()

    rows = []
    for done in sorted(args.results_dir.glob("*.done.txt")):
        design = done.name[: -len(".done.txt")]
        score_jsons = sorted(args.results_dir.glob(f"{design}_scores_rank_001_*.json"))
        pae_jsons = sorted(args.results_dir.glob(f"{design}_predicted_aligned_error_v1.json"))
        pdbs = sorted(args.results_dir.glob(f"{design}_unrelaxed_rank_001_*.pdb"))
        if not (score_jsons and pae_jsons and pdbs):
            print(f"SKIP {design}: missing outputs")
            continue
        s = json.loads(score_jsons[0].read_text())
        plddt = s["plddt"]
        pae = json.loads(pae_jsons[0].read_text())
        mat = pae["predicted_aligned_error"]
        n = len(mat)
        inter = [mat[i][j] for i in range(n) for j in range(n)
                 if i // CHAIN_LEN != j // CHAIN_LEN]
        target, sample = design.rsplit("_sample", 1)
        rows.append({
            "design_id": design,
            "source_target": target,
            "sample": f"sample={sample}",
            "plddt_mean": f"{mean(plddt):.2f}",
            "linker_plddt_mean": f"{mean([plddt[i] for i in LINKER]):.2f}",
            "linker_flank_plddt_mean": f"{mean([plddt[i] for i in FLANK]):.2f}",
            "ptm": f"{s['ptm']:.4f}",
            "iptm": f"{s['iptm']:.4f}",
            "pae_mean": f"{mean([v for row in mat for v in row]):.2f}",
            "inter_chain_pae_mean": f"{mean(inter):.2f}",
            "pdb_path": pdbs[0].relative_to(STEP.parent),
            "score_json_path": score_jsons[0].relative_to(STEP.parent),
        })

    rows.sort(key=lambda r: -float(r["iptm"]))
    cols = ["design_id", "source_target", "sample", "plddt_mean", "linker_plddt_mean",
            "linker_flank_plddt_mean", "ptm", "iptm", "pae_mean",
            "inter_chain_pae_mean", "pdb_path", "score_json_path"]
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)

    n_pass = sum(1 for r in rows if float(r["iptm"]) >= 0.65)
    print(f"wrote {args.out}: {len(rows)} designs, {n_pass} with ipTM>=0.65")
    for r in rows:
        print(f"  {r['design_id']:<40} ipTM={r['iptm']} pTM={r['ptm']} "
              f"pLDDT={r['plddt_mean']} linker_pLDDT={r['linker_plddt_mean']} "
              f"interPAE={r['inter_chain_pae_mean']}")


if __name__ == "__main__":
    main()
