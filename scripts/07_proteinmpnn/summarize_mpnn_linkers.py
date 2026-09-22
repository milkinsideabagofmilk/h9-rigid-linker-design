#!/usr/bin/env python3
"""Aggregate MPNN C3-tied linker sequences into mpnn_linker_sequences.csv.

Reads output/mpnn_out_c3_tied/seqs/*.fa (native + 2 designed entries per
backbone), extracts the designed linker (chain-internal residues 47-66) from
each of the A/B/C chains, and verifies C3 symmetry (three linkers identical).
CSV format matches the pre-relax version:
target,sample,score,global_score,linker_20aa,chainA_linker,chainB_linker,chainC_linker
"""
import argparse
import csv
import re
from pathlib import Path

STEP = Path(__file__).resolve().parents[1]


def parse_fa(path):
    entries, header, chunks = [], None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                entries.append((header, "".join(chunks)))
            header, chunks = line[1:], []
        else:
            chunks.append(line.strip())
    if header is not None:
        entries.append((header, "".join(chunks)))
    return entries


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seqs-dir", type=Path, default=STEP / "output" / "mpnn_out_c3_tied" / "seqs")
    ap.add_argument("--out", type=Path, default=STEP / "output" / "mpnn_linker_sequences.csv")
    args = ap.parse_args()

    rows, asymmetric = [], []
    for fa in sorted(args.seqs_dir.glob("*.fa")):
        for header, seq in parse_fa(fa):
            m = re.search(r"sample=(\d+)", header)
            if not m:
                continue  # skip native sequence entry
            chains = seq.split("/")
            assert len(chains) == 3, f"{fa.name}: expected 3 chains, got {len(chains)}"
            linkers = [c[46:66] for c in chains]  # chain-internal 47-66
            if len(set(linkers)) != 1:
                asymmetric.append((fa.stem, m.group(1)))
            rows.append({
                "target": fa.stem,
                "sample": f"sample={m.group(1)}",
                "score": re.search(r"score=([\d.]+)", header).group(1),
                "global_score": re.search(r"global_score=([\d.]+)", header).group(1),
                "linker_20aa": linkers[0],
                "chainA_linker": linkers[0],
                "chainB_linker": linkers[1],
                "chainC_linker": linkers[2],
            })

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["target", "sample", "score", "global_score",
                                           "linker_20aa", "chainA_linker",
                                           "chainB_linker", "chainC_linker"])
        w.writeheader()
        w.writerows(rows)

    scores = [float(r["score"]) for r in rows]
    print(f"wrote {args.out}: {len(rows)} designed sequences from {len(set(r['target'] for r in rows))} backbones")
    print(f"score range {min(scores):.4f}..{max(scores):.4f}")
    if asymmetric:
        print(f"WARNING: {len(asymmetric)} entries violate C3-tied symmetry: {asymmetric}")
    else:
        print("C3 symmetry: all three chain linkers identical in every entry")


if __name__ == "__main__":
    main()
