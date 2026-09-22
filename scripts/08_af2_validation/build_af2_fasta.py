#!/usr/bin/env python3
"""Build the AF2-multimer input FASTA from MPNN C3-tied designs.

Each entry = 3 identical chains (MPNN output already contains the C-terminal
His6 tag inside the 193 aa chain) joined by ':', named <target>_sample<N>.
Format identical to the pre-relax version (verified against the archived
deprecated_unrelaxed_top10/input/ck_xt_24_mpnn_c3_tied_20.fasta).
"""
import argparse
from pathlib import Path

STEP = Path(__file__).resolve().parents[1]
MPNN_SEQS = STEP.parent / "07_proteinmpnn" / "output" / "mpnn_out_c3_tied" / "seqs"


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
    ap.add_argument("--seqs-dir", type=Path, default=MPNN_SEQS)
    ap.add_argument("--out", type=Path,
                    default=STEP / "input" / "ck_xt_24_mpnn_c3_tied_20.fasta")
    args = ap.parse_args()

    n = 0
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        for fa in sorted(args.seqs_dir.glob("*.fa")):
            for header, seq in parse_fa(fa):
                if "sample=" not in header:
                    continue  # native entry
                sample = header.split("sample=")[1].split(",")[0].split()[0]
                chains = seq.split("/")
                assert len(chains) == 3, f"{fa.name}: {len(chains)} chains"
                assert chains[0] == chains[1] == chains[2], (
                    f"{fa.name} sample{sample}: designed chains differ")
                assert len(chains[0]) == 193, f"{fa.name}: chain len {len(chains[0])}"
                assert chains[0].endswith("HHHHHH"), f"{fa.name}: no His6 tag"
                fh.write(f">{fa.stem}_sample{sample}\n")
                fh.write(":".join([chains[0]] * 3) + "\n")
                n += 1
    print(f"wrote {args.out} with {n} trimer entries (3x193 aa + His6, ':'-joined)")


if __name__ == "__main__":
    main()
