#!/usr/bin/env python3
"""FastRelax + linker/flank energy rescoring of RFdiffusion backbones.

For each watch-passed PDB: score linker regions, FastRelax with CA coordinate
constraints (sd 0.5, repeats 1), then score again in a single rosetta_scripts
run. Aggregates a CSV and compares the relaxed ranking against the old
unrelaxed ranking (06_rosetta_scoring/output/ck_xt_24_rosetta_linker_scores.csv).
"""
import argparse
import csv
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

# Rosetta is not redistributed; point this at a local release-427 (or later)
# build, e.g. via the ROSETTA_BIN environment variable or --rosetta.
ROSETTA = Path(os.environ.get(
    "ROSETTA_BIN",
    "/path/to/rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease",
))
STEP_DIR = Path(__file__).resolve().parents[1]  # 06_rosetta_scoring/
REPO = STEP_DIR.parent
PROTOCOL = STEP_DIR / "protocol" / "linker_score_relaxed.xml"
PASS_LIST = REPO / "05_watch_prefilter" / "output" / "ck_xt_24_pass_all_rigid_linkers.txt"
OLD_CSV = STEP_DIR / "output" / "ck_xt_24_rosetta_linker_scores.csv"

METRICS = [
    "pre_linker_flank_total_energy",
    "pre_linker_total_energy",
    "post_linker_flank_total_energy",
    "post_linker_total_energy",
]


def remap_pdb(old_path: str) -> Path:
    """old 'outputs/rank_x/y.pdb' -> REPO/04_rfdiffusion/output/rank_x/y.pdb

    Pass-list entries under outputs/ck_xt_24_rosetta_linker*_top10_pdbs/ are
    byte-identical 'NN_'-prefixed copies of rank_* designs; canonicalize them
    back to the rank_* source (verified with cmp).
    """
    p = Path(old_path)
    if p.is_absolute() and p.exists():
        return p
    parts = p.parts
    if parts and parts[0] == "outputs":
        parts = parts[1:]
    m = re.match(r"^\d+_(rank_\d+_sample_\d+_\d+)\.pdb$", parts[-1])
    if m:
        stem = m.group(1)
        rank_dir = stem.rsplit("_", 1)[0]
        return REPO / "04_rfdiffusion" / "output" / rank_dir / f"{stem}.pdb"
    return REPO / "04_rfdiffusion" / "output" / Path(*parts)


def parse_score_sc(score_sc: Path) -> dict:
    lines = score_sc.read_text().splitlines()
    header = None
    for line in lines:
        if line.startswith("SCORE:"):
            cols = line.split()[1:]
            if header is None:
                header = cols  # first SCORE: line is the column header
            else:
                return dict(zip(header, cols))  # second SCORE: line is the data row
    raise ValueError(f"no data row in {score_sc}")


def run_one(pdb: Path, jobs_dir: Path) -> dict:
    job = jobs_dir / pdb.stem
    job.mkdir(parents=True, exist_ok=True)
    out_pdb = job / f"{pdb.stem}_0001.pdb"
    score_sc = job / "score.sc"
    cmd = [
        str(ROSETTA),
        "-s", str(pdb),
        "-parser:protocol", str(PROTOCOL),
        "-no_fconfig",
        "-ignore_unrecognized_res",
        "-out:file:scorefile", str(score_sc),
        "-out:path:all", str(job),
        "-overwrite",
    ]
    row = {"pdb": str(pdb.relative_to(REPO)), "name": pdb.stem,
           "status": "PASS", "error": "", "wall_sec": ""}
    try:
        if score_sc.exists() and out_pdb.exists():  # resume: reuse finished job
            scores = parse_score_sc(score_sc)
            row["wall_sec"] = "cached"
        else:
            t0 = time.time()
            log = job / "rosetta.log"
            with log.open("w") as fh:
                proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT)
            row["wall_sec"] = f"{time.time() - t0:.0f}"
            if proc.returncode != 0:
                raise RuntimeError(f"rosetta exit {proc.returncode}")
            scores = parse_score_sc(score_sc)
        for m in METRICS:
            row[m] = scores[m]
        row["rosetta_pdb"] = str(out_pdb.relative_to(REPO)) if out_pdb.exists() else ""
    except Exception as e:  # keep batch going
        row["status"] = "FAIL"
        row["error"] = str(e)[:200]
    return row


def ranks(values):
    order = sorted(range(len(values)), key=lambda i: values[i])
    r = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r


def spearman(xs, ys):
    rx, ry = ranks(xs), ranks(ys)
    n = len(xs)
    mx, my = sum(rx) / n, sum(ry) / n
    cov = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    vx = sum((a - mx) ** 2 for a in rx)
    vy = sum((b - my) ** 2 for b in ry)
    return cov / (vx * vy) ** 0.5 if vx and vy else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=26)
    ap.add_argument("--outdir", type=Path, default=STEP_DIR / "output_relaxed")
    ap.add_argument("--limit", type=int, default=0, help="only first N pdbs (test)")
    ap.add_argument("--pdb", type=Path, default=None, help="single pdb (test)")
    ap.add_argument("--pdb-list", type=Path, default=None,
                    help="text file of pdb paths (absolute or repo-relative)")
    args = ap.parse_args()

    if args.pdb:
        pdbs = [args.pdb.resolve()]
    elif args.pdb_list:
        pdbs = []
        for line in args.pdb_list.read_text().splitlines():
            if not line.strip():
                continue
            p = Path(line.strip())
            pdbs.append(p if p.is_absolute() else (REPO / p).resolve())
    else:
        seen = set()
        pdbs = []
        for line in PASS_LIST.read_text().splitlines():
            if not line.strip():
                continue
            p = remap_pdb(line.strip())
            if p not in seen:  # 12 pass-list entries are top10 copies of rank_* designs
                seen.add(p)
                pdbs.append(p)
        missing = [p for p in pdbs if not p.exists()]
        if missing:
            sys.exit(f"missing {len(missing)} pdbs, e.g. {missing[0]}")
        if args.limit:
            pdbs = pdbs[: args.limit]

    outdir = args.outdir.resolve()  # absolute: rosetta prepends -out:path:all to relative scorefile paths
    jobs_dir = outdir / "jobs"
    jobs_dir.mkdir(parents=True, exist_ok=True)
    print(f"scoring {len(pdbs)} structures with {args.jobs} workers", flush=True)

    rows = []
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = {ex.submit(run_one, p, jobs_dir): p for p in pdbs}
        for i, fut in enumerate(as_completed(futs), 1):
            row = fut.result()
            rows.append(row)
            if i % 10 == 0 or i == len(pdbs):
                print(f"  {i}/{len(pdbs)} done", flush=True)

    rows.sort(key=lambda r: float(r["post_linker_flank_total_energy"]) if r["status"] == "PASS" else 1e18)
    for rank, r in enumerate([r for r in rows if r["status"] == "PASS"], 1):
        r["rank_relaxed"] = rank

    out_csv = outdir / "ck_xt_24_rosetta_linker_scores_relaxed.csv"
    cols = ["rank_relaxed", "name", "pdb", "status", "error", "wall_sec"] + METRICS + [
        "delta_linker_flank_total", "rosetta_pdb"]
    for r in rows:
        if r["status"] == "PASS":
            r["delta_linker_flank_total"] = (
                float(r["post_linker_flank_total_energy"]) - float(r["pre_linker_flank_total_energy"]))
    with out_csv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out_csv}", flush=True)

    # comparison against the old unrelaxed ranking
    old = {}
    if OLD_CSV.exists():
        with OLD_CSV.open() as fh:
            for r in csv.DictReader(fh):
                old[Path(r["pdb"]).stem] = r
    common = [r for r in rows if r["status"] == "PASS" and r["name"] in old]
    if common:
        old_e = [float(old[r["name"]]["linker_flank_total_energy"]) for r in common]
        pre_e = [float(r["pre_linker_flank_total_energy"]) for r in common]
        post_e = [float(r["post_linker_flank_total_energy"]) for r in common]
        print(f"n_common={len(common)}")
        print(f"spearman(old_unrelaxed, pre_relax)  = {spearman(old_e, pre_e):.4f}  (protocol sanity, ~1 expected)")
        print(f"spearman(old_unrelaxed, post_relax) = {spearman(old_e, post_e):.4f}  (rank stability)")
        old_top10 = {n for n, _ in sorted(((n, float(old[n]['linker_flank_total_energy'])) for n in [r['name'] for r in common]), key=lambda t: t[1])[:10]}
        new_top10 = {r["name"] for r in sorted(common, key=lambda r: float(r["post_linker_flank_total_energy"]))[:10]}
        print(f"old top10: {sorted(old_top10)}")
        print(f"new top10: {sorted(new_top10)}")
        print(f"top10 overlap: {len(old_top10 & new_top10)}/10")
    n_fail = sum(1 for r in rows if r["status"] != "PASS")
    print(f"PASS={len(rows) - n_fail} FAIL={n_fail}")


if __name__ == "__main__":
    main()
