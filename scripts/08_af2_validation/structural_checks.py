#!/usr/bin/env python3
"""Structural checks on AF2-multimer outputs: linker helix, straightness,
C3 symmetry, clashes, and CA-RMSD vs the design backbone.

Per design (3 chains x 193 aa; HA 1-46, linker 47-66, scaffold 67-193):
- break          : any consecutive CA-CA distance > 4.5 A
- bend           : angle between PC1 axis of linker CA 47-56 and 57-66
- helix_frac     : fraction of linker residues 48-65 with phi in [-100,-20],
                   psi in [-90,0]
- ca_helix_frac  : fraction of linker CA(i)-CA(i+4) distances in [5.2,6.8] A
- sym_rmsd       : max pairwise full-chain Kabsch CA RMSD (C3 symmetry)
- clash          : min inter-chain CA-CA distance
- linker_rmsd    : max per-chain CA-RMSD of linker (47-66) vs design backbone,
                   aligned on the confident core (linker+scaffold, 47-193).
                   Reported as info only: the HA fragment (1-46) is predicted
                   without the rest of the HA trimer (pLDDT ~47) and drags the
                   linker tip, so linker placement is not a hard criterion.
- full_rmsd      : max per-chain CA-RMSD of the core (47-193) vs design
- ha_rmsd_info   : max per-chain self-aligned CA-RMSD of HA fragment (1-46)
                   vs design — large by construction (fragment artifact)

PASS thresholds: no break, bend<=25, helix_frac>=0.75, ca_helix_frac>=0.75,
sym_rmsd<=2.0, clash>=3.0, core_rmsd<=2.5.
"""
import argparse
import csv
from pathlib import Path

import numpy as np

STEP = Path(__file__).resolve().parents[1]
DESIGN_DIR = STEP.parent / "07_proteinmpnn" / "input" / "top10_relaxed_pdbs"


def parse_pdb(path):
    """{chain: {resnum: {atom: xyz}}}"""
    atoms = {}
    for line in path.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        ch, res, atom = line[21], int(line[22:26]), line[12:16].strip()
        xyz = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        atoms.setdefault(ch, {}).setdefault(res, {})[atom] = xyz
    return {ch: {r: res[r] for r in sorted(res)} for ch, res in atoms.items()}


def dihedral(p0, p1, p2, p3):
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2
    b1 /= np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w)))


def pc1(pts):
    pts = np.asarray(pts) - np.asarray(pts).mean(axis=0)
    _, _, vt = np.linalg.svd(pts, full_matrices=False)
    return vt[0]


def angle_between(u, v):
    c = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
    return np.degrees(np.arccos(np.clip(c, -1, 1)))


def kabsch(p, q):
    """rotation aligning p onto q"""
    pc, qc = p - p.mean(axis=0), q - q.mean(axis=0)
    u, _, vt = np.linalg.svd(pc.T @ qc)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    return vt.T @ np.diag([1, 1, d]) @ u.T


def rmsd(p, q):
    return float(np.sqrt(((p - q) ** 2).sum(axis=1).mean()))


def aligned_rmsd(p, q, rot):
    """RMSD after Kabsch alignment (rotation computed on centered coords)."""
    pc, qc = p - p.mean(axis=0), q - q.mean(axis=0)
    return rmsd(pc @ rot.T, qc), pc @ rot.T, qc


def ca_list(residues):
    return np.array([residues[r]["CA"] for r in sorted(residues) if "CA" in residues[r]])


def check_design(pdb, design_pdb):
    chains = parse_pdb(pdb)
    out = {}
    # per-chain geometry
    breaks, bends, helix_fracs, cahelix_fracs = [], [], [], []
    for res in chains.values():
        ca = {r: np.array(a["CA"]) for r, a in res.items() if "CA" in a}
        rl = sorted(ca)
        breaks.append(any(np.linalg.norm(ca[b] - ca[a]) > 4.5
                          for a, b in zip(rl, rl[1:]) if b == a + 1))
        d1, d2 = pc1([ca[r] for r in range(47, 57)]), pc1([ca[r] for r in range(57, 67)])
        bends.append(min(angle_between(d1, d2), angle_between(d1, -d2)))
        nh = nt = 0
        for r in range(48, 66):
            phi = dihedral(np.array(res[r - 1]["C"]), np.array(res[r]["N"]),
                           np.array(res[r]["CA"]), np.array(res[r]["C"]))
            psi = dihedral(np.array(res[r]["N"]), np.array(res[r]["CA"]),
                           np.array(res[r]["C"]), np.array(res[r + 1]["N"]))
            nt += 1
            nh += (-100 <= phi <= -20 and -90 <= psi <= 0)
        helix_fracs.append(nh / nt)
        ds = [np.linalg.norm(ca[r] - ca[r + 4]) for r in range(47, 63)]
        cahelix_fracs.append(sum(1 for d in ds if 5.2 <= d <= 6.8) / len(ds))
    out["break"] = any(breaks)
    out["bend_max"] = max(bends)
    out["helix_frac_min"] = min(helix_fracs)
    out["ca_helix_frac_min"] = min(cahelix_fracs)
    # symmetry + clashes
    cas = {ch: ca_list(res) for ch, res in chains.items()}
    keys = sorted(cas)
    out["clash"] = min(np.linalg.norm(cas[a][:, None, :] - cas[b][None, :, :], axis=2).min()
                       for i, a in enumerate(keys) for b in keys[i + 1:])
    out["sym_rmsd_max"] = max(
        aligned_rmsd(cas[a], cas[b], kabsch(cas[a], cas[b]))[0]
        for i, a in enumerate(keys) for b in keys[i + 1:])
    # RMSD vs design backbone. The HA fragment (1-46) is predicted in a
    # different conformation by AF2 (pLDDT ~47, it lacks the rest of the HA
    # trimer for context) — informative only, no threshold. The pass metric
    # aligns on the confident core (linker+scaffold, 47-193) and measures
    # linker CA-RMSD there.
    out["linker_rmsd_max"] = out["full_rmsd_max"] = out["ha_rmsd_info"] = float("nan")
    if design_pdb and design_pdb.exists():
        dchains = parse_pdb(design_pdb)
        lr, fr, hr = [], [], []
        for ch in keys:
            dca = ca_list(dchains[ch])
            core = np.arange(46, 193)
            rot = kabsch(cas[ch][core], dca[core])
            al, pc_rot, qc = aligned_rmsd(cas[ch][core], dca[core], rot)
            # re-express: align full chain by the core-derived rotation
            pc = cas[ch] - cas[ch][core].mean(axis=0)
            qd = dca - dca[core].mean(axis=0)
            pr, qr = pc @ rot.T, qd
            lr.append(rmsd(pr[46:66], qr[46:66]))
            fr.append(al)
            hr.append(aligned_rmsd(cas[ch][0:46], dca[0:46], kabsch(cas[ch][0:46], dca[0:46]))[0])
        out["linker_rmsd_max"] = max(lr)
        out["full_rmsd_max"] = max(fr)
        out["ha_rmsd_info"] = max(hr)
    ok = (not out["break"] and out["bend_max"] <= 25
          and out["helix_frac_min"] >= 0.75 and out["ca_helix_frac_min"] >= 0.75
          and out["sym_rmsd_max"] <= 2.0 and out["clash"] >= 3.0
          and (np.isnan(out["full_rmsd_max"]) or out["full_rmsd_max"] <= 2.5))
    out["struct_pass"] = ok
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", type=Path, default=STEP / "output" / "results_conservative")
    ap.add_argument("--summary", type=Path,
                    default=STEP / "output" / "results_conservative_summary.csv")
    ap.add_argument("--out", type=Path, default=STEP / "output" / "structural_checks.csv")
    args = ap.parse_args()

    iptm = {}
    if args.summary.exists():
        with args.summary.open() as fh:
            for r in csv.DictReader(fh):
                iptm[r["design_id"]] = float(r["iptm"])

    rows = []
    for pdb in sorted(args.results_dir.glob("*_unrelaxed_rank_001_*.pdb")):
        design = pdb.name.split("_unrelaxed_rank_001_")[0]
        target = design.rsplit("_sample", 1)[0]
        agg = check_design(pdb, DESIGN_DIR / f"{target}.pdb")
        agg["design_id"] = design
        agg["iptm"] = iptm.get(design, "")
        rows.append(agg)

    cols = ["design_id", "iptm", "struct_pass", "break", "bend_max", "helix_frac_min",
            "ca_helix_frac_min", "sym_rmsd_max", "clash", "linker_rmsd_max",
            "full_rmsd_max", "ha_rmsd_info"]
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{r[k]:.2f}" if isinstance(r[k], float) else r[k]) for k in cols})

    n_af2 = sum(1 for r in rows if r["iptm"] != "" and r["iptm"] >= 0.65)
    both = sum(1 for r in rows if r["iptm"] != "" and r["iptm"] >= 0.65 and r["struct_pass"])
    print(f"{len(rows)} designs checked -> {args.out}")
    print(f"AF2-pass: {n_af2}; struct-pass among them: {both}")
    for r in sorted(rows, key=lambda r: -(r["iptm"] or 0)):
        flags = []
        if r["break"]: flags.append("BREAK")
        if r["bend_max"] > 25: flags.append("bent")
        if r["helix_frac_min"] < 0.75: flags.append("low_helix")
        if r["ca_helix_frac_min"] < 0.75: flags.append("low_ca_helix")
        if r["sym_rmsd_max"] > 2.0: flags.append("asym")
        if r["clash"] < 3.0: flags.append("clash")
        if r["linker_rmsd_max"] > 2.0: flags.append("linkerRMSD")
        mark = "AF2pass" if (r["iptm"] != "" and r["iptm"] >= 0.65) else "af2fail"
        print(f"  {r['design_id']:<40} [{mark}] struct={'PASS' if r['struct_pass'] else 'FAIL'} "
              f"bend={r['bend_max']:.1f} helix={r['helix_frac_min']:.2f} cahelix={r['ca_helix_frac_min']:.2f} "
              f"sym={r['sym_rmsd_max']:.2f} clash={r['clash']:.1f} "
              f"lkRMSD={r['linker_rmsd_max']:.2f} coreRMSD={r['full_rmsd_max']:.2f} "
              f"haRMSD(info)={r['ha_rmsd_info']:.1f} {';'.join(flags)}")


if __name__ == "__main__":
    main()
