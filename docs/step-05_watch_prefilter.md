# Step 05 — linker helix geometry watch pre-filter (180 → 164)

Purpose: run a per-linker helix geometry check on the 180 RFdiffusion backbones from
step 04, removing designs with chain breaks, bends, low helix content, or clashes.
Designs passing on all three chains proceed to step 06 Rosetta scoring.

## Layout

- `output/`
  - `ck_xt_24_rigid_linker_watch.csv` — per PDB × per chain linker metrics:
    helix_ca_fraction, rama_helix_fraction, bend_angle, chain break, <2 Å clash
  - `ck_xt_24_rigid_linker_watch_summary.csv` — summary: 164 PASS / 28 FAIL (per
    chain; failure causes: break×12, bent×10, break|bent×5, low_ca_helix|bent×1)
  - `ck_xt_24_pass_all_rigid_linkers.txt` — list of 164 PDB paths passing on all
    three linkers (the input list for step 06)
  - `deprecated_top30/` — `ck_xt_24_structural_prefilter_top30.{csv,txt}`: a
    discontinued pure-geometry top-30 interim scheme, **do not use**; the project
    switched to Rosetta scoring for top-10 selection (step 06)

## How to run / key parameters

The watch script `watch_rigid_linkers.py` was not archived into this repository (it
lived in the RFdiffusion workspace `scripts/`; usage is described in §6 of the work
summary in the private archive). The command used at the time:

```bash
python scripts/watch_rigid_linkers.py \
  --root outputs \
  --chain-ranges A:47-66,B:240-259,C:433-452 \
  --out outputs/ck_xt_24_rigid_linker_watch.csv \
  --summary-out outputs/ck_xt_24_rigid_linker_watch_summary.csv \
  --pass-list outputs/ck_xt_24_pass_all_rigid_linkers.txt \
  --poll-seconds 30 --stable-seconds 2
```

Linker residue ranges (pose numbering): A:47-66, B:240-259, C:433-452.

**Note**: 12 of the 164 pass-list entries point to `NN_`-prefixed copies under
`outputs/ck_xt_24_rosetta_linker*_top10_pdbs/` — they were double-counted after the
top-10 copies were placed into the watch scan root. `cmp` verified the copies are
byte-identical to the `rank_*` originals, so the **number of unique structures = 152**
(matching the original Rosetta scoring coverage; there are no "12 unscored new
structures"). The step-06 relaxed scoring script has this dedup built in.

## Completion criteria

Pass-list line count = number of three-chain-PASS PDBs in the watch summary (164);
FAIL reason labels are tallied in the summary.
