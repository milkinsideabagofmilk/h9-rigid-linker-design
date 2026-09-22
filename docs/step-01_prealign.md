# Step 01 — C3-axis pre-alignment (prealign)

Purpose: hold the HA trimer (`CK_XT_24`, residues 1–525) fixed and apply a rigid-body
transform to the ChainABC scaffold trimer so that the two C3 axes are collinear in a
face-to-face arrangement (HA C-terminus facing the scaffold N-terminus), leaving room
for insertion of a 20 aa rigid linker. This is the first step of the pipeline; its
output is the direct input of step 02 (`02_pose_sampling`).

## Layout

- `scripts/prealign_ck_xt_24_to_scaffold.py` — custom script. Fits both C3 axes by
  SVD, aligns the scaffold axis antiparallel to the HA axis, and scans twist to pick
  the optimum.
- `input/`
  - `CK_XT_24.cif` — original HA C3 trimer (AlphaFold Server output); main input `--ha`
  - `ChainABC_518_644.cif` — scaffold trimer (exported with PyMOL 3.1.6.1); main input `--scaffold`
  - `ChainABC_519_644.cif` — scaffold variant (starting at 519); alternate default for the step-02 sampling script
  - `CK_XT_24_1_525.cif` — HA 1–525 truncation (a by-product of this step, `--truncated-out`; also an alternate input for step 02)
  - `CKXT24_GS_FUSION.cif` — legacy GS-linker fusion reference structure; not an input of the current run
- `output/CK_XT_24_scaffold_prealigned.cif` — pre-aligned merged structure: scaffold
  renumbered from 526, centroid distance 19 Å, twist 353°
- `failure/ChainABC_1_516.cif` — leftover from a failed early scaffold truncation
  (1–516, 2026-05-14); archive only, do not use as input

## How to run / key parameters

```bash
python scripts/prealign_ck_xt_24_to_scaffold.py   # defaults are exactly what this run used
```

Script defaults (i.e. the values actually used this run): `--chains A,B,C`,
`--ha-end-residue 525`, `--target-distance 19.0` (Å, centroid distance between the two
trimers), `--twist-step-deg 1.0` (twist scan step around the axis),
`--renumber-scaffold-start 526`. Note: the path constants hardcoded in the script still
point at the pre-reorganization layout (`REPO_DIR`/`ALIGN_DIR`); when rerunning, pass
`--ha/--scaffold/--combined-out` explicitly.

## Completion criteria

In `output/CK_XT_24_scaffold_prealigned.cif` the three scaffold chains are renumbered
526–652, the two C3 axes are collinear (axis angle ≈ 0°), the centroid distance is
≈ 19 Å, and the orientation is face-to-face.
