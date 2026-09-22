# Step 02 — C3-axis pose sampling

Purpose: on top of the step-01 pre-aligned structure, hold the HA fixed and
exhaustively generate candidate relative geometries by translating the scaffold along
the common C3 axis and rotating it around the axis, for the step-03 geometric screen.
Pure geometric sampling — no screening or design logic.

## Layout

- `scripts/sample_c3_axis.py` — custom script (depends only on numpy). Fits each C3
  axis by SVD → aligns the scaffold axis antiparallel to the HA axis → radial
  phase alignment (chain A onto chain A) → applies z-translation + rotation about the axis.
- `output/`
  - `sample_0001.cif … sample_0144.cif` — 144 mmCIF files (≈ 1.5 MB each, containing
    all atoms of the HA plus the moved scaffold; header comments record
    z_shift/rotation/orientation). Numbering increments with z_shift outermost and
    rotation innermost (sample_0001 = z20/0°).
  - `summary.csv` — per-sample z_shift/rotation and three-chain endpoint distance
    statistics (endpoint mean/min/max and per-chain values, range ≈ 20.2–23.7 Å).

## How to run / key parameters

```bash
python scripts/sample_c3_axis.py \
  --prealigned <01_prealign/output/CK_XT_24_scaffold_prealigned.cif> \
  --outdir <this directory's output/>
```

Parameters used this run (i.e. the script defaults): `--z-shifts 20,24,28,32` (4
levels, Å along the axis) × `--rotations 0:350:10` (36 levels; chain identity is
fixed, so a 120° rotation is not equivalent and the full 360° must be covered) =
**144 samples**; `--x-shifts 0 --y-shifts 0`; `--orientation face-to-face`; phase
alignment on (`--no-phase-align` not used); `--max-samples 1000`; endpoint atoms: HA
C-terminus / scaffold N-terminus. Completed in a single run on 2026-05-18, no retries.

## Completion criteria

`output/` contains 144 sample CIFs + 1 summary.csv (145 lines including the header);
the parameters in each CIF's header comments match the numeric file index.
