# Step 03 — CA-tangent geometric filter (144 → 18)

Purpose: apply a hard geometric screen to the 144 pose samples from step 02, keeping
geometries in which the gap between the HA C-terminus and the scaffold N-terminus is
suitable for a 20 aa rigid helix linker. Outputs 18 candidates for step 04
(RFdiffusion).

## Layout

- `scripts/filter_linker_candidates.py` — custom script. For each sample it computes
  endpoint distances, endpoint angles, tangent bend, C3-axis consistency, etc.,
  filters by thresholds and ranks the survivors.
- `output/`
  - `screened_summary.csv` — screening record for all 144 samples
  - `candidates.csv` — the 18 passing samples (19 lines including the header);
    rank 1 = sample_0108 (z28/350°)
  - `cif/` — CIFs of the 18 passing samples, renamed `rank_NNN_sample_MMMM.cif`
    (this naming convention persists through all later steps)

## How to run / key parameters

```bash
python scripts/filter_linker_candidates.py \
  --samples-dir <02_pose_sampling/output/> \
  --outdir <this directory's output/>
```

Screening logic `--orientation-mode ca-tangent` (the default, and what this run
used): requires the CA524→CA525 direction at the HA end and the CA526→CA527 direction
at the scaffold end to be tangent-continuous. Hard thresholds (script defaults):
endpoint distance 22–38 Å (`--min/max-distance`, target 30 Å), inter-chain distance
spread `--max-distance-spread 1.0`, mean endpoint angle ≤ 70°, single endpoint ≤ 75°,
tangent bend ≤ 25°, C3 axis angle ≤ 2°, axis offset ≤ 1 Å.

**History**: the older `terminal-bond` logic produced unreasonable results, was
deprecated, deleted, and the screen was rerun (see §5 of
`docs/CK_XT_24_rigid_linker_work_summary.md` in the private archive); all 18 current
candidates come from the ca-tangent logic.

## Completion criteria

`candidates.csv` contains exactly 18 passing samples, and the 18 files in `cif/`
correspond one-to-one with the ranks in candidates.csv.
