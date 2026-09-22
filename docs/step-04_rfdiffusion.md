# Step 04 — RFdiffusion generation of 20 aa linker backbones (18 × 10 = 180)

Purpose: on the 18 candidate geometries, use RFdiffusion C3-symmetric motif
scaffolding to generate 20 aa rigid helix linker backbones connecting HA (A480-525)
to the scaffold (A526-652). 10 designs per candidate, 180 in total.

## Layout

- `scripts/prepare_rfdiffusion_candidates.py` — input preparation: CIF → PDB,
  canonicalize the C3 axis onto the Z axis (`REMARK canonicalized_c3_to_z True` in
  the PDB header), generate contigs and a command draft. Key defaults:
  `--contig-ha-start 480`, `--linker-spec 20-20`, `--num-designs 100` (10 were
  actually run, see below).
- `input/`
  - `rank_001…018_*.pdb` — 18 canonicalized input PDBs
  - `index.csv` — candidate → input file / command mapping
  - `rfdiffusion_commands.sh` — the 18 commands actually executed this run
  - `ck_xt_24_20aa_copy/` — another copy of the same 18 PDBs (originally
    RFdiffusion `inputs/ck_xt_24_20aa/`); duplicates of the PDBs above
- `output/`
  - `rank_001_sample_0108/ … rank_018_sample_0039/` — 18 output groups, each with
    `*_0.pdb…*_9.pdb` + same-named `.trb` (config/contig metadata) + `traj/`
    (Xt-1 and pX0 trajectories)
  - `hydra_runs/` — 18 hydra run directories (originally
    `2026-05-20/<timestamp>/`, containing `.hydra/` configs and `run_inference.log`)

## How to run / key parameters

Inside the RFdiffusion installation (this run: a local RFdiffusion checkout, RTX
5090), run once per input (all 18 commands in `input/rfdiffusion_commands.sh`):

```bash
./scripts/run_inference.py --config-name=symmetry \
  inference.symmetry="C3" \
  inference.input_pdb=inputs/ck_xt_24_20aa/rank_001_sample_0108.pdb \
  inference.output_prefix=outputs/rank_001_sample_0108/rank_001_sample_0108 \
  inference.num_designs=10 \
  'contigmap.contigs=[A480-525/20-20/A526-652/0 B480-525/20-20/B526-652/0 C480-525/20-20/C526-652/0]'
```

Other settings (recorded in the hydra config.yaml): checkpoint `Base_ckpt.pt`, T=50
steps, `model_runner=SelfConditioning`, `symmetric_self_cond=true`, `final_step=1`.
Executed sequentially 2026-05-20 11:41 → 16:23, ≈ 16 min per candidate, 180/180
generated; sampled-motif RMSD ≈ 0.24–0.43.

**Numbering**: in the output pose numbering the linker is A:47-66 / B:240-259 /
C:433-452 (46+20+127 = 193 aa per chain). Rework history: the contig start was once
mistakenly written as A1-525; it was corrected to A480-525 and the old products were
removed.

## Completion criteria

Each of the 18 `output/rank_*/` directories contains 10 PDBs + 10 TRBs; the step-05
watch pre-filters exactly these 180 PDBs.
