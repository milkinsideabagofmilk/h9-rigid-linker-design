# H9 rigid-linker design & screening pipeline

Computational pipeline for designing a **20 aa rigid helix linker** that fuses the H9
hemagglutinin (HA) C3 trimer (`CK_XT_24`, residues 1–525) to a C3-symmetric scaffold
trimer (`ChainABC`, renumbered 526–652) in a face-to-face (C-to-N) arrangement,
preserving C3 symmetry. Each HA chain is joined to the corresponding scaffold chain by
one linker, giving a covalently linked, three-chain fusion construct
(193 aa per chain, 579 aa per trimer).

The linker is intended to act as a rigid spacer that holds the scaffold at a defined
distance and orientation below the HA trimer. All screening criteria, thresholds and
counts recorded here reflect the scripts and outputs actually in this repository.

## Pipeline overview

```
CK_XT_24.cif (HA trimer, AlphaFold Server)          ChainABC_518_644.cif (scaffold, PyMOL 3.1.6.1)
        └────────────────── 01 prealign ──────────────────┘
                          (C3 axes collinear, centroid 19 Å)
        02 pose sampling            144 poses (4 z-shifts × 36 rotations)
        03 geometric filter         144 → 18   (ca-tangent hard thresholds)
        04 RFdiffusion              18 × 10 = 180 linker backbones
        05 helix-geometry watch     180 → 164 PASS rows (152 unique structures)
        06 Rosetta relaxed scoring  top 10 by post-FastRelax linker+flank energy
        07 ProteinMPNN              linker-only, C3-tied design, 10 × 2 = 20 sequences
        08 AF2-multimer validation  16/20 ipTM ≥ 0.65 → 10 final candidates
```

Per-step documentation (purpose, exact parameters, commands, completion criteria) is
in [`docs/`](docs): `docs/step-01_prealign.md` through `docs/step-08_af2_validation.md`.

## Repository layout

```
README.md
docs/            per-step documentation (step-01 … step-08)
scripts/         original code, grouped by step (01–04, 06–08)
protocols/       Rosetta protocol XML (linker_score_relaxed.xml)
results/         authoritative result tables (CSV) + ranking comparison report
examples/        small representative inputs/outputs (see below)
```

This repository intentionally contains **representative examples, not full design
outputs** (the complete run produced ~5.7 GB of structures). `examples/` holds:

- `input_structures/` — HA trimer, scaffold trimer, and the pre-aligned complex (mmCIF)
- `pose_samples/` — 2 of 144 enumerated poses
- `geometric_candidates/` — 3 of 18 screened candidate geometries
- `rfdiffusion_outputs/` — one representative generated backbone, the 18-run command list
- `top10_relaxed_backbones/` — the 10 FastRelax-scored backbones that entered sequence design
- `af2_final_candidates/` — AF2-multimer predictions (PDB) for the 10 final candidates,
  plus pLDDT/PAE/coverage plots for the top 2
- `ck_xt_24_mpnn_c3_tied_20.fasta` — the 20 designed trimer sequences validated in step 08

## Input formats

- **Structures**: mmCIF (steps 01–03), PDB (steps 04–08). Step 04 inputs carry
  `REMARK canonicalized_c3_to_z True` and use Rosetta/pose continuous numbering.
- **Sequences**: FASTA, 193 aa per chain (designed sequence already includes a
  C-terminal His6 tag); step 08 concatenates three identical chains with `:` into one
  579-aa multimer query.
- **Constraints**: ProteinMPNN JSONL files (parsed chains, fixed positions, C3-tied
  positions) were generated per run and are described in `docs/step-07_proteinmpnn.md`.

## Dependencies

Versions are taken from run records, not guessed:

- Custom steps 01–03, 08 helper scripts: Python 3 + `numpy` only (stdlib otherwise).
- PyMOL 3.1.6.1 — scaffold structure export (step 01 input provenance).
- RFdiffusion (`run_inference.py`, symmetry config, `Base_ckpt.pt`, T=50,
  `model_runner=SelfConditioning`, `symmetric_self_cond=true`) — step 04.
- Rosetta release-427 (`rosetta_scripts.default.linuxgccrelease`), ref2015 score
  function — steps 06, 08 rescore. Set `ROSETTA_BIN` (or pass `--rosetta`) to a local
  build; Rosetta is **not** redistributed here.
- ProteinMPNN weights `v_48_020` — step 07.
- ColabFold 1.6.1 / AlphaFold2-multimer v3 (`colabfold_batch`, `--num-models 1
  --num-recycle 3 --rank multimer`) — step 08; MSA via the public mmseqs2 server.

Heavy steps (04, 06, 07, 08) ran on a single local RTX 5090 under WSL2. For WSL2,
`--disable-unified-memory` is required for ColabFold (XLA unified-memory spilling
otherwise OOMs the VM; see `docs/step-08_af2_validation.md`).

## How to run

Run steps in order; each stage doc gives the exact command and parameters. Notes:

- The step 01 script's default path constants still point at the original workspace
  layout — pass `--ha/--scaffold/--combined-out` explicitly when rerunning.
- Steps 05 (watch script) and the old unrelaxed step 06 driver were not archived;
  their behavior and commands are documented in the respective step docs, and the
  archived `score_relaxed_linkers.py` supersedes the old driver. Some driver scripts
  still reference the original workspace directory layout and are provided as a
  faithful record of how the runs were executed.
- All results in `results/` were produced by the recorded runs (2026-05 design
  campaign; 2026-09-18/20 relax + rerun); nothing in this repository has been
  re-executed during archiving.

## Results (authoritative, 2026-09-20)

Tables in `results/`:

| File | Content |
| --- | --- |
| `pose_sampling_summary.csv` | 144 poses: z-shift, rotation, per-chain endpoint distances |
| `geometric_candidates.csv` / `geometric_screen_all144.csv` | 18 passing candidates / full 144-row screen record |
| `helix_watch_summary.csv` / `helix_watch_pass_list.txt` | step-05 helix-geometry prefilter (164 PASS rows, 152 unique) |
| `rosetta_relaxed_scores.csv` / `rosetta_relaxed_vs_unrelaxed_comparison.md` | FastRelax rescoring of 152 backbones + old-vs-new ranking comparison |
| `mpnn_linker_sequences.csv` | 20 C3-tied linker sequences (10 backbones × 2) |
| `af2_results_summary.csv` / `af2_structural_checks.csv` | AF2-multimer metrics and structural checks for all 20 designs |
| `final_candidates.csv` | **10 final candidates** (ipTM ≥ 0.65 + structural pass) |
| `final_candidates_rosetta_rescore.csv` | Rosetta rescore of the 16 AF2-passing designs (real sequences) |

Headline numbers:

- Geometric screen: 144 poses → 18 candidates; rank 1 = `sample_0108` (z 28 Å,
  twist 350°). Orientation logic is `ca-tangent` (endpoint CA tangents must be
  nearly collinear); the older `terminal-bond` logic was abandoned and must not be
  revived.
- RFdiffusion: 180/180 backbones generated; sampled-motif RMSD ≈ 0.24–0.43 Å.
- Watch prefilter: 164 PASS rows = 152 unique structures (12 rows were byte-identical
  top-10 copies re-counted by the scanner; dedup is built into the step 06 driver).
- Rosetta relaxed scoring: 152/152 PASS; post-FastRelax linker+flank energy
  −53.6 to −1.4; old unrelaxed ranking agreed only 1/10 with the relaxed top 10, so
  the relaxed ranking is authoritative.
- ProteinMPNN: 20 sequences, linker-only design, C3-tied across chains; sequences are
  E/K/A/L-rich helix patterns (EAAEK-like).
- AF2-multimer: 16/20 pass ipTM ≥ 0.65 (pLDDT 72–84); structural checks pass 10;
  best overall: `05_rank_016_sample_0075_3_sample1` (ipTM 0.72, flankE −224.1) and
  `03_rank_004_sample_0109_5_sample1` (ipTM 0.71, flankE −224.0); best structural
  fidelity: `04_rank_003_sample_0073_1_sample1` (core RMSD 0.75 Å).

## Numbering conventions (a historical source of bugs)

- RFdiffusion/Rosetta outputs use pose continuous numbering: linker =
  A:47-66, B:240-259, C:433-452.
- ProteinMPNN and per-chain analyses use within-chain numbering: every chain is
  193 aa and its linker is 47–66. Constraint JSONL files must use within-chain
  residue numbers.

## Limitations

- All step 06 Rosetta scores are backbone-level evaluations with placeholder/
  repacked side chains, computed **before** sequence design; they are a pre-filter,
  not a stability verdict. The step 08 rescore (real sequences) is the comparable one.
- The construct validated in step 08 contains only the HA 480–525 fragment (46 aa),
  not full-length HA. AF2 predicts this fragment poorly in isolation (pLDDT ≈ 47,
  ~12 Å deviation from native) and it weakly drags the linker terminus; the hard
  fidelity threshold is therefore core RMSD ≤ 2.5 Å on linker+scaffold, with linker
  terminal displacement reported but not thresholded.
- Energetic ranking does not replace sequence-level validation, and no experimental
  validation (expression, SEC, EM) exists yet; in particular, EM particle images
  would not by themselves prove linker rigidity or correct assembly.
- Next scientific step: full-length validation of top candidates
  (HA 1–525 + linker + scaffold, ~672 aa/chain, ~2016 aa trimer) by AF2-multimer
  or experiment.

## Provenance & licensing

Original computation May 2026 in a personal RFdiffusion workspace; archived and
documented September 2026. Third-party tools and model weights (RFdiffusion,
Rosetta, ProteinMPNN, ColabFold/AlphaFold) are **not** redistributed in this
repository and remain under their own licenses; obtain them from their respective
sources. Code under `scripts/` is original to this project except where noted in the
step documentation.
