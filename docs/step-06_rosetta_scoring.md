# Step 06 — Rosetta linker+flank energy scoring → top 10

Purpose: score the linker regions of the step-05 watch-passed backbones with Rosetta,
rank by `linker_flank_total_energy` ascending, and take the top 10 for step 07
ProteinMPNN sequence design.

**The current authoritative ranking is in `output_relaxed/` (FastRelax version,
2026-09-18)**: the old unrelaxed scoring (`output/`) proved to be a clash-dominated
pre-filter; the old and new top-10 overlap by only 1 design. Downstream steps must use
the relaxed ranking. Details in `output_relaxed/comparison_report.md`.

## Relaxed rescoring (`output_relaxed/`, authoritative)

- Protocol `protocol/linker_score_relaxed.xml`: within a single run it records both
  pre-relax (`pre_*`) and post-relax (`post_*`) linker / linker+flank
  TotalEnergyMetric (ref2015, coordinate_constraint weight explicitly set to 0); in
  between, `AddConstraintsToCurrentConformationMover` (CA_only, coord_dev 0.5) +
  `FastRelax` (repeats=1, constraint weight 1.0).
- Run: `python3 scripts/score_relaxed_linkers.py --jobs 26` (reads the step-05
  pass-list, maps old `outputs/` paths onto the reorganized layout, drops the 12
  duplicate top-10 copies → 152 unique structures; finished jobs are skipped, so the
  run is resumable).
- Results: 152 PASS / 0 FAIL; post-relax flank energy −53.6 ~ −1.4 (all negative);
  Spearman(old, post-relax) = 0.625, old/new top-10 overlap 1/10. New top 1 =
  `rank_003_sample_0073_0` (−53.56). Full comparison in
  `output_relaxed/comparison_report.md`.
- Products: `output_relaxed/jobs/<name>/` (relaxed PDB + score.sc),
  `ck_xt_24_rosetta_linker_scores_relaxed.csv`, `comparison_report.md`.
- Known boundary: this is still a backbone + Rosetta-repacked placeholder-side-chain
  level evaluation (scoring precedes MPNN sequence design); the one old top-10 member
  surviving in the new top 10 later failed in AF2 — energy ranking does not replace
  sequence-level validation.

## Old unrelaxed scoring (`output/`, superseded, archive only)

Old protocol `protocol/linker_score.xml`: ref2015 + two `TotalEnergyMetric`s
(`linkers` = 47-66,240-259,433-452; `linker_flanks` = 44-69,237-262,430-455), pure
scoring with no mover; missing side chains were only repacked at import
(`pack_missing_sidechains`). All scores were positive (fa_rep clash-dominated) and
reflect only whether the backbone/junction has obvious clashes, not a stability
ranking.

- `output/jobs/` — 152 job directories (each with `score.sc`)
- `output/ck_xt_24_rosetta_linker_scores.csv` — 152 PASS summary (old ranking; pre-relax
  energies are consistent with the `pre_*` values of the new version, Spearman 0.9975)
- `output/ck_xt_24_rosetta_linker_top10.txt` + `output/top10_pdbs/` — old top 10
  (the actual step-07 MPNN input, produced **before** the relaxed version existed)
- `test_run/` — protocol development / trial-run products (including the full
  `linker_score.out` log)

The old driver script `score_rosetta_linkers.py` was not archived; the single-structure
command format is the same as the relaxed version (see the cmd construction in
`scripts/score_relaxed_linkers.py`). Rosetta was release-427
`rosetta_scripts.default.linuxgccrelease`.

## Completion criteria

`output_relaxed/ck_xt_24_rosetta_linker_scores_relaxed.csv` covers all 152 unique
structures with PASS status; the new top 10 is ranked by
`post_linker_flank_total_energy` ascending.
