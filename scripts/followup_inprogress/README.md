# Follow-up scripts (in progress, unpublished)

Scripts in this folder are part of the **unpublished follow-up project**: integrating AlphaGenome predictions with the published gene-level ASE results to prioritize candidate causal non-coding regulatory variants at the five neurodegeneration risk loci.

These analyses are **not yet finalized**. Layout, methodology, and outputs may change.

## Goal

Take the variants and ASE evidence generated for the paper, score them with AlphaGenome (and eventually Evo2), and produce a **ranked shortlist of candidate non-coding regulatory variants** that likely drive the haplotype-specific expression observed at the gene level.

## Contents

| # | Script | What it does | Status |
|---|--------|-------------|--------|
| 01 | `01_liftover_for_alphagenome.py` | Lift T2T het SNVs to hg38 via per-region minimap2 alignment; cross-validate against UCSC dipcall VCF. | Working |
| 02 | `02_score_alphagenome.py` | Submit each SNV to the AlphaGenome API (1 Mb context, 6 modalities). Requires `ALPHAGENOME_API_KEY`. | Working — already run for 2,631 SNVs |
| 03 | `03_summarize_alphagenome.py` | Reduce 22 GB raw scores → per-variant summary, brain-tissue subset, top-effects subset. | Working |
| 04 | `04_integrate_ase_alphagenome.py` | Join AG scores + ASE counts → priority ranking. | **Placeholder — needs to be reconstructed** (see docstring; output `results/followup_inprogress/integrated_variant_analysis.csv` exists from a prior run) |
| 05 | `05_fig_alphagenome_results.py` | 4-panel AlphaGenome overview figure. | Scratch — figure is rough |
| 06 | `06_fig_ucsc_tracks.py` | UCSC-browser-style variant landscape plots. | Scratch |
| 07 | `07_fig_variant_annotation.py` | CADD + consequence-class + GWAS-hit map. | Scratch — most polished of the three |
| 08 | `08_fig_pipeline.py` | Methods/pipeline schematic. | Scratch |

## Outputs → results/followup_inprogress/

- `integrated_variant_analysis.csv` — the headline ranked shortlist (preserved from a prior run; rebuild script pending)
- `snvs_alphagenome_scored.csv` — intermediate AG-only join (no ASE yet)
- `alphagenome_brain_scores.csv` — brain/neural-tissue–filtered scores
- `alphagenome_top_effects.csv` — quantile ≥ 0.9 subset
- `tf_differential_binding.csv` — FIMO TF-disruption predictions

## What's missing / known issues

1. **`04_integrate_ase_alphagenome.py` is a placeholder** — the original build script was lost. The output table is preserved under `results/followup_inprogress/integrated_variant_analysis.csv`; the rebuild needs to reproduce the same columns and `priority_score` formula.
2. **Hard-coded paths** — several scripts reference absolute paths under `~/reganalysis/`. Path-portability refactor pending.
3. **Evo2 scoring not yet run** — sequences are extracted (in `master_variants.csv`) but scoring needs GPU. Will be added as a parallel scoring step alongside AlphaGenome.
4. **SMN1 hg38 liftover gap** — only 81/573 SMN1 SNVs have hg38 coordinates because of segmental-duplication divergence between our minimap2 alignment and UCSC dipcall. Either deferred or needs a SMN1-specific alignment fix.
