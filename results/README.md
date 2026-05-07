# Results

Pipeline outputs are partitioned by project:

```
results/
├── paper_published/           ← outputs that map to Sawarkar et al. 2026
│   └── master_variants.csv
└── followup_inprogress/       ← outputs of the unpublished AlphaGenome follow-up
    ├── integrated_variant_analysis.csv
    ├── snvs_alphagenome_scored.csv
    ├── alphagenome_brain_scores.csv
    ├── alphagenome_top_effects.csv
    └── tf_differential_binding.csv
```

## paper_published/

### `master_variants.csv` (7.2 MB)

Single source of truth: 4,515 het variants between H9 hap1 and hap2 across 5 genes (SNCA, HTT, LRRK2, GBA, SMN1/2). Native T2T coordinates, PAF-verified positions, gene-body classification from CAT GTF, and per-variant ASE read counts across 4 lineages × 2 haplotypes. Includes Evo2 context sequences (±512 bp) for downstream scoring.

Built by: `scripts/paper_published/01_build_master_variants.py` + `03_count_ase_per_variant.py`.

## followup_inprogress/

### `integrated_variant_analysis.csv` (6.3 MB) — **headline output**

Priority-ranked variant table combining AlphaGenome scores + ASE counts + VEP/CADD annotations + gene-body residence into a per-variant `priority_score`.

Built by: `scripts/followup_inprogress/04_integrate_ase_alphagenome.py` (script to be reconstructed; output is preserved here from a prior run on 2026-03-18).

Top SNCA candidate: `chr4:89,836,694 C>T` (gene body, AlphaGenome astrocyte q = 0.9995).

### `snvs_alphagenome_scored.csv` (710 KB)

Intermediate file: 2,717 SNVs joined with their AlphaGenome variant-summary scores (max_quantile, max_raw, n_high_scores) before the ASE integration step.

### `alphagenome_brain_scores.csv` (9.8 MB)

Brain/neural-tissue-filtered subset of AlphaGenome scores: 30,490 rows across 3,049 variants × ~10 brain-related tracks each (UBERON brain regions + neural cell types). Use this for tissue-aware ranking instead of the saturated `max_quantile`.

### `alphagenome_top_effects.csv` (4.8 MB)

Tracks with `quantile_score ≥ 0.9` from the full AlphaGenome output, kept as the top 5 per variant.

### `tf_differential_binding.csv` (190 KB)

Predicted TF-binding-motif disruption between hap1 and hap2: 4,078 entries from FIMO scanning against JASPAR 2024.

## Files NOT included

- **Raw AlphaGenome score table** (`alphagenome_scores.csv` / `.csv.gz`, 22 GB raw / 2 GB gzipped) — too large for git. Will be deposited at Zenodo when finalized; see [`../data/EXTERNAL_DATA.md`](../data/EXTERNAL_DATA.md). The brain-filtered and top-effects subsets here are the analysis-ready summaries derived from it.
