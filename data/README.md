# Data

Inputs to the pipeline are partitioned by which project they belong to.

```
data/
├── small_inputs/
│   ├── paper_published/         ← inputs for already-published Sawarkar et al. 2026
│   │   └── ase_gene_level_counts.csv
│   └── followup_inprogress/     ← inputs for unpublished AlphaGenome follow-up
│       ├── snvs_for_alphagenome.csv
│       ├── alphagenome_variant_summary.csv
│       └── het_snvs_from_ucsc_dipcall_vcf.csv
├── EXTERNAL_DATA.md             ← large files hosted externally (Zenodo / OSF / UCSC)
└── README.md                    ← this file
```

## paper_published/

### `ase_gene_level_counts.csv`

Gene-level allele-specific expression counts produced by `scripts/paper_published/02_count_ase_gene_level.sh`. **This is the data behind Fig 3 of the paper.**

| column | meaning |
|--------|---------|
| Gene | SNCA, HTT, LRRK2, GBA1, SMN1, SMN2 |
| CellType | iPSC, FPP, astrocyte, microglia |
| Hap1_ASE | reads aligned to hap1 in haplotype-specific gene-body interval |
| Hap2_ASE | reads aligned to hap2 in haplotype-specific gene-body interval |
| ASE_Total | Hap1_ASE + Hap2_ASE |
| Full_H1 | total reads in hap1 region (no allele filter) |
| Full_H2 | total reads in hap2 region (no allele filter) |
| Full_Total | Full_H1 + Full_H2 |
| Pct_Informative | ASE_Total / Full_Total × 100 |

## followup_inprogress/

### `snvs_for_alphagenome.csv`

Output of `scripts/followup_inprogress/01_liftover_for_alphagenome.py`. T2T het SNVs lifted to hg38 via per-region minimap2 alignment, with hg38 REF base from `samtools faidx` and `hap_alt` determined per haplotype.

Per-gene UCSC dipcall validation rates: SNCA 99 %, HTT 85 %, LRRK2 98 %, GBA 98 %, SMN1 0 % (segdup-driven divergence).

### `alphagenome_variant_summary.csv`

Compact summary of AlphaGenome scores (one row per variant, 3,049 variants). Generated from the 22 GB raw output by `scripts/followup_inprogress/03_summarize_alphagenome.py`.

| column | meaning |
|--------|---------|
| our_variant_id | `chr_pos_ref_alt` (hg38) |
| gene | SNCA / HTT / LRRK2 / GBA / SMN1 |
| in_gene_body | True if variant lies within H9_HAP1.gtf transcript span |
| max_quantile | max quantile_score across all tissue × assay tracks (saturates ≥ 0.95 — use brain-specific scores for ranking) |
| max_raw | max abs(raw_score) across all tracks |
| n_high_scores | tracks where quantile ≥ 0.9 |
| n_scores | total tracks scored for this variant |

### `het_snvs_from_ucsc_dipcall_vcf.csv`

Authoritative hg38 het SNV positions from the UCSC dipcall diploid VCF (`hg38.t2t_h9_v01.dip.vcf.gz`). Used as the truth set for cross-validating the minimap2-based liftover. 3,272 SNVs across the 5 focus loci.

## Large data NOT in repo

See [`EXTERNAL_DATA.md`](EXTERNAL_DATA.md) for genome FASTAs, BAM files, GTF annotations, and the 22 GB raw AlphaGenome score table.
