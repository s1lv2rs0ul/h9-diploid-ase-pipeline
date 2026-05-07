# Data

## What's in this folder

```
data/
├── small_inputs/             # < 1 MB files, committed to repo
│   ├── ase_gene_level_counts.csv          # Fig 3 ASE counts (5 genes × 4 lineages)
│   ├── alphagenome_variant_summary.csv    # per-variant max_quantile summary (3,049 variants)
│   ├── het_snvs_from_ucsc_dipcall_vcf.csv # 3,272 het SNVs from UCSC dipcall VCF (hg38 truth)
│   └── snvs_for_alphagenome.csv           # 3,252 SNVs lifted T2T → hg38 for AlphaGenome
└── EXTERNAL_DATA.md           # large files hosted externally (Zenodo / OSF)
```

## File descriptions

### `small_inputs/ase_gene_level_counts.csv`

Gene-level allele-specific expression counts produced by `scripts/05_count_ase_gene_level.sh`.

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

### `small_inputs/alphagenome_variant_summary.csv`

Compact summary of AlphaGenome scores (one row per variant). Generated from the 22 GB raw output by `scripts/10_summarize_alphagenome.py`.

| column | meaning |
|--------|---------|
| our_variant_id | `chr_pos_ref_alt` (hg38 coords) |
| gene | SNCA / HTT / LRRK2 / GBA / SMN1 |
| in_gene_body | True if variant lies within H9_HAP1.gtf transcript span |
| max_quantile | max quantile_score across ALL tissue × assay tracks (saturates ≥ 0.95 — use brain-specific scores for ranking) |
| max_raw | max abs(raw_score) across all tracks |
| n_high_scores | number of tracks where quantile ≥ 0.9 |
| n_scores | total number of tracks scored for this variant |

**Note:** `max_quantile` is a noisy upper bound (saturates for most variants). For ranking, use brain-tissue-filtered scores in `results/alphagenome_brain_scores.csv`.

### `small_inputs/het_snvs_from_ucsc_dipcall_vcf.csv`

Authoritative hg38 het SNV positions from the UCSC dipcall diploid VCF (`hg38.t2t_h9_v01.dip.vcf.gz`). Used as the truth set for cross-validating our minimap2-based liftover. 3,272 SNVs across the 5 focus loci.

### `small_inputs/snvs_for_alphagenome.csv`

Output of `scripts/08_liftover_for_alphagenome.py`. T2T het SNVs lifted to hg38 via per-region minimap2 alignment, with the hg38 REF base looked up from the genome FASTA and `hap_alt` determined per haplotype.

| column | meaning |
|--------|---------|
| variant_id | T2T native variant ID |
| chrom, pos, ref, alt | hg38 coordinates |
| hap_alt | hap1 / hap2 / both_differ (tri-allelic vs hg38) |
| in_gene_body | as above |
| ucsc_validated | ucsc_exact (matches dipcall VCF) / our_only |
| ... | |

Per-gene UCSC dipcall validation rates: SNCA 99 %, HTT 85 %, LRRK2 98 %, GBA 98 %, SMN1 0 % (segdup-driven divergence).

## Large data NOT in repo

See [`EXTERNAL_DATA.md`](EXTERNAL_DATA.md) for genome FASTAs, BAM files, GTF annotations, and the 22 GB raw AlphaGenome score table.
