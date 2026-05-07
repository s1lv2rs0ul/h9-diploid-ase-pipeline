# Paper-published scripts

Scripts in this folder reproduce analyses **already published** in:

> Sawarkar et al. *The diploid reference genome of a human embryonic stem cell line.* bioRxiv 2026.03.26.714432. <https://www.biorxiv.org/content/10.64898/2026.03.26.714432v1>

These are the upstream variant-calling and gene-level ASE counting steps. The methods are described in the paper; this folder ships the actual code so the analysis is reproducible.

## Contents

| # | Script | What it does |
|---|--------|--------------|
| 01 | `01_build_master_variants.py` | CIGAR-parse hap1↔hap2 PAF alignments → het variant catalog (4,515 variants across SNCA, HTT, LRRK2, GBA, SMN1/2). Classifies `in_gene_body` from CAT GTF transcript spans. |
| 02 | `02_count_ase_gene_level.sh` | Count haplotype-specific reads in gene-body intervals across iPSC, FPP, astrocyte, microglia from `.ase.bam` files. Output is the data behind **Fig 3** of the paper. |
| 03 | `03_count_ase_per_variant.py` | Count haplotype reads at every individual het position (per-variant ASE table) — generated as part of the paper data products, not used as a published headline result. |
| — | `lib/variant_analysis_library.py` | Shared helpers: position arithmetic (`rel_to_abs`), gene-body membership tests, samtools wrappers. |

## Output → results/paper_published/

- `master_variants.csv` — the 4,515-variant catalog
- gene-level ASE counts → see `data/small_inputs/paper_published/ase_gene_level_counts.csv`

## How to re-run

These scripts require the raw data described in [`../../data/EXTERNAL_DATA.md`](../../data/EXTERNAL_DATA.md):
- Both haplotype FASTAs (`t2t_h9_v01_hap{1,2}.fa`)
- 8 `.ase.bam` files (4 lineages × 2 haplotypes)
- 4 diploid BAMs
- CAT GTF annotations (`H9_HAP{1,2}.gtf`)

```bash
python 01_build_master_variants.py
bash   02_count_ase_gene_level.sh
python 03_count_ase_per_variant.py
```
