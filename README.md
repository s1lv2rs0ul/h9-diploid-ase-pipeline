# H9 Diploid ASE Pipeline

**Allele-specific expression pipeline using the diploid reference genome of an embryonic stem cell line**

Companion repository to:

> Sawarkar et al. *The diploid reference genome of a human embryonic stem cell line.* bioRxiv 2026.03.26.714432. <https://www.biorxiv.org/content/10.64898/2026.03.26.714432v1>

This repository contains the scripts, small input files, and final result tables/figures needed to reproduce the **allele-specific expression (ASE) and non-coding regulatory variant prioritization analyses** at five neurodegeneration risk loci (SNCA, HTT, LRRK2, GBA1, SMN1/2) in H9 human embryonic stem cells across four lineages (iPSC, floor-plate progenitors, astrocytes, microglia), built on the haplotype-resolved T2T H9 diploid assembly (`t2t_h9_v01`).

---

## What the pipeline does

```
T2T H9 diploid genome (hap1 + hap2)
        │
        ├── extract gene regions (PAF-verified offsets)
        ├── hap1 vs hap2 alignment (Winnowmap2)
        └── CIGAR-parsed het variant catalog (4,515 variants)
                │
                ├── per-haplotype RNA-seq (UCSC .ase.bam, 4 lineages)
                │   └── ASE counts: gene level (Fig 3) + per variant
                │
                └── hg38 liftover (minimap2 + UCSC dipcall validation)
                        └── AlphaGenome scoring (2,631 SNVs × 6 modalities)
                                └── Integrated priority ranking (ASE + AG)
                                        └── Top candidate per gene
```

**Headline result:** SNCA shows astrocyte-specific hap2 expression bias (23/56 reads). The AlphaGenome-prioritized non-coding candidate explaining this bias is **`chr4:89,836,694 C>T`** (SNCA gene body, AlphaGenome astrocyte total-RNA-seq quantile = 0.9995, predicted up-regulation).

---

## Repository layout

```
.
├── scripts/                    # numbered pipeline scripts (run in order)
│   ├── 03_build_master_variants.py
│   ├── 05_count_ase_gene_level.sh
│   ├── 06_count_ase_per_variant.py
│   ├── 08_liftover_for_alphagenome.py
│   ├── 09_score_alphagenome.py
│   ├── 10_summarize_alphagenome.py
│   ├── 11_integrate_ase_alphagenome.py    # TO BE REBUILT (see below)
│   ├── 12_fig_alphagenome_results.py
│   ├── 13_fig_ucsc_tracks.py
│   ├── 14_fig_variant_annotation.py
│   ├── 15_fig_pipeline.py
│   └── lib/
│       └── variant_analysis_library.py    # shared helpers
│
├── data/
│   ├── small_inputs/                       # everything < 1 MB ships in repo
│   │   ├── ase_gene_level_counts.csv
│   │   ├── alphagenome_variant_summary.csv
│   │   ├── het_snvs_from_ucsc_dipcall_vcf.csv
│   │   └── snvs_for_alphagenome.csv
│   ├── EXTERNAL_DATA.md                    # Zenodo / OSF DOIs (TBD)
│   └── README.md                           # describes every input file
│
├── results/                                # final tables (committed for inspection)
│   ├── master_variants.csv
│   ├── integrated_variant_analysis.csv     # the priority-ranked output
│   ├── snvs_alphagenome_scored.csv
│   ├── alphagenome_brain_scores.csv
│   ├── alphagenome_top_effects.csv
│   └── tf_differential_binding.csv
│
├── figures/                                # figures (PDF + PNG)
│   ├── Fig3_ASE_combined.{pdf,png}         # published — paper Fig 3
│   ├── scratch/                            # rough first-pass figures (not publication-ready)
│   └── README.md
├── docs/                                   # GitHub Pages source (optional)
├── environment.yml                         # conda env
├── CITATION.cff
└── LICENSE                                 # MIT
```

## Pipeline steps (numbered)

| # | Script | What it does |
|---|--------|--------------|
| 01–02 | (manual, requires raw data) | Extract syntenic gene regions per haplotype with `samtools faidx`, align hap1↔hap2 with Winnowmap2 (`asm5`). |
| 03 | `03_build_master_variants.py` | CIGAR-parse PAFs → het variant catalog (4,515 variants), classify in_gene_body using CAT GTF. |
| 04 | (validation) | Spot-check FASTA bases at variant positions. |
| 05 | `05_count_ase_gene_level.sh` | Count haplotype-specific reads in gene bodies across 4 lineages (`samtools view -c`). |
| 06 | `06_count_ase_per_variant.py` | Per-variant haplotype read counts at every het position. |
| 07 | (manual) | Hap1↔hap2↔hg38 alignment via minimap2. |
| 08 | `08_liftover_for_alphagenome.py` | Lift T2T variants to hg38, cross-validate against UCSC dipcall VCF. |
| 09 | `09_score_alphagenome.py` | AlphaGenome API (1 Mb context, 6 modalities). Requires `ALPHAGENOME_API_KEY`. |
| 10 | `10_summarize_alphagenome.py` | Reduce 22 GB raw scores to brain-filtered + top-effects subsets. |
| 11 | `11_integrate_ase_alphagenome.py` | **TO BE REBUILT** — joins AG scores + ASE counts → priority ranking. The output (`integrated_variant_analysis.csv`) exists and ships here; the build script needs reconstruction from its columns. |
| 12–15 | `1*_fig_*.py` | Figure scripts (ASE, AlphaGenome, UCSC tracks, variant annotation, pipeline diagram). |

## Quick start

### 1. Clone and set up environment

```bash
git clone https://github.com/<your-username>/h9-diploid-ase-pipeline.git
cd h9-diploid-ase-pipeline
conda env create -f environment.yml
conda activate h9-ase
```

### 2. Inspect committed results

The final tables are already in `results/` — you can explore them without rerunning anything:

```bash
head results/integrated_variant_analysis.csv
```

### 3. Re-run figures from committed results

```bash
python scripts/15_fig_pipeline.py
python scripts/12_fig_alphagenome_results.py
```

### 4. Re-run the full pipeline from raw data

You will need:

- **T2T H9 diploid FASTAs** (`t2t_h9_v01_hap{1,2}.fa`) — UCSC, see `data/EXTERNAL_DATA.md`
- **Haplotype-phased RNA-seq BAMs** (`.ase.bam`, 8 files) — UCSC Genome Server Group
- **CAT GTF annotations** (`H9_HAP{1,2}.gtf`) — UCSC
- **UCSC dipcall diploid VCF** (`hg38.t2t_h9_v01.dip.vcf.gz`) — UCSC
- **AlphaGenome API key** — <https://deepmind.google/alphagenome>

See [`data/EXTERNAL_DATA.md`](data/EXTERNAL_DATA.md) for download links and DOIs.

Then run the numbered scripts in order. **Note:** several scripts currently use hard-coded paths (e.g. `/Users/.../reganalysis/`); see open issues for the path-portability refactor.

## Citation

If you use this pipeline, please cite the H9 T2T diploid reference paper (see `CITATION.cff`).

## License

[MIT](LICENSE)

## Contact

Issues and questions: please open a GitHub issue.
