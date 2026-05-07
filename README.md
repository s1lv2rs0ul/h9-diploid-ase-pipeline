# H9 Diploid ASE Pipeline

**Allele-specific expression pipeline using the diploid reference genome of an embryonic stem cell line**

Companion repository to:

> Sawarkar et al. *The diploid reference genome of a human embryonic stem cell line.* bioRxiv 2026.03.26.714432. <https://www.biorxiv.org/content/10.64898/2026.03.26.714432v1>

This repository contains the scripts, small input files, and result tables/figures for **two related projects** at five neurodegeneration risk loci (SNCA, HTT, LRRK2, GBA1, SMN1/2) in H9 human embryonic stem cells across four lineages (iPSC, floor-plate progenitors, astrocytes, microglia), built on the haplotype-resolved T2T H9 diploid assembly (`t2t_h9_v01`):

1. **`paper_published/`** — variant catalog and gene-level ASE analysis as reported in the Sawarkar et al. paper.
2. **`followup_inprogress/`** — unpublished follow-up: AlphaGenome-based non-coding regulatory variant prioritization, building on the paper's data.

Every script, data file, result, and figure folder is partitioned into these two sections so the boundary between published and in-progress work is explicit.

---

## Repository layout

```
.
├── scripts/
│   ├── paper_published/               ← variant calling + gene-level ASE (Sawarkar et al. 2026)
│   │   ├── 01_build_master_variants.py
│   │   ├── 02_count_ase_gene_level.sh
│   │   ├── 03_count_ase_per_variant.py
│   │   ├── lib/variant_analysis_library.py
│   │   └── README.md
│   └── followup_inprogress/           ← AlphaGenome regulatory variant prioritization (unpublished)
│       ├── 01_liftover_for_alphagenome.py
│       ├── 02_score_alphagenome.py
│       ├── 03_summarize_alphagenome.py
│       ├── 04_integrate_ase_alphagenome.py    # placeholder — rebuild needed
│       ├── 05_fig_alphagenome_results.py
│       ├── 06_fig_ucsc_tracks.py
│       ├── 07_fig_variant_annotation.py
│       ├── 08_fig_pipeline.py
│       └── README.md
│
├── data/
│   ├── small_inputs/
│   │   ├── paper_published/           ← ase_gene_level_counts.csv (Fig 3 data)
│   │   └── followup_inprogress/       ← AlphaGenome inputs, hg38 liftover, dipcall validation
│   ├── EXTERNAL_DATA.md                ← Zenodo / OSF / UCSC links for big files
│   └── README.md
│
├── results/
│   ├── paper_published/               ← master_variants.csv (4,515-variant catalog)
│   └── followup_inprogress/           ← AG outputs, integrated priority shortlist, FIMO TF disruption
│       └── README.md  (in results/)
│
├── figures/
│   ├── paper_published/               ← Fig3_ASE_combined (paper Fig 3)
│   └── followup_inprogress/scratch/   ← rough first-pass figures (not publication-ready)
│       └── README.md  (in figures/)
│
├── environment.yml                    ← conda env
├── CITATION.cff
└── LICENSE                            ← MIT
```

## Pipeline overview

### Paper-published (already in Sawarkar et al. 2026)

```
T2T H9 diploid genome (hap1 + hap2 FASTAs)
        │
        ├─ Winnowmap2 hap1↔hap2 alignment, CIGAR-parsed
        │       → het variant catalog (4,515 variants)
        │
        └─ haplotype-phased RNA-seq (.ase.bam, 4 lineages)
                → gene-level ASE counts ────────► Paper Fig 3
                → per-variant ASE counts (used as input to follow-up)
```

### Follow-up (in progress)

```
het variant catalog (T2T) ──► hg38 liftover (minimap2 + UCSC dipcall validation)
                                       │
                                       ▼
                             AlphaGenome API scoring (2,631 SNVs × 6 modalities)
                                       │
                              ┌────────┴─────────┐
                              ▼                  ▼
                  per-variant AG summary   brain-tissue subset
                              │                  │
                              └────────┬─────────┘
                                       ▼
                       integrate AG + ASE → priority_score
                                       │
                                       ▼
                       per-gene shortlist of candidate causal non-coding variants
```

**Headline follow-up result (preserved from a prior run):** SNCA shows astrocyte-specific hap2 expression bias (23/56 reads). The AlphaGenome-prioritized non-coding candidate explaining this bias is **`chr4:89,836,694 C>T`** (SNCA gene body, AlphaGenome astrocyte total-RNA-seq quantile = 0.9995, predicted up-regulation).

## Quick start

### 1. Clone and set up environment

```bash
git clone https://github.com/<username>/h9-diploid-ase-pipeline.git
cd h9-diploid-ase-pipeline
conda env create -f environment.yml
conda activate h9-ase
```

### 2. Inspect committed results

Final tables are already in `results/` — explore without rerunning anything:

```bash
head results/followup_inprogress/integrated_variant_analysis.csv
head results/paper_published/master_variants.csv
```

### 3. Re-run figures from committed results

```bash
python scripts/followup_inprogress/05_fig_alphagenome_results.py
python scripts/followup_inprogress/08_fig_pipeline.py
```

### 4. Re-run the full pipeline from raw data

You will need (see [`data/EXTERNAL_DATA.md`](data/EXTERNAL_DATA.md) for sources):

- **T2T H9 diploid FASTAs** (`t2t_h9_v01_hap{1,2}.fa`) — UCSC
- **Haplotype-phased RNA-seq BAMs** (`.ase.bam`, 8 files) — UCSC Genome Server Group
- **CAT GTF annotations** (`H9_HAP{1,2}.gtf`) — UCSC
- **UCSC dipcall diploid VCF** (`hg38.t2t_h9_v01.dip.vcf.gz`) — UCSC
- **AlphaGenome API key** — <https://deepmind.google/alphagenome>

Run the paper-published scripts first (in `scripts/paper_published/`, numbered 01–03), then the follow-up scripts (in `scripts/followup_inprogress/`, numbered 01–08).

**Note:** several scripts currently use hard-coded paths (e.g. `/Users/.../reganalysis/`); a path-portability refactor is on the to-do list.

## Citation

If you use this pipeline, please cite the H9 T2T diploid reference paper (see [`CITATION.cff`](CITATION.cff)).

## License

[MIT](LICENSE)

## Contact

Issues and questions: please open a GitHub issue.
