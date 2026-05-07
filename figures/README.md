# Figures

Figures are partitioned by project:

```
figures/
├── paper_published/             ← figures already in Sawarkar et al. 2026
│   ├── Fig3_ASE_combined.pdf
│   └── Fig3_ASE_combined.png
└── followup_inprogress/         ← figures from unpublished AlphaGenome follow-up
    └── scratch/                 ← rough first-pass exercises (not publication-ready)
        ├── Fig_AlphaGenome_results.pdf / .png
        ├── Fig_UCSC_Tracks.pdf / .png
        └── Fig_Variant_Annotation.pdf / .png
```

## paper_published/

| File | What it shows |
|------|---------------|
| `Fig3_ASE_combined.pdf` / `.png` | Gene-level allele-specific expression: hap1 vs hap2 read counts and hap1 fractions across SNCA, HTT, LRRK2, GBA1, SMN1, SMN2 in iPSC, FPP, astrocyte, microglia. **Figure 3 in Sawarkar et al. bioRxiv 2026.03.26.714432.** |

## followup_inprogress/scratch/

Quick first-pass figures generated as an exercise to see what the data looks like — **not finalized and not formal scientific claims**. Layout, labeling, and panel selection are rough; kept here for transparency about what has been computed.

| File | What it shows | Status |
|------|---------------|--------|
| `scratch/Fig_AlphaGenome_results.pdf` / `.png` | 4-panel AlphaGenome overview: max-quantile scatter, TF-disruption heatmap, regulatory-mechanism breakdown, haplotype burden vs ASE | Scratch — panel D labeling incomplete |
| `scratch/Fig_UCSC_Tracks.pdf` / `.png` | UCSC-browser-style variant landscape (±500 kb) per gene, CADD-binned variants | Scratch — styling not polished |
| `scratch/Fig_Variant_Annotation.pdf` / `.png` | 3-panel: CADD distributions, consequence-class breakdown, gene maps with GWAS hits | Scratch — most polished of the three; intended panel D not yet implemented |

## Regenerating figures

Figure scripts live in [`../scripts/followup_inprogress/`](../scripts/followup_inprogress/):

```bash
python ../scripts/followup_inprogress/05_fig_alphagenome_results.py
python ../scripts/followup_inprogress/06_fig_ucsc_tracks.py
python ../scripts/followup_inprogress/07_fig_variant_annotation.py
python ../scripts/followup_inprogress/08_fig_pipeline.py
```

Note: scripts currently reference hard-coded paths under `~/reganalysis/`; portability refactor is on the to-do list.
