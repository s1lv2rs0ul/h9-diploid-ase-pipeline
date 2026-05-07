# Figures

## Published / final figures (top level)

These figures are publication-ready and appear (or are intended to appear) in the associated paper.

| File | What it shows |
|------|---------------|
| `Fig3_ASE_combined.pdf` / `.png` | Gene-level allele-specific expression: hap1 vs hap2 read counts and hap1 fractions across SNCA, HTT, LRRK2, GBA1, SMN1, SMN2 in iPSC, FPP, astrocyte, microglia. **Published as Figure 3 in Patil et al. bioRxiv 2026.03.26.714432.** |

## Scratch figures (`scratch/`)

Quick first-pass figures generated as an exercise to see what the data looks like — **not finalized and not formal scientific claims**. Layout, labeling, and panel selection are rough; kept here for transparency about what has been computed.

| File | What it shows | Status |
|---|---|---|
| `scratch/Fig_AlphaGenome_results.pdf` / `.png` | 4-panel AlphaGenome overview: max-quantile scatter, TF-disruption heatmap, regulatory-mechanism breakdown, haplotype burden vs ASE | Scratch — panel D labeling incomplete |
| `scratch/Fig_UCSC_Tracks.pdf` / `.png` | UCSC-browser-style variant landscape (±500 kb) per gene, CADD-binned variants | Scratch — styling not polished |
| `scratch/Fig_Variant_Annotation.pdf` / `.png` | 3-panel: CADD distributions, consequence-class breakdown, gene maps with GWAS hits | Scratch — most polished of the three; intended panel D not yet implemented |

## Regenerating figures

Figure scripts live in [`../scripts/`](../scripts/):

```bash
python ../scripts/12_fig_alphagenome_results.py     # Fig_AlphaGenome_results
python ../scripts/13_fig_ucsc_tracks.py             # Fig_UCSC_Tracks
python ../scripts/14_fig_variant_annotation.py      # Fig_Variant_Annotation
python ../scripts/15_fig_pipeline.py                # pipeline diagram
```

Note: scripts currently reference hard-coded paths under `~/reganalysis/`; portability refactor is on the to-do list.
