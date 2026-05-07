# Figures

## Published / final figures (top level)

These figures are publication-ready and appear (or are intended to appear) in the associated paper.

| File | What it shows |
|------|---------------|
| `Fig3_ASE_combined.pdf` / `.png` | Gene-level allele-specific expression: hap1 vs hap2 read counts and hap1 fractions across SNCA, HTT, LRRK2, GBA1, SMN1, SMN2 in iPSC, FPP, astrocyte, microglia. **Published as Figure 3 in Patil et al. bioRxiv 2026.03.26.714432.** |

## Exploratory figures (`exploratory/`)

These figures were generated during downstream analysis exploration and are **not finalized**. Layout, labeling, and panel selection are draft-quality and may change. They are kept here for transparency about what has been computed, not as final scientific claims.

| File | What it shows | Status |
|---|---|---|
| `exploratory/Fig_AlphaGenome_results.pdf` / `.png` | 4-panel AlphaGenome overview: max-quantile scatter, TF-disruption heatmap, regulatory-mechanism breakdown, haplotype burden vs ASE | Exploratory — panel D labeling incomplete |
| `exploratory/Fig_UCSC_Tracks.pdf` / `.png` | UCSC-browser-style variant landscape (±500 kb) per gene, CADD-binned variants | Exploratory — styling not polished |
| `exploratory/Fig_Variant_Annotation.pdf` / `.png` | 3-panel: CADD distributions, consequence-class breakdown, gene maps with GWAS hits | Exploratory — most polished of the three; intended panel D not yet implemented |

## Regenerating figures

Figure scripts live in [`../scripts/`](../scripts/):

```bash
python ../scripts/12_fig_alphagenome_results.py     # Fig_AlphaGenome_results
python ../scripts/13_fig_ucsc_tracks.py             # Fig_UCSC_Tracks
python ../scripts/14_fig_variant_annotation.py      # Fig_Variant_Annotation
python ../scripts/15_fig_pipeline.py                # pipeline diagram
```

Note: scripts currently reference hard-coded paths under `~/reganalysis/`; portability refactor is on the to-do list.
