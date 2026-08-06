# Notebooks

Interactive Jupyter/Colab notebooks that let a biologist run the ASE analysis without local setup or coding.

## Available notebooks

### `hg38_to_T2T_ASE.ipynb` — Main analysis notebook

Open in Google Colab (no local install required):

<a href="https://colab.research.google.com/github/s1lv2rs0ul/h9-diploid-ase-pipeline/blob/main/notebooks/hg38_to_T2T_ASE.ipynb" target="_blank"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

**What it does:** takes any RNA-seq BAM (bulk or 10x scRNA-seq), re-aligns just the reads at the 5 focus genes to the H9 T2T diploid reference, and produces gene-level allele-specific expression counts + verification views.

**Who it's for:** biologists who understand the biology but don't want to write pipeline code. Every parameter is a form field; every step has a plain-English "why" explanation; every result is verifiable in UCSC Genome Browser.

**Sections:**
- A. Why we need to re-align (biology intro)
- B. Setup (one cell — installs tools, downloads T2T reference)
- C. Point at your data (form widgets + fuzzy cell-type matching)
- D. Auto-detect assay type (bulk vs 10x — adapts pipeline accordingly)
- E. Load and view the region (UCSC-style track — Phase 2)
- F. Extract → re-align → count ASE (Phase 2)
- G1. Base-pair verification — the sanity check
- G3. Open in UCSC Genome Browser
- G4. Direct links to raw reference FASTAs
- H. Troubleshooting + FAQ

**Status:** Phase 1 complete (setup, config, verification cells). Phase 2 (alignment pipeline) coming after STAR index is validated on the SSD.

**Runtime on Colab free tier:** ~15 min for one BAM at all 5 loci.

## Building / editing

The notebook is generated from a Python script for maintainability. To regenerate after edits:

```bash
python3 scripts/followup_inprogress/build_notebook.py
```
