#!/usr/bin/env python3
"""
═══════════════════════════════════════════════════════════════════════════════
PLACEHOLDER — Integrate ASE counts with AlphaGenome scores
═══════════════════════════════════════════════════════════════════════════════

This script needs to be reconstructed. The output it produced
(`results/integrated_variant_analysis.csv`, dated 2026-03-18) is preserved in
this repo, but the build logic was lost.

What it should do (inferred from output columns):
  1. Load:
       - data/small_inputs/snvs_for_alphagenome.csv         (T2T↔hg38 bridge)
       - data/small_inputs/alphagenome_variant_summary.csv  (max_quantile per variant)
       - results/alphagenome_brain_scores.csv               (brain-tissue scores)
       - results/master_variants.csv                        (per-variant ASE counts)
  2. Join AlphaGenome (hg38 variant_id) with ASE (T2T variant_id) via the bridge.
  3. Compute per-variant features:
       - ase_imbalance_astro  = |hap1 - hap2| / (hap1 + hap2) for astrocyte reads
       - in_gene_body         (boolean from master_variants)
       - all_max_quantile     (from AlphaGenome summary)
  4. Compute priority_score (reverse-engineered formula, ~99% of rows match):
       priority = all_max_quantile * (0.7 if in_gene_body else 0.5)
                + ase_imbalance * 0.3
  5. Apply hap_alt sign correction so AG predicted direction is interpretable
     in haplotype space (hap1-up vs hap2-up).
  6. Sort by priority_score (descending) and write to
     results/integrated_variant_analysis.csv

To verify your reconstruction matches the existing output:
  python scripts/11_integrate_ase_alphagenome.py
  # then diff against the committed results/integrated_variant_analysis.csv
"""

import sys

if __name__ == '__main__':
    print(__doc__)
    print("\nThis script is a placeholder. See docstring above for what to rebuild.")
    sys.exit(0)
