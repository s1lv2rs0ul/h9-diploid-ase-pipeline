#!/usr/bin/env python3
"""
Summarize the 22GB AlphaGenome raw scores into manageable files:

1. Top variant effects (highest quantile scores per variant)
2. Brain/neural-specific effects (filtered by relevant ontology terms)
3. Per-gene summary statistics

Reads in chunks to avoid memory issues.

Usage:
  python3 summarize_alphagenome.py
"""

import os
import pandas as pd
import numpy as np

BASE    = os.path.expanduser('~/reganalysis/variant_catalog')
RAW_CSV = os.path.join(BASE, 'alphagenome_scores.csv')
CHUNK   = 500_000  # rows per chunk

# Brain/neural ontology terms (UBERON and CL)
BRAIN_ONTOLOGIES = {
    'UBERON:0000955',  # brain
    'UBERON:0001890',  # forebrain
    'UBERON:0001891',  # midbrain
    'UBERON:0002298',  # cerebral cortex
    'UBERON:0001882',  # nucleus accumbens
    'UBERON:0001873',  # caudate nucleus
    'UBERON:0001874',  # putamen
    'UBERON:0002037',  # cerebellum
    'UBERON:0002038',  # substantia nigra
    'UBERON:0001876',  # amygdala
    'UBERON:0002421',  # hippocampus
    'UBERON:0009834',  # dorsolateral prefrontal cortex
    'UBERON:0001870',  # frontal cortex
    'CL:0000540',      # neuron
    'CL:0000127',      # astrocyte
    'CL:0000129',      # microglial cell
    'CL:0000128',      # oligodendrocyte
    'CL:0002319',      # neural cell
    'CL:0000030',      # glioblast
    'CL:0002633',      # neural progenitor
    'CL:0000047',      # neuronal stem cell
    'EFO:0002067',     # neural cell (EFO)
    'EFO:0003042',     # brain
}

# Also match by biosample_name keywords
BRAIN_KEYWORDS = [
    'brain', 'neural', 'neuron', 'astrocyte', 'glia', 'cortex',
    'cerebellum', 'hippocampus', 'substantia nigra', 'midbrain',
    'forebrain', 'caudate', 'putamen', 'amygdala', 'frontal',
    'temporal', 'occipital', 'parietal', 'cerebral', 'dorsolateral',
    'neuroblastoma', 'iPSC', 'iPS',
]


def is_brain_related(row):
    """Check if a score row is brain/neural-related."""
    ont = str(row.get('ontology_curie', ''))
    bio = str(row.get('biosample_name', '')).lower()

    if ont in BRAIN_ONTOLOGIES:
        return True
    for kw in BRAIN_KEYWORDS:
        if kw.lower() in bio:
            return True
    return False


def main():
    print("Summarizing AlphaGenome scores (reading in chunks) ...")
    print(f"Input: {RAW_CSV}")
    print(f"Size: {os.path.getsize(RAW_CSV)/1e9:.1f} GB")

    # Accumulators
    top_per_variant = []    # top N scores per variant
    brain_scores = []       # brain-related scores
    gene_stats = []         # per-gene aggregates

    chunk_num = 0
    total_rows = 0

    for chunk in pd.read_csv(RAW_CSV, chunksize=CHUNK, low_memory=False):
        chunk_num += 1
        total_rows += len(chunk)
        if chunk_num % 10 == 1:
            print(f"  Chunk {chunk_num}: {total_rows:,} rows processed ...")

        # Convert scores to numeric
        for col in ['raw_score', 'quantile_score']:
            if col in chunk.columns:
                chunk[col] = pd.to_numeric(chunk[col], errors='coerce')

        # --- 1. Top scores per variant (quantile_score > 0.9) ---
        if 'quantile_score' in chunk.columns:
            high = chunk[chunk['quantile_score'] >= 0.9].copy()
            if len(high) > 0:
                top_per_variant.append(high)

        # --- 2. Brain-related scores ---
        brain_mask = chunk.apply(is_brain_related, axis=1)
        brain_chunk = chunk[brain_mask].copy()
        if len(brain_chunk) > 0:
            brain_scores.append(brain_chunk)

    print(f"\nTotal rows processed: {total_rows:,}")

    # --- Combine and save ---

    # 1. Top variant effects
    if top_per_variant:
        df_top = pd.concat(top_per_variant, ignore_index=True)
        # Keep top 5 per variant by quantile_score
        df_top = (df_top.sort_values('quantile_score', ascending=False)
                  .groupby('our_variant_id').head(5)
                  .reset_index(drop=True))
        out1 = os.path.join(BASE, 'alphagenome_top_effects.csv')
        df_top.to_csv(out1, index=False)
        print(f"\nTop effects (quantile >= 0.9): {len(df_top):,} rows")
        print(f"  Unique variants with high scores: {df_top['our_variant_id'].nunique()}")
        print(f"  Saved: {out1}")
        print(f"  Size: {os.path.getsize(out1)/1e6:.1f} MB")

        # Per-gene summary of top effects
        print("\n  Per-gene high-impact variants:")
        for gene in ['SNCA', 'HTT', 'LRRK2', 'GBA', 'SMN1']:
            g = df_top[df_top['gene'] == gene]
            print(f"    {gene}: {g['our_variant_id'].nunique()} variants "
                  f"with high scores")
    else:
        print("\nNo high-scoring variants found (quantile >= 0.9)")

    # 2. Brain-related scores
    if brain_scores:
        df_brain = pd.concat(brain_scores, ignore_index=True)
        # Keep top 10 per variant for brain
        if 'quantile_score' in df_brain.columns:
            df_brain = (df_brain.sort_values('quantile_score', ascending=False)
                       .groupby('our_variant_id').head(10)
                       .reset_index(drop=True))
        out2 = os.path.join(BASE, 'alphagenome_brain_scores.csv')
        df_brain.to_csv(out2, index=False)
        print(f"\nBrain/neural scores: {len(df_brain):,} rows")
        print(f"  Unique variants: {df_brain['our_variant_id'].nunique()}")
        print(f"  Saved: {out2}")
        print(f"  Size: {os.path.getsize(out2)/1e6:.1f} MB")

        # Top brain-specific variants per gene
        if 'quantile_score' in df_brain.columns:
            print("\n  Top brain-specific variants per gene:")
            for gene in ['SNCA', 'HTT', 'LRRK2', 'GBA', 'SMN1']:
                g = df_brain[df_brain['gene'] == gene]
                if len(g) == 0:
                    print(f"    {gene}: no brain scores")
                    continue
                top = g.nlargest(3, 'quantile_score')
                print(f"    {gene} (top 3):")
                for _, r in top.iterrows():
                    print(f"      {r['our_variant_id']}  "
                          f"q={r['quantile_score']:.3f}  "
                          f"{r.get('biosample_name','')}  "
                          f"{r.get('output_type','')}")
    else:
        print("\nNo brain-related scores found")

    # 3. Per-variant max scores (compact summary)
    print("\nGenerating per-variant max score summary ...")
    max_scores = []
    for chunk in pd.read_csv(RAW_CSV, chunksize=CHUNK, low_memory=False):
        for col in ['raw_score', 'quantile_score']:
            if col in chunk.columns:
                chunk[col] = pd.to_numeric(chunk[col], errors='coerce')

        agg = (chunk.groupby(['our_variant_id', 'gene', 'in_gene_body'])
               .agg(
                   max_quantile=('quantile_score', 'max'),
                   max_raw=('raw_score', lambda x: x.abs().max()),
                   n_high_scores=('quantile_score', lambda x: (x >= 0.9).sum()),
                   n_scores=('quantile_score', 'count'),
               ).reset_index())
        max_scores.append(agg)

    df_max = pd.concat(max_scores, ignore_index=True)
    # Re-aggregate across chunks
    df_summary = (df_max.groupby(['our_variant_id', 'gene', 'in_gene_body'])
                  .agg(
                      max_quantile=('max_quantile', 'max'),
                      max_raw=('max_raw', 'max'),
                      n_high_scores=('n_high_scores', 'sum'),
                      n_scores=('n_scores', 'sum'),
                  ).reset_index()
                  .sort_values('max_quantile', ascending=False))

    out3 = os.path.join(BASE, 'alphagenome_variant_summary.csv')
    df_summary.to_csv(out3, index=False)
    print(f"Per-variant summary: {len(df_summary)} variants")
    print(f"  Saved: {out3}")
    print(f"  Size: {os.path.getsize(out3)/1e6:.1f} MB")

    print("\nDone!")


if __name__ == '__main__':
    main()
