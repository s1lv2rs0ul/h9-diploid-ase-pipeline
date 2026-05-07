#!/usr/bin/env python3
"""
═══════════════════════════════════════════════════════════════════════════════
AlphaGenome Variant Effect Scoring — H9 Het SNVs (hg38)
═══════════════════════════════════════════════════════════════════════════════

Scores SNVs from het_variants_hg38_official.csv using the AlphaGenome API.
Predicts variant effects on gene expression, chromatin accessibility,
histone modifications, and transcription factor binding.

Input:  vep_redo/het_variants_hg38_official.csv (hg38 coordinates)
Output: variant_catalog/alphagenome_scores.csv

Setup:
  pip install alphagenome
  export ALPHAGENOME_API_KEY="your-key-here"

Usage:
  python3 score_alphagenome.py
  python3 score_alphagenome.py --gene SNCA          # score one gene
  python3 score_alphagenome.py --gene SNCA --dry-run # preview without API calls
"""

import os
import sys
import csv
import argparse
import time

# ════════════════════════════════════════════════════════════════
# PATHS
# ════════════════════════════════════════════════════════════════
BASE     = os.path.expanduser('~/reganalysis')
HET_CSV  = os.path.join(BASE, 'vep_redo/het_variants_hg38_official.csv')
OUT_DIR  = os.path.join(BASE, 'variant_catalog')
OUT_CSV  = os.path.join(OUT_DIR, 'alphagenome_scores.csv')

GENES = ['SNCA', 'HTT', 'LRRK2', 'GBA', 'SMN1']


def load_snvs(gene_filter=None):
    """Load SNVs from het_variants_hg38_official.csv."""
    variants = []
    with open(HET_CSV) as f:
        reader = csv.DictReader(f)
        for row in reader:
            gene = row['gene']
            if gene not in GENES:
                continue
            if gene_filter and gene != gene_filter:
                continue
            ref = row['ref'].upper()
            alt = row['alt'].upper()
            if len(ref) != 1 or len(alt) != 1 or ',' in alt:
                continue
            variants.append({
                'gene':      gene,
                'chrom':     row['chrom'],
                'pos':       int(row['pos']),
                'ref':       ref,
                'alt':       alt,
                'hap_alt':   row.get('hap_alt', ''),
                'in_gene_body': row.get('in_gene_body', ''),
                'variant_id': f"{row['chrom']}_{row['pos']}_{ref}_{alt}",
            })
    return variants


def load_already_scored(out_csv):
    """Load variant IDs already in the output CSV to enable resuming."""
    scored = set()
    if os.path.exists(out_csv):
        with open(out_csv) as f:
            reader = csv.DictReader(f)
            for row in reader:
                scored.add(row.get('variant_id', ''))
    return scored


def score_variants(variants, api_key, out_csv, seq_length_key='SEQUENCE_LENGTH_1MB'):
    """
    Score variants using AlphaGenome API.
    Appends tidy scores to out_csv incrementally (resumable).
    """
    from alphagenome.data import genome
    from alphagenome.models import dna_client, variant_scorers

    model = dna_client.create(api_key)

    # Select scorers
    all_scorers = variant_scorers.RECOMMENDED_VARIANT_SCORERS
    selected_keys = [
        'rna_seq', 'cage', 'atac', 'dnase',
        'chip_histone', 'chip_tf',
    ]
    selected_scorers = [
        all_scorers[key] for key in all_scorers
        if key.lower() in selected_keys
    ]

    sequence_length = dna_client.SUPPORTED_SEQUENCE_LENGTHS[seq_length_key]
    organism = dna_client.Organism.HOMO_SAPIENS

    # Remove unsupported scorers
    unsupported = [
        s for s in selected_scorers
        if organism.value not in
           variant_scorers.SUPPORTED_ORGANISMS[s.base_variant_scorer]
    ]
    for s in unsupported:
        selected_scorers.remove(s)

    # Check which variants are already scored (using our_variant_id column)
    already_scored = set()
    if os.path.exists(out_csv):
        import pandas as pd
        try:
            existing = pd.read_csv(out_csv)
            if 'our_variant_id' in existing.columns:
                already_scored = set(existing['our_variant_id'].unique())
        except Exception:
            pass
    print(f"  Already scored: {len(already_scored)} variants")

    # Determine if we need to write header
    write_header = not os.path.exists(out_csv) or os.path.getsize(out_csv) == 0

    total = len(variants)
    scored_count = 0
    skipped_count = 0

    for i, v in enumerate(variants):
        variant_id = v['variant_id']

        # Skip if already scored
        if variant_id in already_scored:
            skipped_count += 1
            if skipped_count % 100 == 0:
                print(f"  [{i+1}/{total}] Skipped {skipped_count} already scored ...")
            continue

        if (scored_count + 1) % 10 == 0 or scored_count == 0:
            print(f"  [{i+1}/{total}] Scoring {v['gene']} "
                  f"{v['chrom']}:{v['pos']} {v['ref']}>{v['alt']} ...")

        ag_variant = genome.Variant(
            chromosome=v['chrom'],
            position=v['pos'],
            reference_bases=v['ref'],
            alternate_bases=v['alt'],
            name=variant_id,
        )
        interval = ag_variant.reference_interval.resize(sequence_length)

        try:
            scores = model.score_variant(
                interval=interval,
                variant=ag_variant,
                variant_scorers=selected_scorers,
                organism=organism,
            )

            # score_variant returns a list of AnnData objects
            df = variant_scorers.tidy_scores(scores)
            # Map our gene/body info using variant_id column
            # variant_id format from API: "chr4:89312890:G>A"
            df['gene'] = v['gene']
            df['in_gene_body'] = v['in_gene_body']
            df['our_variant_id'] = variant_id

            # Append to CSV
            df.to_csv(out_csv, mode='a', header=write_header, index=False)
            write_header = False  # only write header once
            scored_count += 1

        except Exception as e:
            print(f"    ERROR scoring {variant_id}: {e}")
            time.sleep(2)

    return scored_count, skipped_count


def main():
    parser = argparse.ArgumentParser(
        description='Score variants with AlphaGenome API')
    parser.add_argument('--gene', type=str, default=None,
                        help='Score only one gene (SNCA, HTT, LRRK2, GBA, SMN1)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Show variants without making API calls')
    parser.add_argument('--seq-length', type=str, default='1MB',
                        choices=['16KB', '100KB', '500KB', '1MB'],
                        help='Sequence context length (default: 1MB)')
    args = parser.parse_args()

    # Load variants
    variants = load_snvs(gene_filter=args.gene)
    print(f"Loaded {len(variants)} SNVs from {HET_CSV}")
    for g in GENES:
        count = sum(1 for v in variants if v['gene'] == g)
        if count > 0:
            print(f"  {g}: {count} SNVs")

    if args.dry_run:
        print("\n[DRY RUN] Would score the above variants. Exiting.")
        for v in variants[:5]:
            print(f"  {v['variant_id']}  gene={v['gene']}  "
                  f"body={v['in_gene_body']}")
        return

    # Get API key
    api_key = os.environ.get('ALPHAGENOME_API_KEY')
    if not api_key:
        print("ERROR: Set ALPHAGENOME_API_KEY environment variable.")
        print("  export ALPHAGENOME_API_KEY='your-key-here'")
        sys.exit(1)

    # Score (appends incrementally to CSV)
    seq_key = f'SEQUENCE_LENGTH_{args.seq_length}'
    print(f"\nScoring with AlphaGenome (seq_length={args.seq_length}) ...")
    print(f"Output: {OUT_CSV}")
    scored, skipped = score_variants(variants, api_key, OUT_CSV,
                                     seq_length_key=seq_key)

    print(f"\nDone! Scored {scored} new variants, skipped {skipped} already scored.")
    if os.path.exists(OUT_CSV):
        with open(OUT_CSV) as f:
            total_rows = sum(1 for _ in f) - 1
        print(f"Total rows in {OUT_CSV}: {total_rows}")


if __name__ == '__main__':
    main()
