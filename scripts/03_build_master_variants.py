#!/usr/bin/env python3
"""
Build the single-source-of-truth master variant file.

Inputs (TRUSTED ONLY):
  - results/SNCA_variants_aligned.csv          (832 variants, hap1↔hap2 minimap2 alignment)
  - results/haplotype_variants_aligned.csv     (2826 variants: HTT, LRRK2, GBA)
  - results/SMN_variants_aligned.csv           (857 variants: SMN1/SMN2)
  - annotations/H9_HAP1.gtf                   (CAT gene annotation, hap1)
  - annotations/H9_HAP2.gtf                   (CAT gene annotation, hap2)

Reference genome: t2t_h9_v01 diploid (hap1.fa, hap2.fa)
  — This IS the genome of the cells. No liftover needed.

Output:
  - clean_analysis/variants/master_variants.csv

Coordinate system: Native H9 T2T haplotype coordinates throughout.
"""

import os
import sys
import pandas as pd
import numpy as np

BASE = os.path.expanduser('~/reganalysis')
OUT = os.path.join(BASE, 'clean_analysis/variants/master_variants.csv')

# ── Trusted input files ──
VARIANT_FILES = {
    'SNCA': os.path.join(BASE, 'results/SNCA_variants_aligned.csv'),
    'HTT_LRRK2_GBA': os.path.join(BASE, 'results/haplotype_variants_aligned.csv'),
    'SMN': os.path.join(BASE, 'results/SMN_variants_aligned.csv'),
}

GTF_HAP1 = os.path.join(BASE, 'annotations/H9_HAP1.gtf')
GTF_HAP2 = os.path.join(BASE, 'annotations/H9_HAP2.gtf')

# ── Gene extraction windows (native haplotype FASTA coords) ──
# These define the chromosome and FASTA start offset for each gene's extraction region.
#
# CRITICAL: Offsets are from PAF query/target names — the ground truth for what
# region was actually extracted. DO NOT derive from GTF gene body coordinates.
#
# Source PAF files:
#   SNCA:  results/winnowmap2_all/SNCA/winnowmap.paf
#   HTT:   results/winnowmap2_all/HTT/winnowmap.paf
#   LRRK2: results/winnowmap2_all/LRRK2/winnowmap.paf
#   GBA:   results/winnowmap2_all/GBA/winnowmap.paf
#   SMN1:  results/SMN_hap1_L2_vs_hap2.paf  (L2 locus, NOT L1!)
#
# Verified 2026-03-13: 100/100 samtools faidx spot-checks passed.
# Previous wrong offsets moved to DEPRECATED_wrong_offsets/
GENE_CONFIG = {
    'SNCA':  {'hap1_chrom': 'chr4_hap1',  'hap1_fasta_start': 92310499,   # PAF: chr4_hap1:92310499-93310499
              'hap2_chrom': 'chr4_hap2',  'hap2_fasta_start': 90246908},  # PAF: chr4_hap2:90246908-91246908
    'HTT':   {'hap1_chrom': 'chr4_hap1',  'hap1_fasta_start': 2585235,    # PAF: chr4_hap1:2585235-3585235
              'hap2_chrom': 'chr4_hap2',  'hap2_fasta_start': 2568744},   # PAF: chr4_hap2:2568744-3568744
    'LRRK2': {'hap1_chrom': 'chr12_hap1', 'hap1_fasta_start': 39837429,   # PAF: chr12_hap1:39837429-40837429
              'hap2_chrom': 'chr12_hap2', 'hap2_fasta_start': 40532905},  # PAF: chr12_hap2:40532905-41532905
    'GBA':   {'hap1_chrom': 'chr1_hap1',  'hap1_fasta_start': 155881588,  # PAF: chr1_hap1:155881588-156881588
              'hap2_chrom': 'chr1_hap2',  'hap2_fasta_start': 148778140}, # PAF: chr1_hap2:148778140-149778140
    'SMN1':  {'hap1_chrom': 'chr5_hap1',  'hap1_fasta_start': 74701987,   # PAF: chr5_hap1:74701987-75701987 (L2!)
              'hap2_chrom': 'chr5_hap2',  'hap2_fasta_start': 70969144},  # PAF: chr5_hap2:70969144-71969144
}


def parse_gtf_gene_boundaries(gtf_path, hap_label):
    """Extract gene start/end from CAT GTF for our target genes.

    CAT GTF has no 'gene' feature type — only 'transcript' and 'exon'.
    We derive gene boundaries from the union of all transcript spans.
    Gene name is in the attributes: gene_name "SNCA";
    """
    target_genes = {'SNCA', 'HTT', 'LRRK2', 'GBA1', 'SMN1', 'SMN2',
                    'SNCA-AS1', 'HTT-AS', 'LRRK2-DT', 'SMN-AS1'}
    genes = {}

    with open(gtf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != 'transcript':
                continue

            attrs = parts[8]
            gene_name = None
            for attr in attrs.split(';'):
                attr = attr.strip()
                if attr.startswith('gene_name '):
                    gene_name = attr.split('"')[1]
                    break

            if gene_name not in target_genes:
                continue

            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]

            if gene_name not in genes:
                genes[gene_name] = {
                    'chrom': chrom, 'start': start, 'end': end,
                    'strand': strand, 'hap': hap_label,
                    'n_transcripts': 1,
                }
            else:
                genes[gene_name]['start'] = min(genes[gene_name]['start'], start)
                genes[gene_name]['end'] = max(genes[gene_name]['end'], end)
                genes[gene_name]['n_transcripts'] += 1

    return genes


def load_snca():
    """Load SNCA variants (no gene column in source)."""
    df = pd.read_csv(VARIANT_FILES['SNCA'])
    df['gene'] = 'SNCA'
    cfg = GENE_CONFIG['SNCA']
    df['hap1_chrom'] = cfg['hap1_chrom']
    df['hap2_chrom'] = cfg['hap2_chrom']
    df['hap1_abs_pos'] = cfg['hap1_fasta_start'] + df['hap1_pos']
    df['hap2_abs_pos'] = cfg['hap2_fasta_start'] + df['hap2_pos']
    return df


def load_multi_gene():
    """Load HTT, LRRK2, GBA variants."""
    df = pd.read_csv(VARIANT_FILES['HTT_LRRK2_GBA'])
    rows = []
    for _, r in df.iterrows():
        gene = r['gene']
        cfg = GENE_CONFIG[gene]
        rows.append({
            **r.to_dict(),
            'hap1_chrom': cfg['hap1_chrom'],
            'hap2_chrom': cfg['hap2_chrom'],
            'hap1_abs_pos': cfg['hap1_fasta_start'] + r['hap1_pos'],
            'hap2_abs_pos': cfg['hap2_fasta_start'] + r['hap2_pos'],
        })
    return pd.DataFrame(rows)


def load_smn():
    """Load SMN variants (no gene column — assigned as SMN1 region)."""
    df = pd.read_csv(VARIANT_FILES['SMN'])
    df['gene'] = 'SMN1'  # All variants in this file are from the SMN1/2 region
    cfg = GENE_CONFIG['SMN1']
    df['hap1_chrom'] = cfg['hap1_chrom']
    df['hap2_chrom'] = cfg['hap2_chrom']
    df['hap1_abs_pos'] = cfg['hap1_fasta_start'] + df['hap1_pos']
    df['hap2_abs_pos'] = cfg['hap2_fasta_start'] + df['hap2_pos']
    return df


def classify_gene_body(df, gtf_hap1, gtf_hap2):
    """Classify each variant as in_gene_body using GTF boundaries."""
    results = []
    for _, r in df.iterrows():
        gene = r['gene']

        # Map gene name to GTF gene name
        gtf_gene = 'GBA1' if gene == 'GBA' else gene

        in_body_hap1 = False
        in_body_hap2 = False

        if gtf_gene in gtf_hap1:
            g1 = gtf_hap1[gtf_gene]
            if r['hap1_chrom'] == g1['chrom']:
                in_body_hap1 = g1['start'] <= r['hap1_abs_pos'] <= g1['end']

        if gtf_gene in gtf_hap2:
            g2 = gtf_hap2[gtf_gene]
            if r['hap2_chrom'] == g2['chrom']:
                in_body_hap2 = g2['start'] <= r['hap2_abs_pos'] <= g2['end']

        # A variant is "in gene body" if it falls within GTF boundaries on EITHER haplotype
        results.append({
            'in_gene_body_hap1': in_body_hap1,
            'in_gene_body_hap2': in_body_hap2,
            'in_gene_body': in_body_hap1 or in_body_hap2,
        })

    return pd.DataFrame(results)


def main():
    print("Building master variant file from trusted sources")
    print("=" * 60)
    print(f"Reference: t2t_h9_v01 diploid (native haplotype coordinates)")
    print()

    # ── Load GTF gene boundaries ──
    print("Parsing GTF gene boundaries ...")
    gtf_hap1 = parse_gtf_gene_boundaries(GTF_HAP1, 'hap1')
    gtf_hap2 = parse_gtf_gene_boundaries(GTF_HAP2, 'hap2')

    print("  HAP1 genes found:", {k: f"{v['chrom']}:{v['start']}-{v['end']}" for k, v in gtf_hap1.items()})
    print("  HAP2 genes found:", {k: f"{v['chrom']}:{v['start']}-{v['end']}" for k, v in gtf_hap2.items()})
    print()

    # ── Load variants ──
    print("Loading SNCA ...")
    df_snca = load_snca()
    print(f"  {len(df_snca)} variants")

    print("Loading HTT, LRRK2, GBA ...")
    df_multi = load_multi_gene()
    print(f"  {len(df_multi)} variants")

    print("Loading SMN1/2 ...")
    df_smn = load_smn()
    print(f"  {len(df_smn)} variants")

    # ── Merge ──
    # Standardize columns
    common_cols = ['gene', 'hap1_chrom', 'hap1_abs_pos', 'hap2_chrom', 'hap2_abs_pos',
                   'hap1_rel', 'hap1_base', 'hap2_base', 'type']

    for c in common_cols:
        for df in [df_snca, df_multi, df_smn]:
            if c not in df.columns:
                df[c] = np.nan

    df_all = pd.concat([
        df_snca[common_cols],
        df_multi[common_cols],
        df_smn[common_cols],
    ], ignore_index=True)

    # Rename type → variant_type for clarity
    df_all = df_all.rename(columns={'type': 'variant_type'})

    print(f"\nMerged: {len(df_all)} total variants")
    print(f"  Per gene: {df_all['gene'].value_counts().to_dict()}")
    print(f"  Variant types: SNV={( df_all['variant_type']=='SNV').sum()}, "
          f"Indel={( df_all['variant_type']!='SNV').sum()}")

    # ── Classify gene body using GTF ──
    print("\nClassifying gene body from CAT GTF ...")
    gb = classify_gene_body(df_all, gtf_hap1, gtf_hap2)
    df_all = pd.concat([df_all, gb], axis=1)

    for gene in ['SNCA', 'HTT', 'LRRK2', 'GBA', 'SMN1']:
        g = df_all[df_all['gene'] == gene]
        n_body = g['in_gene_body'].sum()
        n_h1 = g['in_gene_body_hap1'].sum()
        n_h2 = g['in_gene_body_hap2'].sum()
        print(f"  {gene:8s}  total={len(g):4d}  in_body={n_body:4d}  (hap1={n_h1}, hap2={n_h2})")

    # ── Save ──
    df_all.to_csv(OUT, index=False)
    print(f"\nSaved: {OUT}")
    print(f"  {len(df_all)} variants, {len(df_all.columns)} columns")
    print(f"  Columns: {list(df_all.columns)}")
    print(f"\nDone. All coordinates are native H9 T2T (t2t_h9_v01).")


if __name__ == '__main__':
    main()
