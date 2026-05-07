#!/usr/bin/env python3
"""
Lift trusted native-coord SNVs to hg38 for AlphaGenome scoring.

Strategy (UCSC-first approach):
  1. Load our corrected master_variants.csv (native T2T coordinates)
  2. Load UCSC diploid VCF het SNVs (authoritative hg38 coordinates)
  3. For each gene, align T2T hap1 region → hg38 via minimap2 to build
     a position mapping (CIGAR-based)
  4. Map each of our SNVs to hg38 via the minimap2 alignment
  5. Cross-validate: check if our mapped hg38 position matches a UCSC
     diploid VCF het SNV at the same position with the same alleles
  6. Tag each variant with match status (ucsc_match / ucsc_close / our_only)

This gives us:
  - UCSC-validated hg38 coordinates for ~90%+ of variants
  - Minimap2-derived hg38 coordinates for the remainder (our-only variants)
  - Confidence flag for downstream filtering

Cases for each het SNV (hap1_base vs hap2_base vs hg38_ref):
  Case 1: hap1_base == hg38_ref → hap_alt = "hap2"
  Case 2: hap2_base == hg38_ref → hap_alt = "hap1"
  Case 3: Neither matches hg38 → hap_alt = "both_differ" (tri-allelic)
    For AlphaGenome: need TWO scores (hap1 vs ref AND hap2 vs ref)

Inputs:
  - clean_analysis/variants/master_variants.csv  (4,515 variants, corrected coords)
  - clean_analysis/variants/het_snvs_from_ucsc_diploid_vcf.csv (3,272 UCSC het SNVs)
  - genomes/t2t_h9_v01_hap{1,2}.fa
  - genomes/hg38.fa

Output:
  - clean_analysis/variants/snvs_for_alphagenome.csv

Requires: minimap2, samtools, pandas
"""

import subprocess
import os
import re
import pandas as pd
import numpy as np

BASE = os.path.expanduser('~/reganalysis')
MASTER = os.path.join(BASE, 'clean_analysis/variants/master_variants.csv')
UCSC_HETS = os.path.join(BASE, 'clean_analysis/variants/het_snvs_from_ucsc_diploid_vcf.csv')
HAP1_FA = os.path.join(BASE, 'genomes/t2t_h9_v01_hap1.fa')
HAP2_FA = os.path.join(BASE, 'genomes/t2t_h9_v01_hap2.fa')
HG38_FA = os.path.join(BASE, 'genomes/hg38.fa')
OUT = os.path.join(BASE, 'clean_analysis/variants/snvs_for_alphagenome.csv')

# CORRECT offsets from PAF query names (verified against FASTA)
GENE_REGIONS_HAP1 = {
    'SNCA':  ('chr4_hap1',  92310499, 93310499),
    'HTT':   ('chr4_hap1',  2585235,  3585235),
    'LRRK2': ('chr12_hap1', 39837429, 40837429),
    'GBA':   ('chr1_hap1',  155881588, 156881588),
    'SMN1':  ('chr5_hap1',  74701987, 75701987),
}
GENE_REGIONS_HAP2 = {
    'SNCA':  ('chr4_hap2',  90246908, 91246908),
    'HTT':   ('chr4_hap2',  2568744,  3568744),
    'LRRK2': ('chr12_hap2', 40532905, 41532905),
    'GBA':   ('chr1_hap2',  148778140, 149778140),
    'SMN1':  ('chr5_hap2',  70969144, 71969144),
}


def align_region_to_hg38(chrom, start, end, hap_fa, label):
    """Align a T2T region to hg38 via minimap2. Returns CIGAR alignment blocks."""
    region = f"{chrom}:{start}-{end}"
    base_chrom = chrom.split('_')[0]

    # Extract hg38 chromosome for target
    hg38_chr_fa = f"/tmp/hg38_{base_chrom}.fa"
    if not os.path.exists(hg38_chr_fa):
        os.system(f"samtools faidx {HG38_FA} {base_chrom} > {hg38_chr_fa}")

    cmd = f"samtools faidx {hap_fa} {region} | minimap2 -a --cs {hg38_chr_fa} - 2>/dev/null"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    for line in result.stdout.split('\n'):
        if line.startswith('@') or not line.strip():
            continue
        parts = line.split('\t')
        if int(parts[1]) & 4:
            continue

        hg38_chrom = parts[2]
        hg38_start = int(parts[3])
        cigar = parts[5]
        mapq = int(parts[4])

        ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
        blocks = []
        q_pos = 0
        r_pos = 0
        for length_str, op in ops:
            length = int(length_str)
            if op in ('M', '=', 'X'):
                blocks.append((start + q_pos, hg38_start + r_pos, length))
                q_pos += length
                r_pos += length
            elif op == 'I':
                q_pos += length
            elif op == 'D':
                r_pos += length
            elif op == 'S':
                q_pos += length

        print(f"  {label}: {region} → {hg38_chrom}:{hg38_start} "
              f"(MAPQ={mapq}, {len(blocks)} blocks)")
        return blocks, hg38_chrom

    print(f"  {label}: {region} → UNMAPPED")
    return [], None


def query_to_ref(blocks, query_pos):
    """Convert T2T position to hg38 position using alignment blocks."""
    for q_start, r_start, block_len in blocks:
        if q_start <= query_pos < q_start + block_len:
            return r_start + (query_pos - q_start)
    return None


def batch_get_ref_bases(chrom, positions, hg38_fa):
    """Get hg38 reference bases for a batch of positions."""
    if not positions:
        return {}
    with open('/tmp/regions.txt', 'w') as f:
        for pos in sorted(positions):
            f.write(f"{chrom}:{pos}-{pos}\n")

    cmd = f"samtools faidx {hg38_fa} -r /tmp/regions.txt"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    ref_bases = {}
    current_pos = None
    for line in result.stdout.strip().split('\n'):
        if line.startswith('>'):
            coords = line[1:].split(':')[1]
            current_pos = int(coords.split('-')[0])
        elif current_pos is not None:
            ref_bases[current_pos] = line.strip().upper()
    return ref_bases


def main():
    print("=" * 70)
    print("Liftover: Native T2T → hg38 (UCSC-validated)")
    print("=" * 70)

    # ── Load data ──
    df = pd.read_csv(MASTER)
    snvs = df[df['variant_type'] == 'SNV'].copy()
    print(f"Master variants: {len(df)} total, {len(snvs)} SNVs")

    ucsc = pd.read_csv(UCSC_HETS)
    print(f"UCSC diploid VCF het SNVs: {len(ucsc)}")
    print()

    # Build UCSC lookup: (gene, hg38_pos) → row
    ucsc_lookup = {}
    for _, r in ucsc.iterrows():
        key = (r['gene'], int(r['hg38_pos']))
        ucsc_lookup[key] = r

    # ── Align gene regions to hg38 ──
    print("Aligning T2T regions to hg38 via minimap2:")
    gene_blocks = {}  # gene → (blocks, hg38_chrom)
    genes = ['SNCA', 'HTT', 'LRRK2', 'GBA', 'SMN1']

    for gene in genes:
        chrom1, s1, e1 = GENE_REGIONS_HAP1[gene]
        blocks1, hg38_chrom1 = align_region_to_hg38(chrom1, s1, e1, HAP1_FA, f"{gene} hap1")

        chrom2, s2, e2 = GENE_REGIONS_HAP2[gene]
        blocks2, hg38_chrom2 = align_region_to_hg38(chrom2, s2, e2, HAP2_FA, f"{gene} hap2")

        gene_blocks[gene] = {
            'hap1': (blocks1, hg38_chrom1),
            'hap2': (blocks2, hg38_chrom2),
        }
    print()

    # ── Map each SNV to hg38 ──
    print("Mapping SNV positions to hg38...")
    mapped_rows = []
    unmapped = 0

    for _, row in snvs.iterrows():
        gene = row['gene']
        hap1_pos = int(row['hap1_abs_pos'])
        hap2_pos = int(row['hap2_abs_pos'])

        blocks1, chrom1 = gene_blocks[gene]['hap1']
        blocks2, chrom2 = gene_blocks[gene]['hap2']

        # Try hap1 first (primary)
        hg38_pos = query_to_ref(blocks1, hap1_pos)
        source_hap = 'hap1'
        hg38_chrom = chrom1
        if hg38_pos is None:
            # Fallback to hap2
            hg38_pos = query_to_ref(blocks2, hap2_pos)
            source_hap = 'hap2'
            hg38_chrom = chrom2
        if hg38_pos is None:
            unmapped += 1
            continue

        mapped_rows.append({
            'gene': gene,
            'hg38_chrom': hg38_chrom,
            'hg38_pos': int(hg38_pos),
            'hap1_abs_pos': hap1_pos,
            'hap2_abs_pos': hap2_pos,
            'hap1_base': str(row['hap1_base']).upper(),
            'hap2_base': str(row['hap2_base']).upper(),
            'in_gene_body': row['in_gene_body'],
            'source_hap': source_hap,
        })

    print(f"  Mapped: {len(mapped_rows)}, Unmapped: {unmapped}")

    # ── Look up hg38 REF bases ──
    print("\nLooking up hg38 reference bases...")
    df_mapped = pd.DataFrame(mapped_rows)
    positions_by_chrom = {}
    for _, r in df_mapped.iterrows():
        chrom = r['hg38_chrom']
        if chrom not in positions_by_chrom:
            positions_by_chrom[chrom] = set()
        positions_by_chrom[chrom].add(int(r['hg38_pos']))

    all_ref_bases = {}
    for chrom, positions in positions_by_chrom.items():
        refs = batch_get_ref_bases(chrom, positions, HG38_FA)
        all_ref_bases.update({(chrom, pos): base for pos, base in refs.items()})
    print(f"  Retrieved {len(all_ref_bases)} REF bases")

    # ── Assign REF/ALT, hap_alt, and cross-validate against UCSC ──
    print("\nDetermining REF/ALT/hap_alt and cross-validating against UCSC diploid VCF...")
    final_rows = []
    ucsc_exact = 0
    ucsc_pos_only = 0
    our_only = 0

    for _, row in df_mapped.iterrows():
        ref = all_ref_bases.get((row['hg38_chrom'], int(row['hg38_pos'])))
        h1 = row['hap1_base']
        h2 = row['hap2_base']

        if ref is None:
            continue

        # Determine hap_alt
        if ref == h1:
            hg38_alt = h2
            hap_alt = 'hap2'
        elif ref == h2:
            hg38_alt = h1
            hap_alt = 'hap1'
        else:
            hg38_alt = h2  # Convention: use hap2 as primary ALT
            hap_alt = 'both_differ'

        # Cross-validate against UCSC diploid VCF
        ucsc_key = (row['gene'], int(row['hg38_pos']))
        ucsc_row = ucsc_lookup.get(ucsc_key)
        if ucsc_row is not None:
            # Position matches — check alleles
            if (str(ucsc_row['hg38_ref']) == ref and
                str(ucsc_row['hg38_alt']) in (h1, h2)):
                match_status = 'ucsc_exact'
                ucsc_exact += 1
            else:
                match_status = 'ucsc_pos_match'
                ucsc_pos_only += 1
        else:
            match_status = 'our_only'
            our_only += 1

        variant_id = f"{row['gene']}_{row['hg38_chrom']}_{int(row['hg38_pos'])}_{ref}_{hg38_alt}"

        final_rows.append({
            'gene': row['gene'],
            'hg38_chrom': row['hg38_chrom'],
            'hg38_pos': int(row['hg38_pos']),
            'hg38_ref': ref,
            'hg38_alt': hg38_alt,
            'hap_alt': hap_alt,
            'hap1_base': h1,
            'hap2_base': h2,
            'hap1_abs_pos': row['hap1_abs_pos'],
            'hap2_abs_pos': row['hap2_abs_pos'],
            'in_gene_body': row['in_gene_body'],
            'variant_id': variant_id,
            'ucsc_validated': match_status,
        })

    df_final = pd.DataFrame(final_rows)

    # ── Report ──
    print(f"\n{'=' * 70}")
    print("Summary:")
    print(f"  Total SNVs in master:   {len(snvs)}")
    print(f"  Mapped to hg38:         {len(df_final)}")
    print(f"  Unmapped (in indels):    {len(snvs) - len(df_final)}")

    print(f"\nUCSC validation:")
    print(f"  Exact match (pos+alleles): {ucsc_exact} ({100*ucsc_exact/len(df_final):.1f}%)")
    print(f"  Position-only match:       {ucsc_pos_only} ({100*ucsc_pos_only/len(df_final):.1f}%)")
    print(f"  Our-only (no UCSC match):  {our_only} ({100*our_only/len(df_final):.1f}%)")

    print(f"\nhap_alt distribution:")
    for k, v in df_final['hap_alt'].value_counts().items():
        pct = 100 * v / len(df_final)
        label = {
            'hap1': 'hap1 differs from hg38 (hap2=REF)',
            'hap2': 'hap2 differs from hg38 (hap1=REF)',
            'both_differ': 'BOTH differ from hg38 (tri-allelic)',
        }.get(k, k)
        print(f"  {k:15s}  {v:5d} ({pct:5.1f}%)  — {label}")

    print(f"\nPer gene:")
    for gene in genes:
        g = df_final[df_final['gene'] == gene]
        h1 = (g['hap_alt'] == 'hap1').sum()
        h2 = (g['hap_alt'] == 'hap2').sum()
        both = (g['hap_alt'] == 'both_differ').sum()
        ucsc_ok = (g['ucsc_validated'] == 'ucsc_exact').sum()
        print(f"  {gene:8s}  total={len(g):4d}  hap1={h1:4d}  hap2={h2:4d}  "
              f"both_differ={both:3d}  ucsc_validated={ucsc_ok:4d} ({100*ucsc_ok/max(len(g),1):.0f}%)")

    # ── Save ──
    df_final.to_csv(OUT, index=False)
    print(f"\nSaved: {OUT}")
    print(f"  {len(df_final)} variants, {len(df_final.columns)} columns")

    # API call estimate
    simple = (df_final['hap_alt'] != 'both_differ').sum()
    tri = (df_final['hap_alt'] == 'both_differ').sum()
    print(f"\nAlphaGenome API calls needed:")
    print(f"  {simple} simple + {tri}×2 tri-allelic = {simple + tri*2} total")


if __name__ == '__main__':
    main()
