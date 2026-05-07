#!/usr/bin/env python3
"""
UCSC Genome Browser–style track figure — shows variant locations
along each gene's haplotype region with multiple annotation tracks.

Tracks per gene (top→bottom):
  1. Chromosome ideogram / coordinate ruler
  2. Gene body (exons as thick boxes, introns as thin line)
  3. Variant density histogram
  4. Individual variant positions coloured by CADD score
  5. AlphaGenome high-impact variants (quantile ≥ 0.95)
  6. Known GWAS / literature hits
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import LineCollection
import warnings
warnings.filterwarnings('ignore')

BASE = '/Users/shruti.patil/reganalysis'
VEP_CSV  = f'{BASE}/vep_redo/vep_het_official_hg38.csv'
AG_CSV   = f'{BASE}/variant_catalog/alphagenome/scores/alphagenome_variant_summary.csv'
OUT_FIG  = f'{BASE}/variant_catalog/alphagenome/figures/Fig_UCSC_Tracks'

GENE_COLORS = {
    'SNCA': '#E69F00', 'HTT': '#56B4E9', 'LRRK2': '#009E73',
    'GBA': '#F0E442', 'SMN1': '#CC79A7',
}
GENES = ['SNCA', 'HTT', 'LRRK2', 'GBA', 'SMN1']

GENE_BODIES_HG38 = {
    'SNCA':  (89724099, 89838315, 'chr4',  '+'),
    'HTT':   (3074681,  3243960,  'chr4',  '+'),
    'LRRK2': (40224977, 40369285, 'chr12', '+'),
    'GBA':   (155234452,155244699,'chr1',  '-'),
    'SMN1':  (70924941, 70953015, 'chr5',  '+'),
}

# Simplified exon structures (approx positions from UCSC hg38)
EXONS_HG38 = {
    'SNCA': [(89724099, 89724610), (89755175, 89755328), (89770634, 89770780),
             (89800400, 89800570), (89825360, 89838315)],
    'HTT':  [(3074681, 3080920), (3088030, 3088240), (3100000, 3100250),
             (3120000, 3120400), (3140000, 3140500), (3180000, 3180800),
             (3220000, 3220500), (3240000, 3243960)],
    'LRRK2': [(40224977, 40226300), (40240000, 40240500), (40260000, 40260600),
              (40280000, 40280800), (40300000, 40301200), (40320000, 40321000),
              (40340000, 40341500), (40360000, 40369285)],
    'GBA':  [(155234452, 155235200), (155236000, 155236500), (155238000, 155238800),
             (155240000, 155240600), (155242000, 155244699)],
    'SMN1': [(70924941, 70925800), (70930000, 70930400), (70935000, 70935500),
             (70940000, 70940800), (70945000, 70945500), (70950000, 70953015)],
}

KNOWN_GWAS = {
    'rs356219':   {'gene': 'SNCA',  'label': 'PD GWAS lead'},
    'rs356220':   {'gene': 'SNCA',  'label': 'PD risk'},
    'rs2583988':  {'gene': 'SNCA',  'label': 'ncRNA exon'},
    'rs35303786': {'gene': 'LRRK2', 'label': 'Missense'},
    'rs2638869':  {'gene': 'LRRK2', 'label': 'Stop gained'},
    'rs11564230': {'gene': 'LRRK2', 'label': 'Splice donor'},
    'rs11549516': {'gene': 'HTT',   'label': 'Missense'},
}


def load_data():
    vep = pd.read_csv(VEP_CSV)
    snvs = vep[(vep['ref'].str.len() == 1) & (vep['alt'].str.len() == 1) &
               (~vep['alt'].str.contains(',', na=False))].copy()
    ag = pd.read_csv(AG_CSV)
    return snvs, ag


def draw_gene_tracks(fig, gs_row, gene, snvs, ag, row_idx, n_genes):
    """Draw all tracks for one gene in a horizontal strip."""
    gb = GENE_BODIES_HG38[gene]
    gene_start, gene_end, chrom, strand = gb
    center = (gene_start + gene_end) / 2
    window = 500000  # ±500kb
    region_start = int(center - window)
    region_end = int(center + window)

    g_snvs = snvs[snvs['gene'] == gene].copy()
    g_ag = ag[ag['gene'] == gene].copy()

    # Create sub-gridspec: 4 track rows
    inner = gs_row.subgridspec(4, 1, height_ratios=[0.8, 1.2, 1.0, 0.6], hspace=0.05)

    # ── Track 1: Coordinate ruler + gene body ──
    ax1 = fig.add_subplot(inner[0])
    ax1.set_xlim(region_start, region_end)
    ax1.set_ylim(-0.5, 1.5)

    # Gene label
    ax1.text(region_start - (region_end - region_start) * 0.01, 0.5,
             gene, fontsize=12, fontweight='bold', ha='right', va='center',
             color=GENE_COLORS[gene],
             bbox=dict(boxstyle='round,pad=0.3', facecolor=GENE_COLORS[gene],
                       alpha=0.15, edgecolor=GENE_COLORS[gene]))

    # Intron line
    ax1.plot([gene_start, gene_end], [0.5, 0.5], color='#333333', linewidth=1.5, zorder=2)

    # Exons as thick rectangles
    for ex_start, ex_end in EXONS_HG38[gene]:
        ax1.add_patch(mpatches.FancyBboxPatch(
            (ex_start, 0.15), ex_end - ex_start, 0.7,
            boxstyle='round,pad=0', facecolor=GENE_COLORS[gene],
            edgecolor='#333333', linewidth=0.5, zorder=3))

    # Strand arrow
    arrow_x = gene_end + (region_end - region_start) * 0.01 if strand == '+' else gene_start - (region_end - region_start) * 0.01
    ax1.annotate('', xy=(arrow_x, 0.5),
                 xytext=(gene_end if strand == '+' else gene_start, 0.5),
                 arrowprops=dict(arrowstyle='->', color='#333333', lw=1.5))

    # Coordinate ticks
    tick_interval = 200000
    ticks = np.arange(region_start, region_end + 1, tick_interval)
    for t in ticks:
        ax1.plot([t, t], [1.1, 1.3], color='grey', linewidth=0.5)
    ax1.set_xticks(ticks)
    ax1.set_xticklabels([f'{t/1e6:.1f}Mb' for t in ticks], fontsize=6, color='grey')
    ax1.tick_params(axis='x', length=0, pad=1)
    ax1.xaxis.set_ticks_position('top')

    # Chromosome label
    ax1.text(region_start, 1.4, f'{chrom}', fontsize=7, color='grey', va='bottom')

    ax1.set_yticks([])
    ax1.spines['top'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # ── Track 2: Variant density histogram ──
    ax2 = fig.add_subplot(inner[1], sharex=ax1)
    if len(g_snvs) > 0:
        positions = g_snvs['pos'].values
        bins = np.linspace(region_start, region_end, 100)
        counts, edges = np.histogram(positions, bins=bins)
        bin_centers = (edges[:-1] + edges[1:]) / 2
        ax2.fill_between(bin_centers, counts, alpha=0.4, color=GENE_COLORS[gene],
                         step='mid', linewidth=0)
        ax2.step(bin_centers, counts, color=GENE_COLORS[gene], linewidth=0.8, where='mid')

    # Shade gene body
    ax2.axvspan(gene_start, gene_end, alpha=0.08, color=GENE_COLORS[gene], zorder=0)
    ax2.set_ylabel('Density', fontsize=7, labelpad=2)
    ax2.tick_params(axis='y', labelsize=6)
    ax2.set_ylim(bottom=0)
    plt.setp(ax2.get_xticklabels(), visible=False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    # ── Track 3: Individual variants by CADD + AlphaGenome ──
    ax3 = fig.add_subplot(inner[2], sharex=ax1)
    if len(g_snvs) > 0:
        # Bottom row: all variants colored by CADD
        cadd = g_snvs['cadd_phred'].fillna(0).values
        positions = g_snvs['pos'].values

        # Color map: grey → orange → red
        colors = []
        for c in cadd:
            if c >= 20:
                colors.append('#d62728')
            elif c >= 15:
                colors.append('#ff7f0e')
            elif c >= 10:
                colors.append('#ffd92f')
            else:
                colors.append('#aaaaaa')

        ax3.scatter(positions, np.zeros(len(positions)) + 0.2, c=colors,
                    s=4, alpha=0.6, zorder=2, edgecolors='none')

        # Top row: AlphaGenome high-impact (quantile ≥ 0.95)
        if len(g_ag) > 0:
            g_ag['pos'] = g_ag['our_variant_id'].str.split('_').str[1].astype(int)
            high = g_ag[g_ag['max_quantile'] >= 0.95]
            if len(high) > 0:
                ax3.scatter(high['pos'].values, np.ones(len(high)) * 0.7,
                            c='#d62728', s=12, alpha=0.8, marker='D', zorder=3,
                            edgecolors='black', linewidth=0.3)

    # Known GWAS variants
    gene_gwas = {k: v for k, v in KNOWN_GWAS.items() if v['gene'] == gene}
    for rsid, info in gene_gwas.items():
        match = g_snvs[g_snvs['rsid'] == rsid]
        if not match.empty:
            pos = match.iloc[0]['pos']
            ax3.scatter([pos], [0.45], s=80, c='red', marker='*', zorder=5,
                        edgecolors='black', linewidth=0.5)
            ax3.annotate(rsid, xy=(pos, 0.55), fontsize=5.5, ha='center',
                         va='bottom', color='red', fontweight='bold',
                         rotation=30)

    ax3.axvspan(gene_start, gene_end, alpha=0.08, color=GENE_COLORS[gene], zorder=0)
    ax3.set_ylim(-0.1, 1.1)
    ax3.set_yticks([0.2, 0.7])
    ax3.set_yticklabels(['CADD', 'AG q≥.95'], fontsize=6)
    plt.setp(ax3.get_xticklabels(), visible=False)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)

    # ── Track 4: Haplotype bar ──
    ax4 = fig.add_subplot(inner[3], sharex=ax1)
    if len(g_snvs) > 0:
        h1 = g_snvs[g_snvs['hap_alt'] == 'hap1']
        h2 = g_snvs[g_snvs['hap_alt'] == 'hap2']
        ax4.scatter(h1['pos'].values, np.zeros(len(h1)) + 0.3, c='#4393C3',
                    s=3, alpha=0.5, edgecolors='none', label='Hap1 alt')
        ax4.scatter(h2['pos'].values, np.zeros(len(h2)) + 0.7, c='#D6604D',
                    s=3, alpha=0.5, edgecolors='none', label='Hap2 alt')

    ax4.axvspan(gene_start, gene_end, alpha=0.08, color=GENE_COLORS[gene], zorder=0)
    ax4.set_ylim(0, 1)
    ax4.set_yticks([0.3, 0.7])
    ax4.set_yticklabels(['H1', 'H2'], fontsize=6)
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)

    if row_idx == n_genes - 1:
        ax4.set_xlabel('Genomic position (hg38)', fontsize=9)
    else:
        plt.setp(ax4.get_xticklabels(), visible=False)

    # Separator line at bottom
    if row_idx < n_genes - 1:
        ax4.axhline(y=-0.05, color='#cccccc', linewidth=0.5,
                     clip_on=False, xmin=-0.05, xmax=1.05)


def main():
    snvs, ag = load_data()

    fig = plt.figure(figsize=(18, 22))
    outer = fig.add_gridspec(len(GENES), 1, hspace=0.35,
                              left=0.10, right=0.96, top=0.94, bottom=0.03)

    for i, gene in enumerate(GENES):
        draw_gene_tracks(fig, outer[i], gene, snvs, ag, i, len(GENES))

    fig.suptitle('Variant Landscape — UCSC Browser-Style Tracks (±500 kb, hg38)',
                 fontsize=15, fontweight='bold', y=0.97)

    # Global legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#aaa', markersize=5,
               label='CADD < 10'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#ffd92f', markersize=5,
               label='CADD 10–15'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#ff7f0e', markersize=5,
               label='CADD 15–20'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#d62728', markersize=5,
               label='CADD ≥ 20'),
        Line2D([0], [0], marker='D', color='w', markerfacecolor='#d62728', markersize=6,
               markeredgecolor='black', markeredgewidth=0.3,
               label='AlphaGenome q ≥ 0.95'),
        Line2D([0], [0], marker='*', color='w', markerfacecolor='red', markersize=10,
               markeredgecolor='black', markeredgewidth=0.5,
               label='Known GWAS/literature hit'),
        mpatches.Patch(facecolor='#4393C3', alpha=0.5, label='Hap1 carries alt'),
        mpatches.Patch(facecolor='#D6604D', alpha=0.5, label='Hap2 carries alt'),
    ]
    fig.legend(handles=legend_elements, loc='upper right', fontsize=8,
               ncol=2, framealpha=0.9, bbox_to_anchor=(0.96, 0.96))

    for fmt in ['png', 'pdf']:
        fig.savefig(f'{OUT_FIG}.{fmt}', dpi=300, bbox_inches='tight')
        print(f'Saved: {OUT_FIG}.{fmt}')
    plt.close()


if __name__ == '__main__':
    main()
