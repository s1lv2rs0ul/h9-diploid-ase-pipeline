#!/usr/bin/env python3
"""
AlphaGenome Results Figure — Multi-panel summary of variant effect predictions.

Panels:
  A: Manhattan-style plot of variant scores along each gene locus
  B: TF binding disruption — which transcription factors are affected per gene
  C: Regulatory mechanism breakdown (% variants by output type)
  D: Haplotype regulatory burden vs ASE direction
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import csv as csv_mod
import warnings
warnings.filterwarnings('ignore')

# ── Paths ──
BASE = '/Users/shruti.patil/reganalysis'
SUMMARY = f'{BASE}/variant_catalog/alphagenome/scores/alphagenome_variant_summary.csv'
BRAIN   = f'{BASE}/variant_catalog/alphagenome/scores/alphagenome_brain_scores.csv'
TOP     = f'{BASE}/variant_catalog/alphagenome/scores/alphagenome_top_effects.csv'
ASE_CSV = f'{BASE}/ase_v2/ase_all_genes_complete.csv'
HET_CSV = f'{BASE}/vep_redo/het_variants_hg38_official.csv'
OUT_FIG = f'{BASE}/variant_catalog/alphagenome/figures/Fig_AlphaGenome_results'

# ── Colors (Okabe-Ito) ──
GENE_COLORS = {
    'SNCA':  '#E69F00',
    'HTT':   '#56B4E9',
    'LRRK2': '#009E73',
    'GBA':   '#F0E442',
    'SMN1':  '#CC79A7',
}
GENES = ['SNCA', 'HTT', 'LRRK2', 'GBA', 'SMN1']

GENE_BODIES_HG38 = {
    'SNCA':  (89724099, 89838315),
    'HTT':   (3074681, 3243960),
    'LRRK2': (40224977, 40369285),
    'GBA':   (155234452, 155244699),
    'SMN1':  (70924941, 70953015),
}


def load_het_hap_alt():
    het = {}
    with open(HET_CSV) as f:
        for row in csv_mod.DictReader(f):
            ref, alt = row['ref'].upper(), row['alt'].upper()
            if len(ref) == 1 and len(alt) == 1 and ',' not in alt:
                vid = f"{row['chrom']}_{row['pos']}_{ref}_{alt}"
                het[vid] = row['hap_alt']
    return het


def load_data():
    summary = pd.read_csv(SUMMARY)
    brain = pd.read_csv(BRAIN)
    top = pd.read_csv(TOP)
    ase = pd.read_csv(ASE_CSV)
    het = load_het_hap_alt()
    summary['hap_alt'] = summary['our_variant_id'].map(het)
    return summary, brain, top, ase


def panel_a_manhattan(ax, summary):
    """Manhattan-style: variant position vs quantile score."""
    for gene in GENES:
        g = summary[summary['gene'] == gene].copy()
        if g.empty:
            continue
        g['pos'] = g['our_variant_id'].str.split('_').str[1].astype(int)
        gb = GENE_BODIES_HG38[gene]
        center = (gb[0] + gb[1]) / 2
        g['rel_pos_kb'] = (g['pos'] - center) / 1000
        ax.scatter(g['rel_pos_kb'], g['max_quantile'],
                   c=GENE_COLORS[gene], s=8, alpha=0.5,
                   label=gene, edgecolors='none', zorder=2)

    ax.axhline(y=0.95, color='red', linestyle='--', linewidth=0.8, alpha=0.6)
    ax.axvspan(-15, 15, alpha=0.05, color='grey')
    ax.set_xlabel('Distance from gene center (kb)', fontsize=10)
    ax.set_ylabel('Max quantile score', fontsize=10)
    ax.set_title('A. Variant effect scores across gene loci',
                 fontsize=11, fontweight='bold', loc='left')
    ax.legend(fontsize=7, ncol=3, loc='lower right', framealpha=0.8)
    ax.set_ylim(0.92, 1.005)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def panel_b_tf_disruption(ax, top):
    """Heatmap: genes × top disrupted TFs."""
    tf = top[top['output_type'] == 'CHIP_TF'].copy()
    tf = tf.dropna(subset=['transcription_factor'])

    # Get top TFs across all genes
    tf_counts = tf.groupby(['gene', 'transcription_factor'])['our_variant_id'].nunique()
    tf_counts = tf_counts.reset_index(name='n_variants')

    # Top 10 TFs by total disrupted variants
    top_tfs = tf_counts.groupby('transcription_factor')['n_variants'].sum() \
                        .nlargest(10).index.tolist()

    # Pivot
    pivot = tf_counts[tf_counts['transcription_factor'].isin(top_tfs)].pivot_table(
        values='n_variants', index='gene', columns='transcription_factor', fill_value=0
    )
    pivot = pivot.reindex(index=GENES, columns=top_tfs, fill_value=0)

    im = ax.imshow(pivot.values, cmap='Blues', aspect='auto')
    ax.set_xticks(range(len(top_tfs)))
    ax.set_xticklabels(top_tfs, rotation=45, ha='right', fontsize=8)
    ax.set_yticks(range(len(GENES)))
    ax.set_yticklabels(GENES, fontsize=9)
    ax.set_title('B. TF binding sites disrupted (unique variants per TF)',
                 fontsize=11, fontweight='bold', loc='left')

    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            val = pivot.values[i, j]
            if val > 0:
                color = 'white' if val > pivot.values.max() * 0.6 else 'black'
                ax.text(j, i, str(int(val)), ha='center', va='center',
                        fontsize=7, color=color, fontweight='bold')

    plt.colorbar(im, ax=ax, label='# variants disrupting TF', shrink=0.8, pad=0.02)


def panel_c_mechanism_breakdown(ax, top):
    """Stacked bar: % of variants affecting each regulatory mechanism."""
    output_order = ['RNA_SEQ', 'CAGE', 'ATAC', 'DNASE', 'CHIP_HISTONE', 'CHIP_TF']
    output_labels = ['RNA-seq', 'CAGE', 'ATAC', 'DNase', 'Histone ChIP', 'TF ChIP']
    output_colors = {
        'RNA_SEQ': '#1f77b4', 'CAGE': '#ff7f0e', 'ATAC': '#2ca02c',
        'DNASE': '#d62728', 'CHIP_HISTONE': '#9467bd', 'CHIP_TF': '#8c564b',
    }

    # Count unique variants per gene per output type
    gene_totals = {gene: top[top['gene'] == gene]['our_variant_id'].nunique() for gene in GENES}
    counts = top.groupby(['gene', 'output_type'])['our_variant_id'].nunique().unstack(fill_value=0)
    counts = counts.reindex(index=GENES, columns=output_order, fill_value=0)

    # Convert to percentage of scored variants (3049 total per gene from summary)
    n_scored = {gene: len(pd.read_csv(SUMMARY).query(f"gene=='{gene}'")) for gene in GENES}
    pct = counts.copy()
    for gene in GENES:
        if n_scored[gene] > 0:
            pct.loc[gene] = (counts.loc[gene] / n_scored[gene]) * 100

    x = np.arange(len(GENES))
    width = 0.6
    bottom = np.zeros(len(GENES))

    for oi, otype in enumerate(output_order):
        vals = pct[otype].values
        ax.bar(x, vals, width, bottom=bottom, label=output_labels[oi],
               color=output_colors[otype], edgecolor='white', linewidth=0.5)
        bottom += vals

    ax.set_xticks(x)
    ax.set_xticklabels(GENES, fontsize=9)
    ax.set_ylabel('% of SNVs with high-impact score', fontsize=10)
    ax.set_title('C. Regulatory mechanisms affected (% of SNVs, q ≥ 0.9)',
                 fontsize=11, fontweight='bold', loc='left')
    ax.legend(fontsize=7, ncol=2, loc='upper right', framealpha=0.8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def panel_d_haplotype_bias(ax, summary, ase):
    """Haplotype regulatory burden vs observed ASE direction."""
    ase_astro = ase[ase['CellType'] == 'astrocyte'].copy()
    ase_astro['gene_key'] = ase_astro['Gene'].replace({'GBA1': 'GBA'})

    x = np.arange(len(GENES))
    width = 0.35

    h1_burdens, h2_burdens, ase_h1, ase_h2 = [], [], [], []

    for gene in GENES:
        g = summary[summary['gene'] == gene]
        h1_burdens.append(g[g['hap_alt'] == 'hap1']['max_raw'].sum())
        h2_burdens.append(g[g['hap_alt'] == 'hap2']['max_raw'].sum())
        row = ase_astro[ase_astro['gene_key'] == gene]
        ase_h1.append(int(row.iloc[0]['Hap1_ASE']) if not row.empty else 0)
        ase_h2.append(int(row.iloc[0]['Hap2_ASE']) if not row.empty else 0)

    h1_arr, h2_arr = np.array(h1_burdens), np.array(h2_burdens)
    total_burden = h1_arr + h2_arr
    h1_frac = h1_arr / total_burden
    h2_frac = h2_arr / total_burden

    ax.bar(x - width/2, h1_frac, width, label='Hap1 alt burden',
           color='#4393C3', edgecolor='white', linewidth=0.5)
    ax.bar(x + width/2, h2_frac, width, label='Hap2 alt burden',
           color='#D6604D', edgecolor='white', linewidth=0.5)
    ax.axhline(y=0.5, color='grey', linestyle='--', linewidth=0.8, alpha=0.5)

    # Binomial test p-values
    sig_genes = {'SNCA': 0.0003, 'SMN1': 0.041}

    for i, gene in enumerate(GENES):
        total_ase = ase_h1[i] + ase_h2[i]
        if total_ase > 0:
            is_sig = gene in sig_genes
            if is_sig:
                direction = 'H2↑' if ase_h2[i] > ase_h1[i] else 'H1↑'
                color = '#D6604D' if ase_h2[i] > ase_h1[i] else '#4393C3'
                p = sig_genes[gene]
                label = f'ASE: {direction}*\n({ase_h1[i]}/{ase_h2[i]})\np={p:.4f}'
                ax.annotate(label, xy=(i, max(h1_frac[i], h2_frac[i]) + 0.02),
                            ha='center', va='bottom', fontsize=7,
                            fontweight='bold', color=color,
                            bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                                      edgecolor=color, alpha=0.8))
            else:
                label = f'ASE: n.s.\n({ase_h1[i]}/{ase_h2[i]})'
                ax.annotate(label, xy=(i, max(h1_frac[i], h2_frac[i]) + 0.02),
                            ha='center', va='bottom', fontsize=7, color='grey',
                            bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                                      edgecolor='grey', alpha=0.6))

    ax.set_xticks(x)
    ax.set_xticklabels(GENES, fontsize=9)
    ax.set_ylabel('Fraction of regulatory burden', fontsize=10)
    ax.set_title('D. Haplotype regulatory burden vs ASE (astrocyte)',
                 fontsize=11, fontweight='bold', loc='left')
    ax.legend(fontsize=8, loc='upper right', framealpha=0.8)
    ax.set_ylim(0, 0.82)
    ax.annotate('* p < 0.05 (binomial test)', xy=(0.02, 0.97),
                xycoords='axes fraction', fontsize=7, fontstyle='italic', color='grey')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def main():
    summary, brain, top, ase = load_data()

    fig = plt.figure(figsize=(16, 14))
    gs = gridspec.GridSpec(2, 2, hspace=0.4, wspace=0.3,
                           left=0.07, right=0.95, top=0.93, bottom=0.08)

    panel_a_manhattan(fig.add_subplot(gs[0, 0]), summary)
    panel_b_tf_disruption(fig.add_subplot(gs[0, 1]), top)
    panel_c_mechanism_breakdown(fig.add_subplot(gs[1, 0]), top)
    panel_d_haplotype_bias(fig.add_subplot(gs[1, 1]), summary, ase)

    fig.suptitle('AlphaGenome Variant Effect Predictions — H9 T2T Diploid Het SNVs',
                 fontsize=14, fontweight='bold', y=0.97)

    for fmt in ['png', 'pdf']:
        fig.savefig(f'{OUT_FIG}.{fmt}', dpi=300, bbox_inches='tight')
        print(f'Saved: {OUT_FIG}.{fmt}')
    plt.close()


if __name__ == '__main__':
    main()
