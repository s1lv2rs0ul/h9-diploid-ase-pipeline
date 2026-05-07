#!/usr/bin/env python3
"""
Variant Annotation Figure — Known effects, CADD scores, consequence types.

Panels:
  A: CADD score distribution per gene with thresholds
  B: Variant consequences (coding, UTR, regulatory, intronic)
  C: Known GWAS/literature variants highlighted on gene map
  D: Haplotype burden combining AlphaGenome + CADD + VEP
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch
import warnings
warnings.filterwarnings('ignore')

BASE = '/Users/shruti.patil/reganalysis'
VEP_CSV  = f'{BASE}/vep_redo/vep_het_official_hg38.csv'
AG_CSV   = f'{BASE}/variant_catalog/alphagenome/scores/alphagenome_variant_summary.csv'
ASE_CSV  = f'{BASE}/ase_v2/ase_all_genes_complete.csv'
OUT_FIG  = f'{BASE}/variant_catalog/alphagenome/figures/Fig_Variant_Annotation'

GENE_COLORS = {
    'SNCA': '#E69F00', 'HTT': '#56B4E9', 'LRRK2': '#009E73',
    'GBA': '#F0E442', 'SMN1': '#CC79A7',
}
GENES = ['SNCA', 'HTT', 'LRRK2', 'GBA', 'SMN1']

# Known PD-associated variants found in H9
KNOWN_GWAS = {
    'rs356219':   {'gene': 'SNCA',  'desc': 'PD GWAS lead SNP\n(Nalls 2019)'},
    'rs356220':   {'gene': 'SNCA',  'desc': 'PD risk variant'},
    'rs2583988':  {'gene': 'SNCA',  'desc': 'PD risk\n(ncRNA exon)'},
    'rs35303786': {'gene': 'LRRK2', 'desc': 'Missense\n(CADD=21.4)'},
    'rs2638869':  {'gene': 'LRRK2', 'desc': 'Stop gained\n(CADD=29.8)'},
    'rs11564230': {'gene': 'LRRK2', 'desc': 'Splice donor\n(CADD=24.5)'},
    'rs11549516': {'gene': 'HTT',   'desc': 'Missense\n(CADD=17.3)'},
}

GENE_BODIES_HG38 = {
    'SNCA':  (89724099, 89838315, 'chr4'),
    'HTT':   (3074681, 3243960, 'chr4'),
    'LRRK2': (40224977, 40369285, 'chr12'),
    'GBA':   (155234452, 155244699, 'chr1'),
    'SMN1':  (70924941, 70953015, 'chr5'),
}

# Approximate exon boundaries (hg38) for coding/non-coding annotation
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


def is_in_exon(pos, gene):
    """Check if a genomic position falls within an exon."""
    for ex_start, ex_end in EXONS_HG38.get(gene, []):
        if ex_start <= pos <= ex_end:
            return True
    return False


def load_data():
    vep = pd.read_csv(VEP_CSV)
    snvs = vep[(vep['ref'].str.len() == 1) & (vep['alt'].str.len() == 1) &
               (~vep['alt'].str.contains(',', na=False))].copy()
    ag = pd.read_csv(AG_CSV)
    return snvs, ag


def panel_a_cadd_distribution(ax, snvs):
    """Violin/box plot of CADD scores per gene."""
    data = []
    positions = []
    colors = []

    for i, gene in enumerate(GENES):
        g = snvs[(snvs['gene'] == gene) & snvs['cadd_phred'].notna()]
        if len(g) > 0:
            data.append(g['cadd_phred'].values)
            positions.append(i)
            colors.append(GENE_COLORS[gene])

    bp = ax.boxplot(data, positions=positions, widths=0.5, patch_artist=True,
                    showfliers=True, flierprops={'markersize': 2, 'alpha': 0.3})

    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    # Threshold lines
    ax.axhline(y=15, color='orange', linestyle='--', linewidth=1, alpha=0.7,
               label='CADD ≥ 15 (likely deleterious)')
    ax.axhline(y=20, color='red', linestyle='--', linewidth=1, alpha=0.7,
               label='CADD ≥ 20 (highly deleterious)')

    # Annotate top variants
    for gene in GENES:
        g = snvs[(snvs['gene'] == gene) & (snvs['cadd_phred'] > 20)]
        gi = GENES.index(gene)
        for _, r in g.iterrows():
            ax.annotate(f'{r["rsid"]}',
                        xy=(gi + 0.3, r['cadd_phred']),
                        fontsize=5.5, color='red', alpha=0.8)

    ax.set_xticks(range(len(GENES)))
    ax.set_xticklabels(GENES, fontsize=9)
    ax.set_ylabel('CADD Phred Score', fontsize=10)
    ax.set_title('A. CADD deleteriousness scores per gene',
                 fontsize=11, fontweight='bold', loc='left')
    ax.legend(fontsize=7, loc='upper left', framealpha=0.8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def panel_b_consequences(ax, snvs):
    """Stacked bar: variant consequences per gene."""
    # Group consequences
    def classify(c):
        if pd.isna(c):
            return 'Other'
        if 'missense' in c or 'stop' in c or 'frameshift' in c:
            return 'Coding'
        if 'splice' in c:
            return 'Splice'
        if 'UTR' in c:
            return 'UTR'
        if 'non_coding_transcript' in c:
            return 'ncRNA'
        if 'regulatory' in c:
            return 'Regulatory'
        if 'upstream' in c or 'downstream' in c:
            return 'Up/Downstream'
        if 'intron' in c:
            return 'Intronic'
        if 'intergenic' in c:
            return 'Intergenic'
        return 'Other'

    snvs['class'] = snvs['consequence'].apply(classify)

    cat_order = ['Coding', 'Splice', 'UTR', 'ncRNA', 'Up/Downstream', 'Intronic', 'Intergenic', 'Other']
    cat_colors = {
        'Coding': '#d62728', 'Splice': '#ff7f0e', 'UTR': '#2ca02c',
        'ncRNA': '#9467bd', 'Up/Downstream': '#8c564b',
        'Intronic': '#1f77b4', 'Intergenic': '#bcbd22', 'Other': '#999999',
    }

    counts = snvs.groupby(['gene', 'class']).size().unstack(fill_value=0)
    counts = counts.reindex(index=GENES, columns=cat_order, fill_value=0)

    # Convert to percentage
    pct = counts.div(counts.sum(axis=1), axis=0) * 100

    x = np.arange(len(GENES))
    width = 0.6
    bottom = np.zeros(len(GENES))

    for cat in cat_order:
        if cat in pct.columns:
            vals = pct[cat].values
            ax.bar(x, vals, width, bottom=bottom, label=cat,
                   color=cat_colors.get(cat, '#999'), edgecolor='white', linewidth=0.5)
            # Label if > 5%
            for i, v in enumerate(vals):
                if v > 8:
                    ax.text(i, bottom[i] + v / 2, f'{v:.0f}%',
                            ha='center', va='center', fontsize=6, color='white',
                            fontweight='bold')
            bottom += vals

    ax.set_xticks(x)
    ax.set_xticklabels(GENES, fontsize=9)
    ax.set_ylabel('% of SNVs', fontsize=10)
    ax.set_title('B. Variant consequence classes',
                 fontsize=11, fontweight='bold', loc='left')
    ax.legend(fontsize=7, ncol=2, loc='upper right', framealpha=0.8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def panel_c_gwas_map(ax, snvs):
    """Gene locus map highlighting known GWAS/literature variants with coding/non-coding."""
    gene_y = {g: i for i, g in enumerate(reversed(GENES))}

    # Draw gene bodies with exon/intron structure
    for gene in GENES:
        y = gene_y[gene]
        gb = GENE_BODIES_HG38[gene]
        center = (gb[0] + gb[1]) / 2

        # Background region (±500kb)
        ax.barh(y, 1000, left=-500, height=0.12, color='#e8e8e8', zorder=1)

        # Intron line (thin)
        body_left = (gb[0] - center) / 1000
        body_right = (gb[1] - center) / 1000
        ax.plot([body_left, body_right], [y, y], color='#555555',
                linewidth=1.5, zorder=2)

        # Exons as thick blocks (coding regions)
        for ex_start, ex_end in EXONS_HG38[gene]:
            ex_left = (ex_start - center) / 1000
            ex_width = (ex_end - ex_start) / 1000
            ax.barh(y, ex_width, left=ex_left, height=0.35,
                    color=GENE_COLORS[gene], alpha=0.8, zorder=3,
                    edgecolor='#333333', linewidth=0.5)

        ax.text(-520, y, gene, ha='right', va='center', fontsize=10, fontweight='bold')

        # Label coding/non-coding regions
        ax.text(body_right + 5, y + 0.15, 'exons', fontsize=5.5, color='#555',
                va='center', fontstyle='italic')

    # Plot all variants — color by coding vs non-coding
    for gene in GENES:
        y = gene_y[gene]
        g = snvs[snvs['gene'] == gene].copy()
        gb = GENE_BODIES_HG38[gene]
        center = (gb[0] + gb[1]) / 2
        g['rel_kb'] = (g['pos'] - center) / 1000
        g['cadd'] = g['cadd_phred'].fillna(0)
        g['in_exon'] = g['pos'].apply(lambda p: is_in_exon(p, gene))

        # Non-coding variants (outside exons)
        nc = g[~g['in_exon']]
        nc_low = nc[nc['cadd'] < 15]
        nc_high = nc[nc['cadd'] >= 15]

        ax.scatter(nc_low['rel_kb'], [y - 0.28] * len(nc_low), s=3,
                   c='#b0b0b0', alpha=0.3, zorder=3, marker='o')
        if len(nc_high) > 0:
            ax.scatter(nc_high['rel_kb'], [y - 0.28] * len(nc_high), s=15,
                       c='#ff7f0e', alpha=0.7, zorder=4, edgecolors='black',
                       linewidth=0.3, marker='o')

        # Coding variants (inside exons)
        cd = g[g['in_exon']]
        cd_low = cd[cd['cadd'] < 15]
        cd_high = cd[cd['cadd'] >= 15]

        ax.scatter(cd_low['rel_kb'], [y - 0.28] * len(cd_low), s=12,
                   c='#2ca02c', alpha=0.6, zorder=5, marker='s')
        if len(cd_high) > 0:
            ax.scatter(cd_high['rel_kb'], [y - 0.28] * len(cd_high), s=25,
                       c='#d62728', alpha=0.9, zorder=6, edgecolors='black',
                       linewidth=0.5, marker='s')

    # Highlight GWAS variants with labels
    for rsid, info in KNOWN_GWAS.items():
        gene = info['gene']
        match = snvs[snvs['rsid'] == rsid]
        if match.empty:
            continue
        r = match.iloc[0]
        y = gene_y[gene]
        gb = GENE_BODIES_HG38[gene]
        center = (gb[0] + gb[1]) / 2
        rel_kb = (r['pos'] - center) / 1000

        ax.scatter([rel_kb], [y + 0.25], s=60, c='red', marker='*',
                   zorder=7, edgecolors='black', linewidth=0.5)
        ax.annotate(f'{rsid}\n{info["desc"]}',
                    xy=(rel_kb, y + 0.25), xytext=(rel_kb, y + 0.55),
                    fontsize=5.5, ha='center', va='bottom', color='red',
                    arrowprops=dict(arrowstyle='-', color='red', lw=0.5),
                    bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                              edgecolor='red', alpha=0.8))

    ax.set_xlabel('Distance from gene center (kb)', fontsize=10)
    ax.set_title('C. Known disease-associated variants in H9 het sites (coding vs non-coding)',
                 fontsize=11, fontweight='bold', loc='left')
    ax.set_xlim(-550, 550)
    ax.set_ylim(-0.7, len(GENES) - 0.3)
    ax.set_yticks([])

    # Legend
    from matplotlib.lines import Line2D
    import matplotlib.patches as mpatches
    legend_elements = [
        mpatches.Patch(facecolor=GENE_COLORS['SNCA'], alpha=0.8,
                       edgecolor='#333', label='Exon (coding)'),
        Line2D([0], [0], color='#555', linewidth=1.5, label='Intron'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor='#2ca02c',
               markersize=5, label='Coding variant (CADD<15)'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor='#d62728',
               markersize=6, label='Coding variant (CADD≥15)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#b0b0b0',
               markersize=4, label='Non-coding (CADD<15)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#ff7f0e',
               markersize=5, label='Non-coding (CADD≥15)'),
        Line2D([0], [0], marker='*', color='w', markerfacecolor='red',
               markersize=10, label='GWAS/literature hit'),
    ]
    ax.legend(handles=legend_elements, fontsize=6, loc='lower right',
              framealpha=0.9, ncol=2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)


def panel_d_summary_table(ax, snvs):
    """Summary table of key findings."""
    ax.axis('off')
    ax.set_title('D. Summary of variant annotations',
                 fontsize=11, fontweight='bold', loc='left')

    # Build table data
    rows = []
    for gene in GENES:
        g = snvs[snvs['gene'] == gene]
        total = len(g)
        has_cadd = g[g['cadd_phred'].notna()]
        high_cadd = (has_cadd['cadd_phred'] >= 15).sum()
        coding = g[g['consequence'].str.contains('missense|stop|splice', na=False, regex=True)]
        utr = g[g['consequence'].str.contains('UTR', na=False)]
        reg = g[(g['n_regulatory'] > 0) | (g['n_motif'] > 0)]
        gwas = len([k for k, v in KNOWN_GWAS.items()
                    if v['gene'] == gene and not g[g['rsid'] == k].empty])

        rows.append([gene, str(total), str(len(coding)), str(len(utr)),
                     str(len(reg)), str(high_cadd), str(gwas)])

    col_labels = ['Gene', 'SNVs', 'Coding/\nSplice', 'UTR', 'Regulatory\nAnnotated',
                  'CADD≥15', 'Known\nGWAS']

    table = ax.table(cellText=rows, colLabels=col_labels,
                     cellLoc='center', loc='center',
                     bbox=[0.05, 0.1, 0.9, 0.8])

    table.auto_set_font_size(False)
    table.set_fontsize(9)

    # Style header
    for j in range(len(col_labels)):
        cell = table[0, j]
        cell.set_facecolor('#2C3E50')
        cell.set_text_props(color='white', fontweight='bold', fontsize=8)

    # Style gene column with colors
    for i, gene in enumerate(GENES):
        table[i + 1, 0].set_facecolor(GENE_COLORS[gene])
        table[i + 1, 0].set_text_props(fontweight='bold')
        # Highlight high values
        if int(rows[i][5]) > 0:  # CADD >= 15
            table[i + 1, 5].set_facecolor('#FFCCCC')
        if int(rows[i][6]) > 0:  # GWAS hits
            table[i + 1, 6].set_facecolor('#FFCCCC')

    for key, cell in table.get_celld().items():
        cell.set_edgecolor('#CCCCCC')


def main():
    snvs, ag = load_data()

    fig = plt.figure(figsize=(16, 14))
    gs = gridspec.GridSpec(2, 2, hspace=0.4, wspace=0.3,
                           left=0.07, right=0.95, top=0.93, bottom=0.06)

    panel_a_cadd_distribution(fig.add_subplot(gs[0, 0]), snvs)
    panel_b_consequences(fig.add_subplot(gs[0, 1]), snvs)
    panel_c_gwas_map(fig.add_subplot(gs[1, :]), snvs)

    fig.suptitle('Variant Annotation & Known Disease Associations — H9 Het SNVs',
                 fontsize=14, fontweight='bold', y=0.97)

    for fmt in ['png', 'pdf']:
        fig.savefig(f'{OUT_FIG}.{fmt}', dpi=300, bbox_inches='tight')
        print(f'Saved: {OUT_FIG}.{fmt}')
    plt.close()


if __name__ == '__main__':
    main()
