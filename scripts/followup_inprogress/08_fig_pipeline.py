#!/usr/bin/env python3
"""
Pipeline/methods figure — biologist-friendly overview.
Shows the analysis workflow from genome assembly to regulatory variant candidates.
"""

import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
import matplotlib.patches as mpatches

plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 8,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

BASE = os.path.expanduser('~/reganalysis')

fig, ax = plt.subplots(figsize=(7.5, 10))
ax.set_xlim(0, 10)
ax.set_ylim(0, 14)
ax.axis('off')

# ── Colors ──
INPUT_COLOR = '#E3F2FD'
PROCESS_COLOR = '#FFF3E0'
OUTPUT_COLOR = '#E8F5E9'
VALIDATE_COLOR = '#FCE4EC'
RESULT_COLOR = '#F3E5F5'

def draw_box(ax, x, y, w, h, text, color, fontsize=8, bold=False, subtext=None):
    """Draw a rounded box with text."""
    box = FancyBboxPatch((x, y), w, h,
                         boxstyle="round,pad=0.15",
                         facecolor=color, edgecolor='#333333',
                         linewidth=0.8)
    ax.add_patch(box)
    weight = 'bold' if bold else 'normal'
    ax.text(x + w/2, y + h/2 + (0.12 if subtext else 0),
           text, ha='center', va='center', fontsize=fontsize,
           fontweight=weight, wrap=True)
    if subtext:
        ax.text(x + w/2, y + h/2 - 0.18, subtext,
               ha='center', va='center', fontsize=6, color='#555555',
               style='italic')

def draw_arrow(ax, x1, y1, x2, y2, color='#333333'):
    """Draw an arrow between boxes."""
    arrow = FancyArrowPatch((x1, y1), (x2, y2),
                            arrowstyle='->', mutation_scale=12,
                            color=color, linewidth=1.2)
    ax.add_patch(arrow)

def draw_validation(ax, x, y, text):
    """Draw a small validation checkmark box."""
    box = FancyBboxPatch((x, y), 2.2, 0.35,
                         boxstyle="round,pad=0.08",
                         facecolor=VALIDATE_COLOR, edgecolor='#C62828',
                         linewidth=0.6, linestyle='--')
    ax.add_patch(box)
    ax.text(x + 1.1, y + 0.175, text, ha='center', va='center',
           fontsize=5.5, color='#C62828')


# ════════════════════════════════════════════════════════════════
# TITLE
# ════════════════════════════════════════════════════════════════
ax.text(5, 13.6, 'Allele-Specific Regulatory Variant Discovery Pipeline',
       ha='center', fontsize=12, fontweight='bold')
ax.text(5, 13.25, 'H9 hESC — T2T Diploid Genome (t2t_h9_v01)',
       ha='center', fontsize=9, color='#555555')


# ════════════════════════════════════════════════════════════════
# ROW 1: Inputs (top)
# ════════════════════════════════════════════════════════════════
# T2T diploid genome
draw_box(ax, 0.3, 12.0, 2.8, 0.8,
         'T2T H9 Diploid\nGenome (Verkko)', INPUT_COLOR,
         fontsize=8, bold=True, subtext='hap1.fa + hap2.fa')

# RNA-seq data
draw_box(ax, 3.6, 12.0, 2.8, 0.8,
         'RNA-seq\n4 cell types', INPUT_COLOR,
         fontsize=8, bold=True, subtext='iPSC, FPP, astro, micro')

# Gene annotations
draw_box(ax, 6.9, 12.0, 2.8, 0.8,
         'CAT Gene\nAnnotations', INPUT_COLOR,
         fontsize=8, bold=True, subtext='H9_HAP1/2.gtf (UCSC)')


# ════════════════════════════════════════════════════════════════
# ROW 2: FASTA extraction + alignment
# ════════════════════════════════════════════════════════════════
draw_arrow(ax, 1.7, 12.0, 1.7, 11.3, '#2166AC')
draw_box(ax, 0.3, 10.4, 4.2, 0.85,
         'Extract syntenic gene regions (1Mb)\nfrom each haplotype', PROCESS_COLOR,
         fontsize=8, bold=True, subtext='samtools faidx — 5 gene loci x 2 haplotypes')

draw_validation(ax, 4.8, 10.55, 'Verify: gene body inside\nextraction window')

draw_arrow(ax, 2.4, 10.4, 2.4, 9.7)

draw_box(ax, 0.3, 8.85, 4.2, 0.8,
         'Align hap1 vs hap2\n(Winnowmap2 asm5)', PROCESS_COLOR,
         fontsize=8, bold=True, subtext='CIGAR string records every base difference')

draw_validation(ax, 4.8, 8.95, 'Check: aligned fraction >95%\nTs/Tv ratio ~2.0')

# RNA-seq path
draw_arrow(ax, 5.0, 12.0, 5.0, 9.7)
draw_box(ax, 5.5, 8.85, 4.2, 0.8,
         'Allele-specific alignment\n(UCSC .ase.bam)', PROCESS_COLOR,
         fontsize=8, bold=True, subtext='Reads assigned to hap1 or hap2 by UCSC')


# ════════════════════════════════════════════════════════════════
# ROW 3: Variant calling + ASE counting
# ════════════════════════════════════════════════════════════════
draw_arrow(ax, 2.4, 8.85, 2.4, 8.15)

draw_box(ax, 0.3, 7.3, 4.2, 0.8,
         'Parse CIGAR → call het variants\nbetween haplotypes', PROCESS_COLOR,
         fontsize=8, bold=True, subtext='4,515 variants (SNVs + indels)')

draw_validation(ax, 4.8, 7.4, 'FASTA spot-check: base at\nabs position matches allele')

draw_arrow(ax, 7.6, 8.85, 7.6, 8.15)

draw_box(ax, 5.5, 7.3, 4.2, 0.8,
         'Count allele-specific reads\nat each het variant position', PROCESS_COLOR,
         fontsize=8, bold=True, subtext='samtools view -c per hap per cell type')


# ════════════════════════════════════════════════════════════════
# ROW 4: Coordinate liftover + integration
# ════════════════════════════════════════════════════════════════
draw_arrow(ax, 2.4, 7.3, 2.4, 6.6)
draw_arrow(ax, 7.6, 7.3, 5.8, 6.6)

draw_box(ax, 0.8, 5.65, 4.5, 0.9,
         'Liftover to hg38 via\nUCSC diploid VCF (dipcall)', PROCESS_COLOR,
         fontsize=8, bold=True, subtext='Match T2T variants → hg38 positions\n2,717 SNVs with validated hg38 coords')

draw_validation(ax, 5.6, 5.95, 'Cross-validate: UCSC dip VCF\nmatch rate 86-99%')
draw_validation(ax, 5.6, 5.5, 'PAF-verified offsets (NOT\ngene_body - 500kb)')


# ════════════════════════════════════════════════════════════════
# ROW 5: AlphaGenome scoring
# ════════════════════════════════════════════════════════════════
draw_arrow(ax, 3.05, 5.65, 3.05, 5.0)

draw_box(ax, 0.8, 4.1, 4.5, 0.85,
         'AlphaGenome regulatory\nvariant effect prediction', PROCESS_COLOR,
         fontsize=8, bold=True, subtext='Cell-type-specific scores (astro, mono, NPC, H9)\n2,631 SNVs scored')

draw_box(ax, 5.8, 4.1, 3.9, 0.85,
         'Evo2 genomic\nlanguage model', PROCESS_COLOR,
         fontsize=8, bold=True, subtext='Log-likelihood ratio: ref vs alt\n(pending GPU)')


# ════════════════════════════════════════════════════════════════
# ROW 6: Integration
# ════════════════════════════════════════════════════════════════
draw_arrow(ax, 3.05, 4.1, 5.0, 3.3)
draw_arrow(ax, 7.75, 4.1, 5.0, 3.3)

draw_box(ax, 2.3, 2.4, 5.4, 0.85,
         'Integrate: ASE + AlphaGenome + Evo2\n→ Priority score per variant', OUTPUT_COLOR,
         fontsize=9, bold=True, subtext='Cell-type-specific regulatory candidates')


# ════════════════════════════════════════════════════════════════
# ROW 7: Key result
# ════════════════════════════════════════════════════════════════
draw_arrow(ax, 5.0, 2.4, 5.0, 1.7)

draw_box(ax, 1.5, 0.6, 7.0, 1.0,
         'SNCA: Astrocyte-specific hap2 bias (23:56 reads)\n'
         'Top candidate: chr4:89,836,694 C>T (AlphaGenome q=0.999)',
         RESULT_COLOR, fontsize=9, bold=True,
         subtext='T allele (hap2) associated with increased SNCA expression in astrocytes')


# ── Legend ──
legend_items = [
    mpatches.Patch(facecolor=INPUT_COLOR, edgecolor='#333', label='Input data'),
    mpatches.Patch(facecolor=PROCESS_COLOR, edgecolor='#333', label='Processing step'),
    mpatches.Patch(facecolor=OUTPUT_COLOR, edgecolor='#333', label='Integrated output'),
    mpatches.Patch(facecolor=RESULT_COLOR, edgecolor='#333', label='Key finding'),
    mpatches.Patch(facecolor=VALIDATE_COLOR, edgecolor='#C62828',
                   linestyle='--', label='Validation checkpoint'),
]
ax.legend(handles=legend_items, loc='lower left', fontsize=7,
         frameon=True, fancybox=True, framealpha=0.9,
         bbox_to_anchor=(0.0, -0.02))


# ── Save ──
outdir = f'{BASE}/clean_analysis/figures'
fig.savefig(f'{outdir}/Fig_Pipeline_Methods.png', dpi=300)
fig.savefig(f'{outdir}/Fig_Pipeline_Methods.pdf')
print(f"Saved: {outdir}/Fig_Pipeline_Methods.png")
plt.close()
