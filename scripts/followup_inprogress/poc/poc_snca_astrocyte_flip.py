#!/usr/bin/env python3
"""
═══════════════════════════════════════════════════════════════════════════════
PROOF OF CONCEPT — SNCA astrocyte ASE flip
═══════════════════════════════════════════════════════════════════════════════

Question:
    SNCA shows balanced hap1/hap2 expression in iPSC, FPP, microglia — but
    astrocyte-specifically flips to hap2-biased (23:56 reads, ratio 0.29).
    Can we identify the non-coding variant(s) responsible?

Data used (all already on disk, no new BAMs / no new AlphaGenome calls):
    - ase_v2/ase_all_genes_complete.csv           (gene-level ASE, 4 cell types)
    - clean_analysis/variants/integrated_variant_analysis.csv
                                                  (per-variant: T2T coords, hg38 coords,
                                                   ASE reads per cell type, AG scores
                                                   per cell type, priority_score)
    - Documents/T2T/variant_catalog/alphagenome/scores/alphagenome_brain_scores.csv
                                                  (per-tissue AG scores, 30,490 rows)

Analysis flow:
    Panel A — statistical validation: is SNCA astrocyte flip real?
        Binomial test on each cell type; Fisher's exact for astro vs pooled others.
    Panel B — AG astrocyte-specific candidates in SNCA gene body:
        Volcano-style plot of AG astrocyte raw_score vs quantile.
    Panel C — Cell-type specificity of top candidate:
        Bar plot of AG max_quantile in astrocyte / NPC / H9 / monocyte tracks.
    Panel D — Direction concordance for top candidate:
        Does predicted direction (AG raw_score sign × hap_alt) match observed
        hap2-bias in astrocyte?

Outputs:
    - results/followup_inprogress/poc/poc_snca_summary.txt   (numeric report)
    - results/followup_inprogress/poc/poc_snca_top_candidates.csv
    - figures/followup_inprogress/scratch/POC_SNCA_astrocyte_flip.{png,pdf}
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ══════════════════════════════════════════════════════════════
# PATHS
# ══════════════════════════════════════════════════════════════
HOME = os.path.expanduser('~')
GENE_ASE   = f'{HOME}/reganalysis/ase_v2/ase_all_genes_complete.csv'
INTEGRATED = f'{HOME}/reganalysis/clean_analysis/variants/integrated_variant_analysis.csv'
BRAIN_AG   = f'{HOME}/Documents/T2T/variant_catalog/alphagenome/scores/alphagenome_brain_scores.csv'

REPO       = f'{HOME}/h9-diploid-ase-pipeline'
OUT_DIR    = f'{REPO}/results/followup_inprogress/poc'
FIG_DIR    = f'{REPO}/figures/followup_inprogress/scratch'
SUMMARY    = f'{OUT_DIR}/poc_snca_summary.txt'
TOP_CSV    = f'{OUT_DIR}/poc_snca_top_candidates.csv'
FIG        = f'{FIG_DIR}/POC_SNCA_astrocyte_flip'

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

# ══════════════════════════════════════════════════════════════
# 1. Statistical validation of the SNCA astrocyte flip
# ══════════════════════════════════════════════════════════════
def validate_astrocyte_flip():
    """Test whether SNCA ASE in astrocyte differs from other cell types."""
    gene = pd.read_csv(GENE_ASE)
    snca = gene[gene['Gene'] == 'SNCA'].copy()
    snca['hap1_frac'] = snca['Hap1_ASE'] / snca['ASE_Total']

    # Per-cell-type binomial test vs H0: p(hap1) = 0.5
    per_ct_stats = []
    for _, r in snca.iterrows():
        h1, tot = int(r['Hap1_ASE']), int(r['ASE_Total'])
        res = stats.binomtest(h1, tot, p=0.5, alternative='two-sided')
        per_ct_stats.append({
            'CellType':  r['CellType'],
            'Hap1_ASE':  h1,
            'Hap2_ASE':  int(r['Hap2_ASE']),
            'Total':     tot,
            'hap1_frac': h1/tot,
            'p_vs_0.5':  res.pvalue,
            'signif':    'yes' if res.pvalue < 0.05 else 'no',
        })
    per_ct = pd.DataFrame(per_ct_stats)

    # Astro vs pooled others (Fisher's exact)
    astro  = snca[snca['CellType'] == 'astrocyte'].iloc[0]
    others = snca[snca['CellType'] != 'astrocyte']
    others_h1 = int(others['Hap1_ASE'].sum())
    others_h2 = int(others['Hap2_ASE'].sum())
    contingency = np.array([[int(astro['Hap1_ASE']), int(astro['Hap2_ASE'])],
                            [others_h1,             others_h2]])
    odds, pval = stats.fisher_exact(contingency, alternative='two-sided')
    astro_vs_others = {
        'astro_h1':   int(astro['Hap1_ASE']),
        'astro_h2':   int(astro['Hap2_ASE']),
        'others_h1':  others_h1,
        'others_h2':  others_h2,
        'odds_ratio': odds,
        'p_value':    pval,
    }
    return per_ct, astro_vs_others


# ══════════════════════════════════════════════════════════════
# 2. AlphaGenome candidates in SNCA — astrocyte-focused
# ══════════════════════════════════════════════════════════════
def snca_astrocyte_candidates():
    """Rank SNCA gene-body variants by AlphaGenome astrocyte signal."""
    df = pd.read_csv(INTEGRATED)
    snca = df[df['gene'] == 'SNCA'].copy()

    # Focus on gene-body variants with AG signal
    snca_gb = snca[(snca['in_gene_body']) &
                   (snca['has_alphagenome']) &
                   (snca['astrocyte_max_quantile'].notna())].copy()

    # Sort by astrocyte quantile
    snca_gb = snca_gb.sort_values('astrocyte_max_quantile', ascending=False)
    return snca_gb, snca


# ══════════════════════════════════════════════════════════════
# 3. Direction concordance for top candidate
# ══════════════════════════════════════════════════════════════
def direction_concordance(candidates):
    """
    Check whether AG-predicted direction of effect matches the observed
    astrocyte hap2 bias (56 hap2 vs 23 hap1 = hap2 up-regulated).

    Sign convention:
        - AG raw_score is the effect of the ALT allele on the tracked output.
          POSITIVE = ALT up-regulates.
        - hap_alt tells us which haplotype carries the ALT allele on hg38.
          (Determined during liftover step; encoded here via hap1_base/hap2_base
          vs hg38_ref/hg38_alt.)

    A candidate variant explains hap2-biased astrocyte expression if:
        - hap2 carries the ALT allele AND AG raw_score > 0  (predicts hap2 up)
        - OR hap1 carries the ALT allele AND AG raw_score < 0  (predicts hap1 down → hap2 relatively up)

    We compute a signed "predicted hap direction" per variant and compare to
    observed astrocyte read imbalance (astro_hap2_reads - astro_hap1_reads).
    """
    df = candidates.copy()

    def hap_alt(row):
        # Which haplotype carries the hg38 ALT allele?
        if pd.isna(row['hg38_alt']) or pd.isna(row['hap1_base']) or pd.isna(row['hap2_base']):
            return 'unknown'
        if row['hap1_base'] == row['hg38_alt']:
            return 'hap1'
        if row['hap2_base'] == row['hg38_alt']:
            return 'hap2'
        return 'both_differ'

    df['hap_alt'] = df.apply(hap_alt, axis=1)

    def pred_hap_direction(row):
        raw = row['astrocyte_max_raw']
        if pd.isna(raw) or row['hap_alt'] == 'unknown':
            return np.nan
        # Positive AG raw = ALT-allele up
        alt_up = raw > 0
        if row['hap_alt'] == 'hap1':
            # hap1 = ALT → alt_up means hap1_up = +1, alt_down means hap2_up = -1
            return +1 if alt_up else -1
        elif row['hap_alt'] == 'hap2':
            # hap2 = ALT → alt_up means hap2_up = -1 (in our +hap1 convention)
            return -1 if alt_up else +1
        else:
            return np.nan

    df['pred_hap1_direction'] = df.apply(pred_hap_direction, axis=1)  # +1 = hap1 up, -1 = hap2 up

    # Observed: sign of (hap1_reads - hap2_reads)
    df['obs_hap1_direction'] = np.sign(df['astro_hap1_reads'] - df['astro_hap2_reads'])

    # Concordance flag (only meaningful when both defined AND observed direction is non-zero)
    both_ok = df['pred_hap1_direction'].notna() & (df['obs_hap1_direction'] != 0)
    df.loc[both_ok, 'concordant'] = (df.loc[both_ok, 'pred_hap1_direction'] ==
                                     df.loc[both_ok, 'obs_hap1_direction'])
    return df


# ══════════════════════════════════════════════════════════════
# 4. Report
# ══════════════════════════════════════════════════════════════
def write_report(per_ct, astro_vs_others, snca_gb, snca_with_dir):
    lines = []
    lines.append("═" * 78)
    lines.append("PROOF OF CONCEPT — SNCA astrocyte ASE flip")
    lines.append("═" * 78)
    lines.append("")

    lines.append("── 1. Gene-level ASE (H9 T2T diploid, from paper Fig 3) ──────────────")
    lines.append(per_ct.to_string(index=False, float_format='%.4g'))
    lines.append("")
    lines.append(f"Astrocyte vs pooled other cell types (Fisher's exact):")
    for k, v in astro_vs_others.items():
        if isinstance(v, float):
            lines.append(f"    {k}: {v:.3g}")
        else:
            lines.append(f"    {k}: {v}")
    verdict = ("SIGNIFICANT — astrocyte hap2 bias is not seen in other lineages"
               if astro_vs_others['p_value'] < 0.05
               else "NOT SIGNIFICANT — cannot rule out ASE ratio matching other cell types")
    lines.append(f"    Verdict: {verdict}")
    lines.append("")

    lines.append("── 2. AlphaGenome candidates in SNCA gene body ────────────────────────")
    lines.append(f"Total SNCA variants with AG scores: {len(snca_gb)}")
    lines.append(f"SNCA gene-body variants scored: {snca_gb['in_gene_body'].sum()}")
    lines.append(f"    astrocyte_max_quantile ≥ 0.99: {(snca_gb['astrocyte_max_quantile']>=0.99).sum()}")
    lines.append(f"    astrocyte_max_quantile ≥ 0.999: {(snca_gb['astrocyte_max_quantile']>=0.999).sum()}")
    lines.append("")

    lines.append("── 3. Top 10 SNCA gene-body variants by astrocyte AG quantile ─────────")
    top = snca_gb[['variant_id','hg38_pos','hap1_base','hap2_base',
                   'astrocyte_max_quantile','astrocyte_max_raw',
                   'astro_hap1_reads','astro_hap2_reads',
                   'priority_score']].head(10)
    lines.append(top.to_string(index=False, float_format='%.4g'))
    lines.append("")

    lines.append("── 4. Direction concordance ──────────────────────────────────────────")
    dd = snca_with_dir[snca_with_dir['concordant'].notna()]
    if len(dd) > 0:
        n_conc = int(dd['concordant'].sum())
        n_total = len(dd)
        frac = n_conc / n_total if n_total else np.nan
        # Binomial test of concordance vs 50%
        conc_p = stats.binomtest(n_conc, n_total, p=0.5, alternative='greater').pvalue
        lines.append(f"Variants with both AG-direction and observed ASE-direction defined: {n_total}")
        lines.append(f"AG-predicted direction matches observed astrocyte direction: {n_conc}/{n_total} "
                     f"({frac:.1%})")
        lines.append(f"Binomial p(concordance > 50%): {conc_p:.3g}")
    else:
        lines.append("No variants have both AG-direction AND non-zero observed astrocyte reads.")
    lines.append("")

    lines.append("── 5. Headline candidate (verify against original pipeline) ───────────")
    headline_id = 'chr4_89836694_C_T'
    hl = snca_gb[snca_gb['variant_id'] == headline_id]
    if len(hl) > 0:
        r = hl.iloc[0]
        lines.append(f"Variant: {headline_id} (SNCA gene body)")
        lines.append(f"    hg38 position:       chr4:{int(r['hg38_pos']):,}")
        lines.append(f"    hap1_base / hap2_base: {r['hap1_base']} / {r['hap2_base']}")
        lines.append(f"    astrocyte AG quantile: {r['astrocyte_max_quantile']:.4f}")
        lines.append(f"    astrocyte AG raw:      {r['astrocyte_max_raw']:+.4f}")
        lines.append(f"    NPC / H9 / mono qtl:   "
                     f"{r['npc_max_quantile']:.3f} / {r['h9_max_quantile']:.3f} / {r['monocyte_max_quantile']:.3f}")
        lines.append(f"    priority_score:        {r['priority_score']:.4f}")
    else:
        lines.append(f"{headline_id} not found in SNCA gene-body variants with AG scores.")
    lines.append("")

    lines.append("═" * 78)
    text = "\n".join(lines)
    with open(SUMMARY, 'w') as f:
        f.write(text)
    print(text)


# ══════════════════════════════════════════════════════════════
# 5. Figure
# ══════════════════════════════════════════════════════════════
def make_figure(per_ct, astro_vs_others, snca_gb, snca_with_dir):
    plt.rcParams.update({'font.family': 'Arial', 'font.size': 9,
                         'figure.dpi': 200, 'savefig.dpi': 300,
                         'savefig.bbox': 'tight'})

    fig = plt.figure(figsize=(12, 9))
    gs = gridspec.GridSpec(2, 2, hspace=0.35, wspace=0.30,
                           left=0.08, right=0.96, top=0.93, bottom=0.08)

    # ── A: gene-level ASE across lineages ──
    axA = fig.add_subplot(gs[0, 0])
    order = ['iPSC', 'FPP', 'astrocyte', 'microglia']
    per_ct = per_ct.set_index('CellType').reindex(order).reset_index()
    x = np.arange(len(per_ct))
    w = 0.4
    axA.bar(x - w/2, per_ct['Hap1_ASE'], w, label='hap1', color='#3B7DDD')
    axA.bar(x + w/2, per_ct['Hap2_ASE'], w, label='hap2', color='#F19A3E')
    for i, r in per_ct.iterrows():
        annot = '***' if r['p_vs_0.5'] < 1e-3 else '**' if r['p_vs_0.5'] < 0.01 else \
                '*' if r['p_vs_0.5'] < 0.05 else 'ns'
        top = max(r['Hap1_ASE'], r['Hap2_ASE'])
        axA.text(i, top + max(per_ct[['Hap1_ASE','Hap2_ASE']].max())*0.03, annot,
                 ha='center', va='bottom', fontsize=8, color='#333')
    axA.set_xticks(x)
    axA.set_xticklabels(order, fontsize=9)
    axA.set_ylabel('ASE-informative reads')
    axA.set_title('A. SNCA gene-level ASE across lineages',
                  fontsize=11, fontweight='bold', loc='left')
    axA.legend(loc='upper left', fontsize=8)
    axA.spines['top'].set_visible(False); axA.spines['right'].set_visible(False)

    axA.text(0.98, 0.98,
             f"Astro vs others\nFisher p = {astro_vs_others['p_value']:.2g}\nOR = {astro_vs_others['odds_ratio']:.2f}",
             transform=axA.transAxes, ha='right', va='top', fontsize=8,
             bbox=dict(boxstyle='round,pad=0.3', facecolor='#FFF8DC', edgecolor='#888'))

    # ── B: AG volcano for SNCA gene body ──
    axB = fig.add_subplot(gs[0, 1])
    x = snca_gb['astrocyte_max_raw']
    y = snca_gb['astrocyte_max_quantile']
    axB.scatter(x, y, s=18, alpha=0.6, c='#888', edgecolor='none')
    # highlight quantile >= 0.99
    hi = snca_gb[snca_gb['astrocyte_max_quantile'] >= 0.99]
    axB.scatter(hi['astrocyte_max_raw'], hi['astrocyte_max_quantile'],
                s=40, c='#D62728', edgecolor='#333', linewidth=0.5)
    # label headline
    hl = snca_gb[snca_gb['variant_id'] == 'chr4_89836694_C_T']
    if len(hl) > 0:
        r = hl.iloc[0]
        axB.annotate('chr4:89,836,694 C>T\n(headline)',
                     xy=(r['astrocyte_max_raw'], r['astrocyte_max_quantile']),
                     xytext=(r['astrocyte_max_raw']-0.05, r['astrocyte_max_quantile']-0.05),
                     fontsize=7, color='#D62728',
                     arrowprops=dict(arrowstyle='->', color='#D62728', lw=0.7))
    axB.axhline(0.99, ls='--', color='#888', lw=0.5)
    axB.axvline(0,    ls='--', color='#888', lw=0.5)
    axB.set_xlabel('AlphaGenome astrocyte raw score  (+/- = up/down)')
    axB.set_ylabel('AlphaGenome astrocyte quantile score')
    axB.set_title('B. AG astrocyte scores — SNCA gene-body variants',
                  fontsize=11, fontweight='bold', loc='left')
    axB.spines['top'].set_visible(False); axB.spines['right'].set_visible(False)

    # ── C: Cell-type specificity for top candidate ──
    axC = fig.add_subplot(gs[1, 0])
    top1 = snca_gb.head(1).iloc[0] if len(snca_gb) else None
    if top1 is not None:
        ct_labels  = ['astrocyte', 'NPC', 'H9', 'monocyte']
        ct_qtl     = [top1['astrocyte_max_quantile'], top1['npc_max_quantile'],
                      top1['h9_max_quantile'],        top1['monocyte_max_quantile']]
        ct_colors  = ['#D62728', '#1F77B4', '#2CA02C', '#9467BD']
        axC.bar(range(4), ct_qtl, color=ct_colors, edgecolor='#333')
        axC.set_xticks(range(4)); axC.set_xticklabels(ct_labels)
        axC.set_ylim(0.85, 1.0)
        axC.axhline(0.99, ls='--', color='#888', lw=0.5)
        axC.set_ylabel('AG max quantile')
        axC.set_title(f'C. Cell-type specificity — top hit ({top1["variant_id"]})',
                      fontsize=11, fontweight='bold', loc='left')
        axC.spines['top'].set_visible(False); axC.spines['right'].set_visible(False)
        for i, v in enumerate(ct_qtl):
            axC.text(i, v + 0.001, f'{v:.4f}', ha='center', va='bottom', fontsize=7)

    # ── D: Direction concordance ──
    axD = fig.add_subplot(gs[1, 1])
    dd = snca_with_dir[snca_with_dir['concordant'].notna()]
    if len(dd) > 0:
        n_conc = int(dd['concordant'].sum())
        n_disc = len(dd) - n_conc
        axD.bar([0, 1], [n_conc, n_disc], color=['#2CA02C', '#D62728'], edgecolor='#333')
        axD.set_xticks([0, 1])
        axD.set_xticklabels(['AG direction agrees\nwith observed ASE',
                             'AG direction opposes\nobserved ASE'], fontsize=8)
        axD.set_ylabel('# variants')
        axD.set_title(f'D. Direction concordance ({n_conc}/{len(dd)} = {n_conc/len(dd):.0%} agree)',
                      fontsize=11, fontweight='bold', loc='left')
        for i, v in enumerate([n_conc, n_disc]):
            axD.text(i, v + 0.05, str(v), ha='center', va='bottom', fontsize=10, fontweight='bold')
        axD.spines['top'].set_visible(False); axD.spines['right'].set_visible(False)
    else:
        axD.text(0.5, 0.5, 'No variants with both AG-direction and observed reads',
                 transform=axD.transAxes, ha='center', va='center')
        axD.axis('off')

    fig.suptitle('SNCA astrocyte-specific ASE flip — proof of concept',
                 fontsize=13, fontweight='bold', y=0.98)

    for fmt in ('png', 'pdf'):
        fig.savefig(f'{FIG}.{fmt}', dpi=300, bbox_inches='tight')
    print(f'Saved: {FIG}.png / .pdf')


# ══════════════════════════════════════════════════════════════
# main
# ══════════════════════════════════════════════════════════════
def main():
    print('Running SNCA astrocyte flip POC...\n')
    per_ct, astro_vs_others = validate_astrocyte_flip()
    snca_gb, snca_all = snca_astrocyte_candidates()
    snca_with_dir = direction_concordance(snca_all)

    # Save top candidate table
    top10 = snca_gb.head(10)[['variant_id','hg38_pos','hap1_base','hap2_base',
                              'astrocyte_max_quantile','astrocyte_max_raw',
                              'npc_max_quantile','h9_max_quantile','monocyte_max_quantile',
                              'astro_hap1_reads','astro_hap2_reads',
                              'priority_score','in_gene_body']]
    top10.to_csv(TOP_CSV, index=False)
    print(f'Saved: {TOP_CSV}\n')

    write_report(per_ct, astro_vs_others, snca_gb, snca_with_dir)
    make_figure(per_ct, astro_vs_others, snca_gb, snca_with_dir)


if __name__ == '__main__':
    main()
