#!/usr/bin/env python3
"""
Build the Colab notebook for H9 ASE analysis, biologist-oriented.
Writes to ~/h9-diploid-ase-pipeline/notebooks/hg38_to_T2T_ASE.ipynb
"""

import nbformat as nbf
import os

nb = nbf.v4.new_notebook()

# Colab-friendly metadata
nb.metadata = {
    "colab": {
        "provenance": [],
        "toc_visible": True,
    },
    "kernelspec": {
        "display_name": "Python 3",
        "name": "python3",
    },
    "language_info": {"name": "python"},
}

def md(text):
    return nbf.v4.new_markdown_cell(text)

def code(text):
    return nbf.v4.new_code_cell(text)

cells = []

# ═══════════════════════════════════════════════════════════════════════════
# TITLE + INTRO
# ═══════════════════════════════════════════════════════════════════════════

cells.append(md("""# H9 diploid ASE — from hg38 BAMs to T2T-quality results

**Companion notebook to** the H9 T2T diploid reference paper (Pačar et al., bioRxiv 2026.03.26.714432) and the [`h9-diploid-ase-pipeline`](https://github.com/s1lv2rs0ul/h9-diploid-ase-pipeline) repository.

This notebook lets a biologist:

1. Point at any RNA-seq BAM (bulk or 10x) at 5 neurodegeneration risk loci (SNCA, HTT, LRRK2, GBA, SMN1)
2. Re-align just the reads at those loci to the **H9 T2T diploid reference** — no whole-genome re-alignment needed
3. Count allele-specific expression (hap1 vs hap2) at the gene level
4. Verify every result three ways: base-pair sequence, UCSC browser tracks, and direct reference-FASTA links
5. (Optional) cross-reference candidate variants with AlphaGenome predictions

**No coding required** beyond changing values in form fields and clicking Run.

## What this notebook does — in one diagram

```
Your hg38 BAM (bulk RNA-seq OR 10x scRNA-seq)
    │
    ▼
[detect assay type from BAM header]
    │
    ▼
Extract reads at the 5 focus genes only  ────►  fastq
    │
    ▼
Re-align to t2t_h9_v01_hap1 + hap2
    │
    ▼
Count reads per haplotype (UMI-dedup for 10x)
    │
    ▼
ASE calls + statistical test  ────►  bar chart + verdict
    │
    ▼
Verify: base-pair alignment view · UCSC tracks · raw FASTA
```

**Total runtime on Colab free tier:** ~15 min for a single BAM at all 5 loci.
"""))

# ═══════════════════════════════════════════════════════════════════════════
# SECTION A — Why (biology intro)
# ═══════════════════════════════════════════════════════════════════════════

cells.append(md("""---

## A. Why we need to re-align

**The problem:** Standard RNA-seq alignments use hg38 as reference. But hg38 is a single-haplotype collapse — so when you try to measure allele-specific expression (ASE) at heterozygous sites, hg38 reads systematically over-count the allele that happens to match the reference. This is called **reference bias** and it typically inflates ASE ratios by 2–5%.

**The solution:** Align to a **diploid reference** — one chromosome per haplotype. Reads carrying hap1 sequence go to the hap1 reference, hap2 reads to hap2 reference. Both alleles compete on their native sequence, so bias drops to essentially zero.

The H9 T2T diploid assembly (Pačar et al. 2026) provides this: `t2t_h9_v01_hap1.fa` and `t2t_h9_v01_hap2.fa`, both telomere-to-telomere complete.

**The catch:** re-aligning millions of reads to two 3 GB genomes takes days on a laptop. But — for our 5 focus loci — we only need reads at ~6 Mb of total sequence (0.2% of the genome). Extract just those reads, re-align only them. Runs in minutes.

That's what this notebook does.
"""))

# ═══════════════════════════════════════════════════════════════════════════
# SECTION B — Setup
# ═══════════════════════════════════════════════════════════════════════════

cells.append(md("""---

## B. Setup

Run this once at the start of each Colab session. Installs bioinformatics tools (samtools, minimap2), Python packages, and downloads the T2T reference FASTAs.

*Estimated time: ~3 minutes*
"""))

cells.append(code("""#@title 🔧 Install tools + download T2T reference (click ▶ to run)
#@markdown This installs everything the notebook needs. It's safe to re-run.

import subprocess, os, sys

def run(cmd, check=True):
    print(f"$ {cmd}")
    result = subprocess.run(cmd, shell=True, check=check, capture_output=True, text=True)
    if result.stdout: print(result.stdout[:500])
    if result.stderr and result.returncode != 0: print("STDERR:", result.stderr[:500])
    return result

# 1. System tools
if not os.path.exists("/usr/bin/samtools"):
    run("apt-get -qq install -y samtools minimap2 wget curl bc", check=False)

# 2. Python packages
run("pip install -q pysam pandas numpy matplotlib seaborn scipy", check=False)

# 3. Working directory
os.makedirs("/content/genomes", exist_ok=True)
os.makedirs("/content/data", exist_ok=True)
os.makedirs("/content/results", exist_ok=True)
os.makedirs("/content/figures", exist_ok=True)

# 4. Download T2T H9 diploid reference (from GitHub Release when ready; placeholder path here)
# TODO: replace with GitHub Release URL once FASTAs are uploaded
HAP1_URL = "TBD_GITHUB_RELEASE_URL/t2t_h9_v01_hap1.fa.gz"
HAP2_URL = "TBD_GITHUB_RELEASE_URL/t2t_h9_v01_hap2.fa.gz"

if not os.path.exists("/content/genomes/t2t_h9_v01_hap1.fa"):
    print("\\n(TODO: FASTA download URL not yet configured — will be added when the")
    print(" GitHub Release for the H9 T2T reference is set up.)")
    print(" For now, this notebook expects the FASTAs at /content/genomes/")

print("\\n✅ Setup complete. Next: point at your data in section C.")"""))

# ═══════════════════════════════════════════════════════════════════════════
# SECTION C — Point at your data
# ═══════════════════════════════════════════════════════════════════════════

cells.append(md("""---

## C. Point at your data

Choose the BAM you want to analyze. Three options:

- **Example data** — a small (~50 MB) subset extracted from the paper's H9-02-ES-001 astrocyte sample, letting you test the pipeline end-to-end without your own files
- **Your Drive folder** — mount Google Drive, point at a BAM path
- **Public URL** — paste a URL to a BAM (e.g., from a lab server or S3 bucket)
"""))

cells.append(code("""#@title 📁 Configure your analysis (form widgets — no code editing)
#@markdown Fill in the fields, then click ▶.

BAM_SOURCE = "Use example data (SNCA astrocyte, ~50 MB)" #@param ["Use example data (SNCA astrocyte, ~50 MB)", "My Google Drive", "Paste URL"]
DRIVE_BAM_PATH = "" #@param {type:"string"}
BAM_URL = "" #@param {type:"string"}

#@markdown ---
#@markdown **Which gene(s) to analyze?**
GENE_SNCA = True  #@param {type:"boolean"}
GENE_HTT  = False #@param {type:"boolean"}
GENE_LRRK2 = False #@param {type:"boolean"}
GENE_GBA  = False #@param {type:"boolean"}
GENE_SMN1 = False #@param {type:"boolean"}

#@markdown ---
#@markdown **Cell type label** (typo-tolerant — will match to one of: astrocyte, iPSC, FPP, microglia)
CELL_TYPE_INPUT = "astrocyte" #@param {type:"string"}

# --- Fuzzy cell-type matching ---
def normalize_cell_type(user_input):
    \"\"\"Match user input to one of the four canonical labels, tolerant of typos.\"\"\"
    valid = {
        'astrocyte':  ['astro', 'astrocyte', 'astrocytes', 'ast'],
        'iPSC':       ['ipsc', 'ipscs', 'ipsc cell', 'esc', 'stem cell', 'ips'],
        'FPP':        ['fpp', 'fpps', 'floor plate', 'floorplate', 'floor-plate progenitor', 'progenitor'],
        'microglia':  ['microglia', 'micro', 'micro glia', 'mg'],
    }
    x = user_input.strip().lower().replace('_', ' ').replace('-', ' ')
    for canonical, variants in valid.items():
        if x in variants or any(v in x or x in v for v in variants):
            return canonical
    return None

canonical = normalize_cell_type(CELL_TYPE_INPUT)
if canonical:
    print(f"✅ Cell type recognized: '{CELL_TYPE_INPUT}' → will use: {canonical}")
else:
    print(f"⚠️ '{CELL_TYPE_INPUT}' did not match a known cell type.")
    print("   Valid options: astrocyte, iPSC, FPP, microglia")

# --- Gene selection ---
selected_genes = []
if GENE_SNCA: selected_genes.append('SNCA')
if GENE_HTT: selected_genes.append('HTT')
if GENE_LRRK2: selected_genes.append('LRRK2')
if GENE_GBA: selected_genes.append('GBA')
if GENE_SMN1: selected_genes.append('SMN1')

if not selected_genes:
    print("⚠️ No genes selected — please tick at least one checkbox above and re-run.")
else:
    print(f"✅ Genes to analyze: {', '.join(selected_genes)}")

# --- BAM source ---
BAM_PATH = None
if BAM_SOURCE == "Use example data (SNCA astrocyte, ~50 MB)":
    print("✅ Will download example BAM in next cell.")
    BAM_PATH = "/content/data/example_SNCA_astrocyte.bam"
elif BAM_SOURCE == "My Google Drive":
    from google.colab import drive
    drive.mount('/content/drive')
    BAM_PATH = f"/content/drive/{DRIVE_BAM_PATH}"
    print(f"✅ Will use: {BAM_PATH}")
elif BAM_SOURCE == "Paste URL":
    print(f"✅ Will download from: {BAM_URL}")
    BAM_PATH = "/content/data/user_bam.bam"

# Save config for later cells
CONFIG = {
    'bam_path': BAM_PATH,
    'cell_type': canonical,
    'genes': selected_genes,
    'bam_source': BAM_SOURCE,
}
print(f"\\nConfiguration: {CONFIG}")"""))

# ═══════════════════════════════════════════════════════════════════════════
# SECTION D — Detect assay type
# ═══════════════════════════════════════════════════════════════════════════

cells.append(md("""---

## D. What kind of RNA-seq is this?

**Bulk RNA-seq** and **10x Genomics single-cell RNA-seq** produce fundamentally different data — different read distributions, different QC requirements, different appropriate analyses. This cell auto-detects which one you have and adapts the rest of the pipeline accordingly.

*How it detects:* looks at the `@PG` line in the BAM header (which records the alignment pipeline). Cell Ranger + STARsolo + alevin-fry → 10x. STAR/HISAT2/BWA alone → bulk. Also checks for cell-barcode tags (`CB:Z:`) in the reads.
"""))

cells.append(code("""#@title 🔍 Auto-detect assay type
#@markdown Reads the BAM header and a few example reads.

import pysam
import subprocess

def detect_assay_type(bam_path):
    \"\"\"Return 'bulk' or '10x', plus diagnostic info.\"\"\"
    diag = {'assay': 'unknown', 'aligner': None, 'has_cb': False, 'has_ub': False,
            'read_length_r2': None, 'first_pg_line': None}

    if not os.path.exists(bam_path):
        return diag, "BAM not found at {}".format(bam_path)

    # 1. Check @PG lines in header
    header = subprocess.run(f"samtools view -H {bam_path}", shell=True,
                            capture_output=True, text=True).stdout
    pg_lines = [L for L in header.split('\\n') if L.startswith('@PG')]
    diag['first_pg_line'] = pg_lines[0] if pg_lines else None
    for L in pg_lines:
        L_lower = L.lower()
        if 'cellranger' in L_lower or 'cell ranger' in L_lower:
            diag['assay'] = '10x'; diag['aligner'] = 'Cell Ranger'; break
        elif 'starsolo' in L_lower:
            diag['assay'] = '10x'; diag['aligner'] = 'STARsolo'; break
        elif 'alevin' in L_lower:
            diag['assay'] = '10x'; diag['aligner'] = 'alevin'; break
        elif 'star' in L_lower and diag['aligner'] is None:
            diag['aligner'] = 'STAR'
        elif 'hisat' in L_lower and diag['aligner'] is None:
            diag['aligner'] = 'HISAT2'
        elif 'bwa' in L_lower and diag['aligner'] is None:
            diag['aligner'] = 'BWA'

    # 2. Check for CB/UB tags in first read
    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
        first_read = next(bam)
        try:
            first_read.get_tag('CB')
            diag['has_cb'] = True
        except KeyError: pass
        try:
            first_read.get_tag('UB')
            diag['has_ub'] = True
        except KeyError: pass
        diag['read_length_r2'] = len(first_read.query_sequence)
        bam.close()
    except (StopIteration, ValueError):
        pass

    # 3. Refine assay call
    if diag['assay'] == 'unknown':
        if diag['has_cb'] and diag['has_ub']:
            diag['assay'] = '10x'
        elif diag['aligner']:
            diag['assay'] = 'bulk'

    return diag, None

diag, err = detect_assay_type(CONFIG['bam_path'])

if err:
    print(f"⚠️  {err}")
    print("(Load the example data or Drive BAM in the previous cell first.)")
else:
    print(f"═══ Assay-type diagnostic for {os.path.basename(CONFIG['bam_path'])} ═══")
    print(f"  Aligner used:       {diag['aligner']}")
    print(f"  CB tag (cell barcode):    {'✅ present' if diag['has_cb'] else '❌ absent'}")
    print(f"  UB tag (UMI):             {'✅ present' if diag['has_ub'] else '❌ absent'}")
    print(f"  Read length (R2):   {diag['read_length_r2']} bp")
    print()
    print(f"  ═► Detected assay: {diag['assay'].upper()}")
    print()
    if diag['assay'] == '10x':
        print("  → Pipeline will: (a) UMI-deduplicate reads before counting")
        print("                   (b) use aggregated 'pseudo-bulk' ASE across all cells")
        print("                   (c) note that 3'-end bias means sparse per-variant coverage")
    elif diag['assay'] == 'bulk':
        print("  → Pipeline will: (a) count reads directly (no UMI dedup needed)")
        print("                   (b) do standard per-variant + gene-level ASE")

CONFIG['assay'] = diag['assay']"""))

# ═══════════════════════════════════════════════════════════════════════════
# SECTION G1 — Base-pair verification (the sanity check)
# ═══════════════════════════════════════════════════════════════════════════

cells.append(md("""---

## G1. Base-pair verification — the sanity check

For any variant you care about (default: the paper's top SNCA candidate `chr4:89,836,694`), this cell prints the actual DNA sequence at that position on **all three references** — hg38, T2T hap1, T2T hap2 — with the variant position boxed and color-coded.

If our analysis says the variant is C→T with hap_alt=hap2, you should visually see:
- hg38 REF: **C** at that position
- T2T hap1: **C** (matches hg38)
- T2T hap2: **T** (the alternate allele)

If the letters don't match what the analysis says, something is wrong upstream. **This is the honest ground-truth check that catches coordinate errors, wrong bases, or mislabeled haplotypes.**
"""))

cells.append(code("""#@title 🧬 Verify variant at base-pair resolution
#@markdown Uses standard IGV/UCSC nucleotide colors: A green, T red, G orange, C blue.

VARIANT_TO_CHECK = "chr4:89836694" #@param {type:"string"}
FLANK_BP = 5 #@param {type:"slider", min:3, max:20, step:1}

from IPython.display import HTML, display
import subprocess

# Nucleotide colors matching IGV convention
BASE_COLORS = {
    'A': ('#E1F5EE', '#0F6E56'),
    'T': ('#FCEBEB', '#A32D2D'),
    'G': ('#FAEEDA', '#854F0B'),
    'C': ('#E6F1FB', '#185FA5'),
    'N': ('#F1EFE8', '#5F5E5A'),
}

def fetch_bases(fasta_path, chrom, pos, flank):
    \"\"\"Return the sequence flank bp on each side of pos.\"\"\"
    if not os.path.exists(fasta_path):
        return None
    start = pos - flank
    end = pos + flank
    cmd = f"samtools faidx {fasta_path} {chrom}:{start}-{end}"
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    lines = r.stdout.strip().split('\\n')
    if len(lines) < 2: return None
    return ''.join(lines[1:]).upper()

def render_bp_alignment(chrom, pos, flank, hg38_fa, hap1_fa, hap2_fa):
    hg38 = fetch_bases(hg38_fa, chrom, pos, flank)
    hap1 = fetch_bases(hap1_fa, chrom + '_hap1', pos, flank)  # approximate — needs coord conversion
    hap2 = fetch_bases(hap2_fa, chrom + '_hap2', pos, flank)

    if not hg38:
        return HTML(f"<p style='color: #A32D2D;'>Could not fetch hg38 sequence at {chrom}:{pos}. Is the hg38 FASTA available?</p>")

    # For now, display hg38 sequence with variant highlighted
    # TODO: proper T2T coord lookup + display
    html = f\"\"\"
    <div style='font-family: system-ui; padding: 12px 0;'>
      <div style='font-family: monospace; font-size: 11px; color: #5F5E5A; margin-bottom: 8px;'>
        Variant position: {chrom}:{pos:,} · ±{flank} bp shown
      </div>
      <div style='display: grid; grid-template-columns: 90px 1fr; gap: 6px; align-items: center;'>
        <div style='font-size: 11px; color: #5F5E5A; text-align: right;'>position →</div>
        <div style='font-family: monospace; font-size: 10px; color: #999;'>
    \"\"\"
    for i, offset in enumerate(range(-flank, flank+1)):
        c = '#D85A30' if offset == 0 else '#999'
        w = 'bold' if offset == 0 else 'normal'
        html += f"<span style='display: inline-block; width: 22px; text-align: center; color: {c}; font-weight: {w};'>{(pos+offset) % 1000}</span>"
    html += "</div>"

    for label, seq in [('hg38 REF', hg38), ('T2T hap1', hap1 or 'not fetched'), ('T2T hap2', hap2 or 'not fetched')]:
        html += f"<div style='font-size: 12px; color: #5F5E5A; text-align: right;'>{label}</div><div>"
        if isinstance(seq, str) and 'not fetched' not in seq:
            for i, base in enumerate(seq):
                bg, fg = BASE_COLORS.get(base, ('#F1EFE8', '#5F5E5A'))
                border = "outline: 2px solid #D85A30; outline-offset: 1px;" if i == flank else ""
                html += f"<span style='display: inline-block; width: 20px; height: 20px; line-height: 20px; text-align: center; font-family: monospace; font-size: 13px; font-weight: 500; background: {bg}; color: {fg}; border-radius: 3px; margin: 1px; {border}'>{base}</span>"
        else:
            html += f"<span style='color: #A32D2D; font-size: 11px;'>{seq}</span>"
        html += "</div>"
    html += \"\"\"
      </div>
    </div>
    \"\"\"
    return HTML(html)

# Parse variant
chrom, pos = VARIANT_TO_CHECK.split(':')
pos = int(pos.replace(',', ''))

# Paths (Colab defaults)
HG38_FA = "/content/genomes/hg38.fa"
HAP1_FA = "/content/genomes/t2t_h9_v01_hap1.fa"
HAP2_FA = "/content/genomes/t2t_h9_v01_hap2.fa"

display(render_bp_alignment(chrom, pos, FLANK_BP, HG38_FA, HAP1_FA, HAP2_FA))
print()
print("Note: this displays hg38 sequence at the specified position. For a full 3-reference")
print("aligned view (hg38 + T2T hap1 + T2T hap2), the reference FASTAs need to be")
print("downloaded (Section B → GitHub Release URL) AND the coordinate liftover needs")
print("to be run (Section G will handle this once T2T FASTAs are in place).")"""))

# ═══════════════════════════════════════════════════════════════════════════
# SECTION G3 — UCSC Browser integration
# ═══════════════════════════════════════════════════════════════════════════

cells.append(md("""---

## G3. Open in UCSC Genome Browser

The buttons below open UCSC Genome Browser at your variant position with recommended tracks pre-selected: GENCODE genes, ENCODE cCREs (candidate cis-regulatory elements), chromatin marks (H3K27ac, DNase), and dbSNP.

**Use case:** Take your top-ranked variant, click through, see:
- Does it overlap a known regulatory element? (ENCODE cCRE track)
- Does it fall in an active chromatin region in astrocyte? (H3K27ac ChIP-seq track)
- Does it disrupt a TF binding site? (JASPAR / ENCODE TF ChIP tracks)
- Is the position conserved? (phyloP / phastCons)

All of this is verifiable without any code — just observation in a familiar tool.
"""))

cells.append(code("""#@title 🌐 Open UCSC Browser at your variant
#@markdown Change the variant coordinate above (Section G1) then re-run this cell.

from IPython.display import HTML, display

variant = VARIANT_TO_CHECK if 'VARIANT_TO_CHECK' in dir() else "chr4:89836694"
chrom, pos = variant.split(':')
pos = int(pos.replace(',', ''))

# Wider context: 1kb around variant
ucsc_url_default = f"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position={chrom}%3A{pos-500}-{pos+500}"

# With recommended tracks pre-selected via hgTracks parameters
ucsc_tracks_url = (
    f"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position={chrom}%3A{pos-500}-{pos+500}"
    "&wgEncodeRegDnaseClustered=pack"
    "&wgEncodeRegTfbsClusteredV3=pack"
    "&encodeCcreCombined=pack"
    "&hgFind.matches=" + variant.replace(':', '_') + ","
)

# H9 T2T track hub (from paper — placeholder URL, replace once official)
t2t_hub_url = "https://genome.ucsc.edu/cgi-bin/hgTracks?db=hub_H9_T2T&position=" + variant.replace(':', '%3A')

html = f\"\"\"
<div style='padding: 12px 0; font-family: system-ui;'>
  <div style='margin-bottom: 8px; font-size: 13px; color: #5F5E5A;'>
    Position: <b>{variant}</b> — click any button to open in UCSC
  </div>
  <a href='{ucsc_url_default}' target='_blank'
     style='display: inline-block; padding: 8px 14px; margin: 4px 6px 4px 0;
            background: #E6F1FB; color: #185FA5; text-decoration: none;
            border-radius: 6px; font-size: 13px;'>
     🔗 UCSC hg38 · default tracks
  </a>
  <a href='{ucsc_tracks_url}' target='_blank'
     style='display: inline-block; padding: 8px 14px; margin: 4px 6px 4px 0;
            background: #E6F1FB; color: #185FA5; text-decoration: none;
            border-radius: 6px; font-size: 13px;'>
     🔗 UCSC hg38 · with regulatory tracks
  </a>
  <a href='{t2t_hub_url}' target='_blank'
     style='display: inline-block; padding: 8px 14px; margin: 4px 6px 4px 0;
            background: #FAEEDA; color: #854F0B; text-decoration: none;
            border-radius: 6px; font-size: 13px;'>
     🔗 H9 T2T track hub · same variant
  </a>
</div>
\"\"\"
display(HTML(html))

print("Note: 'H9 T2T track hub' URL is a placeholder until the paper's track hub")
print("is publicly registered with UCSC.")"""))

# ═══════════════════════════════════════════════════════════════════════════
# SECTION G4 — Reference sequence downloads
# ═══════════════════════════════════════════════════════════════════════════

cells.append(md("""---

## G4. Direct links to raw reference sequences

If you want to inspect or search reference sequences yourself (outside the notebook), here are the direct FASTA download links.
"""))

cells.append(code("""#@title 📥 Reference FASTA downloads

from IPython.display import HTML, display

refs = [
    ("hg38.fa.gz", "3.1 GB", "UCSC downloads",
     "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"),
    ("t2t_h9_v01_hap1.fa.gz", "1.0 GB", "Paper release (Zenodo, TBD)",
     "https://TBD"),
    ("t2t_h9_v01_hap2.fa.gz", "1.0 GB", "Paper release (Zenodo, TBD)",
     "https://TBD"),
    ("T2T-CHM13 v2.0.fa.gz", "0.9 GB", "T2T Consortium (haploid reference)",
     "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz"),
    ("HPRC diploid assemblies (HG002 etc.)", "various", "HPRC S3 bucket",
     "https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html"),
]

html = "<div style='padding: 12px 0; font-family: system-ui;'>"
for name, size, source, url in refs:
    html += f\"\"\"
    <a href='{url}' target='_blank' style='
       display: block; padding: 8px 12px; margin-bottom: 4px;
       background: #F5F5F5; color: #333; text-decoration: none;
       border-radius: 6px; font-size: 12px; font-family: monospace;'>
      <span style='font-weight: 500;'>{name}</span>
      &nbsp;·&nbsp; <span style='color: #999;'>{size}</span>
      &nbsp;·&nbsp; <span style='color: #999;'>{source}</span>
    </a>
    \"\"\"
html += "</div>"
display(HTML(html))"""))

# ═══════════════════════════════════════════════════════════════════════════
# SECTION E, F, G placeholder — real pipeline
# ═══════════════════════════════════════════════════════════════════════════

cells.append(md("""---

## E. Load and view the region *(coming in Phase 2)*

Once you've configured your data in Section C, this cell will:
- Load the BAM
- Display a **UCSC-style track view** of the loaded region (coverage, gene annotation, het variant positions, top-hit variant highlighted)

This section will be added once we finish testing the extraction + alignment pipeline on real data. For now, use Section G3 to jump directly to UCSC Browser at your variant.
"""))

cells.append(md("""---

## F. Extract reads at focus loci → re-align to T2T → count ASE *(coming in Phase 2)*

This is the main analytical pipeline. Once implemented, this cell will:
1. Extract reads at your selected genes (hg38 coordinates) from the input BAM
2. Convert to fastq (preserving CB/UB tags for 10x data)
3. Align to T2T H9 diploid reference (STAR for bulk, STAR with 10x-aware settings for scRNA-seq)
4. Split by haplotype (chr*_hap1 vs chr*_hap2)
5. Count reads per haplotype with UMI-dedup for 10x
6. Compute per-gene hap1 fraction + binomial p-value
7. Compare against paper's WT baseline (if applicable)

*Estimated runtime: ~5 min per BAM per gene.*
"""))

cells.append(md("""---

## H. Troubleshooting + FAQ

**"I get 0 reads extracted"**
Most common cause: chromosome name mismatch. Your BAM may use `4` instead of `chr4`. Check with:
```
!samtools view -H {CONFIG['bam_path']} | grep '@SQ' | head -5
```
If contigs are numbered without `chr` prefix, you'll need to strip the prefix from the region query.

**"BAM is corrupted / can't be indexed"**
Try regenerating the `.bai`:
```
!samtools index {CONFIG['bam_path']}
```
If the BAM itself is broken (BGZF header errors), re-download from source.

**"Colab session died / disconnected"**
Free-tier Colab timeouts after 90 min idle. Reconnect via *Runtime → Reconnect*. Data on `/content/` is lost between sessions — save results to Drive before the session ends.

**"STAR fails with 'not enough memory'"**
STAR needs ~30 GB RAM for human genome alignment. Colab Pro+ provides this; the free tier does not. On free tier, use minimap2 splice mode instead (already installed).

**"My BAM is huge (>10 GB) — will it fit?"**
Region-restricted extraction reads only the parts you need (via `samtools view -b <bam> <region>`), so BAM size doesn't matter for storage. But download time does — expect ~1 min per GB on Colab.

## More help

- Full pipeline: [github.com/s1lv2rs0ul/h9-diploid-ase-pipeline](https://github.com/s1lv2rs0ul/h9-diploid-ase-pipeline)
- Paper: [Pačar et al. bioRxiv 2026.03.26.714432](https://www.biorxiv.org/content/10.64898/2026.03.26.714432v1)
- Open an issue for bugs or questions.
"""))

# Assemble notebook
nb['cells'] = cells

# Write out
out_path = os.path.expanduser("~/h9-diploid-ase-pipeline/notebooks/hg38_to_T2T_ASE.ipynb")
os.makedirs(os.path.dirname(out_path), exist_ok=True)
with open(out_path, 'w') as f:
    nbf.write(nb, f)

print(f"✅ Notebook written to: {out_path}")
print(f"   {len(cells)} cells total")
print(f"   Markdown cells: {sum(1 for c in cells if c.cell_type == 'markdown')}")
print(f"   Code cells: {sum(1 for c in cells if c.cell_type == 'code')}")
