#!/usr/bin/env bash
# ============================================================================
# ASE Counting: All genes x All cell types x Both haplotypes
# Uses .ase.bam files (gene body, introns included) + diploid BAMs for totals
# ============================================================================

SAMTOOLS="/opt/homebrew/bin/samtools"
ASE_DIR="$HOME/reganalysis/ase"
OUT="$HOME/reganalysis/ase_v2/ase_all_genes_counts.tsv"

echo -e "Gene\tCellType\tHap1_ASE\tHap2_ASE\tDiploid_Total" > "$OUT"

CELLTYPES="iPSC FPP astrocyte microglia"

echo "===== ASE Read Counts (gene body, .ase.bam) ====="
echo ""

for CT in $CELLTYPES; do

  HAP1_BAM="$ASE_DIR/t2t_h9_v01_hap1.H9-02-ES-001_${CT}.ase.bam"
  HAP2_BAM="$ASE_DIR/t2t_h9_v01_hap2.H9-02-ES-001_${CT}.ase.bam"
  DIP_BAM="$ASE_DIR/t2t_h9_v01_diploid.H9-02-ES-001_${CT}.bam"

  echo "--- $CT ---"

  H1=$($SAMTOOLS view -c "$HAP1_BAM" chr4_hap1:92701524-92810571)
  H2=$($SAMTOOLS view -c "$HAP2_BAM" chr4_hap2:90611547-90749628)
  DT=$($SAMTOOLS view -c "$DIP_BAM" chr4_hap1:92701524-92810571 2>/dev/null || echo "NA")
  echo "SNCA       H1=$H1  H2=$H2  Dip=$DT"
  echo -e "SNCA\t$CT\t$H1\t$H2\t$DT" >> "$OUT"

  H1=$($SAMTOOLS view -c "$HAP1_BAM" chr4_hap1:3051762-3254349)
  H2=$($SAMTOOLS view -c "$HAP2_BAM" chr4_hap2:3035268-3237826)
  DT=$($SAMTOOLS view -c "$DIP_BAM" chr4_hap1:3051762-3254349 2>/dev/null || echo "NA")
  echo "HTT        H1=$H1  H2=$H2  Dip=$DT"
  echo -e "HTT\t$CT\t$H1\t$H2\t$DT" >> "$OUT"

  H1=$($SAMTOOLS view -c "$HAP1_BAM" chr12_hap1:40193975-40366656)
  H2=$($SAMTOOLS view -c "$HAP2_BAM" chr12_hap2:40889815-41062145)
  DT=$($SAMTOOLS view -c "$DIP_BAM" chr12_hap1:40193975-40366656 2>/dev/null || echo "NA")
  echo "LRRK2      H1=$H1  H2=$H2  Dip=$DT"
  echo -e "LRRK2\t$CT\t$H1\t$H2\t$DT" >> "$OUT"

  H1=$($SAMTOOLS view -c "$HAP1_BAM" chr1_hap1:156397616-156407863)
  H2=$($SAMTOOLS view -c "$HAP2_BAM" chr1_hap2:149294168-149304415)
  DT=$($SAMTOOLS view -c "$DIP_BAM" chr1_hap1:156397616-156407863 2>/dev/null || echo "NA")
  echo "GBA1       H1=$H1  H2=$H2  Dip=$DT"
  echo -e "GBA1\t$CT\t$H1\t$H2\t$DT" >> "$OUT"

  H1=$($SAMTOOLS view -c "$HAP1_BAM" chr5_hap1:75190828-75219735)
  H2=$($SAMTOOLS view -c "$HAP2_BAM" chr5_hap2:71457984-71486893)
  DT=$($SAMTOOLS view -c "$DIP_BAM" chr5_hap1:75190828-75219735 2>/dev/null || echo "NA")
  echo "SMN1       H1=$H1  H2=$H2  Dip=$DT"
  echo -e "SMN1\t$CT\t$H1\t$H2\t$DT" >> "$OUT"

  H1=$($SAMTOOLS view -c "$HAP1_BAM" chr5_hap1:74694993-74723871)
  echo "SMN2       H1=$H1  H2=0 (no hap2 copy)  Dip=NA"
  echo -e "SMN2\t$CT\t$H1\t0\tNA" >> "$OUT"

  echo ""
done

echo "===== Results saved to $OUT ====="
cat "$OUT"
