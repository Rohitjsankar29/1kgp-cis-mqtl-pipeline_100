#!/bin/bash
# asm_per_sample.sh — allele-specific methylation (ASM) for one ONT sample, chr22.
#   1. extract this sample's PHASED chr22 genotypes from the 1000G panel
#   2. whatshap haplotag the modBAM  -> reads tagged HP:i:1 / HP:i:2
#   3. modkit pileup --partition-tag HP -> one CpG bedMethyl per haplotype
#   4. asm_from_haplotypes.py (hap1 vs hap2 Fisher test) -> SAMPLE.asm.tsv
# Steps 1-3 run in the modkit env; step 4 uses the tensorqtl env python (scipy+pandas).
#
# Usage: bash asm_per_sample.sh SAMPLE MODBAM PANEL REF OUTDIR [THREADS]
set -euo pipefail
SAMPLE=$1; MODBAM=$2; PANEL=$3; REF=$4; OUT=$5; THREADS=${6:-8}
SC=/scratch/cy94/rs4477/cov_v2
PY=/scratch/cy94/rs4477/micromamba/envs/mqtl_tensorqtl/bin/python3   # has scipy+pandas
mkdir -p "$OUT"

echo "[1/4] $(date +%T) extract phased genotypes for $SAMPLE"
bcftools view -s "$SAMPLE" -r chr22 -Oz -o "$OUT/$SAMPLE.chr22.phased.vcf.gz" "$PANEL"
tabix -f -p vcf "$OUT/$SAMPLE.chr22.phased.vcf.gz"

echo "[2/4] $(date +%T) haplotag reads"
whatshap haplotag --ignore-read-groups --skip-missing-contigs --reference "$REF" \
  --output-haplotag-list "$OUT/$SAMPLE.htlist.tsv" \
  -o "$OUT/$SAMPLE.haplotagged.bam" \
  "$OUT/$SAMPLE.chr22.phased.vcf.gz" "$MODBAM"
samtools index "$OUT/$SAMPLE.haplotagged.bam"
echo "  reads per haplotype (want a balanced 1/2 split, not all 'none'):"
sed 1d "$OUT/$SAMPLE.htlist.tsv" | cut -f2 | sort | uniq -c

echo "[3/4] $(date +%T) haplotype-partitioned CpG pileup"
PU="$OUT/pileup_$SAMPLE"; mkdir -p "$PU"
modkit pileup --cpg --combine-strands --partition-tag HP --ref "$REF" \
  -t "$THREADS" "$OUT/$SAMPLE.haplotagged.bam" "$PU"
echo "  pileup files:"; ls -la "$PU"
for hp in "$PU"/*_1.bed "$PU"/*_2.bed; do
  [ -f "$hp" ] && { bgzip -f "$hp"; tabix -f -p bed "$hp.gz"; }
done

echo "[4/4] $(date +%T) call ASM (hap1 vs hap2)"
HP1=$(ls "$PU"/*_1.bed.gz 2>/dev/null | head -1)
HP2=$(ls "$PU"/*_2.bed.gz 2>/dev/null | head -1)
[ -z "$HP1" ] || [ -z "$HP2" ] && { echo "ERROR: haplotype beds not found in $PU"; exit 1; }
"$PY" "$SC/asm_from_haplotypes.py" --hp1 "$HP1" --hp2 "$HP2" --min-cov 5 \
  --out "$OUT/$SAMPLE.asm.tsv"
echo "done -> $OUT/$SAMPLE.asm.tsv"
