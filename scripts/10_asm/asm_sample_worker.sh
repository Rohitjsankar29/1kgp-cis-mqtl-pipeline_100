#!/bin/bash
# asm_sample_worker.sh — full ASM chain for ONE sample, with cleanup.
# Idempotent: skips if the sample's asm_regions.tsv already exists.
# Deletes the big intermediates (raw calls, haplotagged BAM, filtered calls,
# phased VCF) after writing the small per-sample result — so 100 samples don't
# fill scratch.
#
# Usage: asm_sample_worker.sh SAMPLE MODBAM REGIONS PANEL REF ASMROOT [THREADS]
set -euo pipefail
S=$1; BAM=$2; REGIONS=$3; PANEL=$4; REF=$5; ASMROOT=$6; T=${7:-4}
SC=/scratch/cy94/rs4477/cov_v2
MODKIT=/scratch/cy94/rs4477/micromamba/envs/mqtl_modkit/bin/modkit
WHATSHAP=/scratch/cy94/rs4477/micromamba/envs/whatshap/bin/whatshap
BCFTOOLS=/scratch/cy94/rs4477/micromamba/envs/mqtl_modkit/bin/bcftools
SAMTOOLS=/scratch/cy94/rs4477/micromamba/envs/mqtl_modkit/bin/samtools
TABIX=/scratch/cy94/rs4477/micromamba/envs/mqtl_modkit/bin/tabix
PY=/scratch/cy94/rs4477/micromamba/envs/mqtl_tensorqtl/bin/python3

OUT="$ASMROOT/$S"; mkdir -p "$OUT"
RES="$OUT/$S.asm_regions.tsv"
if [ -s "$RES" ]; then echo "[$S] already done -> $RES (skip)"; exit 0; fi
echo "[$S] $(date +%T) start"

# reuse an existing haplotagged BAM + htlist if present (e.g. HG00355)
BAMH="$OUT/$S.haplotagged.bam"; HT="$OUT/$S.htlist.tsv"
if [ ! -s "$BAMH" ] || [ ! -s "$HT" ]; then
  "$BCFTOOLS" view -s "$S" -r chr22 -Oz -o "$OUT/$S.phased.vcf.gz" "$PANEL"
  "$TABIX" -f -p vcf "$OUT/$S.phased.vcf.gz"
  "$WHATSHAP" haplotag --ignore-read-groups --skip-missing-contigs --reference "$REF" \
    --output-haplotag-list "$HT" -o "$BAMH" "$OUT/$S.phased.vcf.gz" "$BAM"
  "$SAMTOOLS" index "$BAMH"
fi

echo "[$S] $(date +%T) extract calls"
"$MODKIT" extract calls --reference "$REF" -t "$T" "$BAMH" "$OUT/$S.calls.raw.tsv"

echo "[$S] $(date +%T) filter+join"
"$PY" "$SC/asm_join_calls.py" "$OUT/$S.calls.raw.tsv" "$HT" "$REGIONS" "$OUT/$S.calls.tsv"

echo "[$S] $(date +%T) ASM caller"
"$PY" "$SC/asm_readlevel.py" --calls "$OUT/$S.calls.tsv" --regions "$REGIONS" --out "$RES"

# cleanup big intermediates (keep result + htlist)
rm -f "$OUT/$S.calls.raw.tsv" "$OUT/$S.calls.tsv" \
      "$OUT/$S.phased.vcf.gz" "$OUT/$S.phased.vcf.gz.tbi"
# also drop the haplotagged BAM to save space (comment out to keep for re-use)
rm -f "$BAMH" "$BAMH.bai"
echo "[$S] $(date +%T) done -> $RES  (intermediates cleaned)"
