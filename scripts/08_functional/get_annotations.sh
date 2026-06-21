#!/bin/bash
# Stage 08a: download + subset chr22 annotation tracks. Run on a login/copyq node
# (needs internet). VERIFY the URLs/column layouts -- UCSC table schemas drift, and
# you may want a newer GENCODE or tissue-specific regulatory maps for PROPHECY.
# Needs gawk (for match(...) with a capture array).
set -euo pipefail
REF=/scratch/cy94/rs4477/reference/annotation; mkdir -p "$REF"; cd "$REF"

# GENCODE genes -> TSS bed (chr22)
wget -qO gencode.gtf.gz \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
zcat gencode.gtf.gz | gawk -F'\t' '$1=="chr22" && $3=="gene"{
  match($9,/gene_name "([^"]+)"/,g); tss=($7=="+")?$4:$5;
  print $1"\t"tss"\t"tss+1"\t"g[1]"\t.\t"$7 }' \
  | sort -k1,1 -k2,2n > gencode_chr22_tss.bed

# UCSC CpG islands (hg38)
wget -qO cpgIslandExt.txt.gz \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz
zcat cpgIslandExt.txt.gz | awk -F'\t' '$2=="chr22"{print $2"\t"$3"\t"$4"\tCpG_island"}' \
  | sort -k1,1 -k2,2n > cpg_islands_chr22.bed

# ENCODE cCREs (hg38)
wget -qO ccre.txt.gz \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/encodeCcreCombined.txt.gz
zcat ccre.txt.gz | awk -F'\t' '$2=="chr22"{print $2"\t"$3"\t"$4"\t"$6}' \
  | sort -k1,1 -k2,2n > ccre_chr22.bed

echo "annotation tracks ready in $REF:"
wc -l gencode_chr22_tss.bed cpg_islands_chr22.bed ccre_chr22.bed
