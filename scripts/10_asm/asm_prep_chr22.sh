#!/bin/bash
# asm_prep_chr22.sh — turn one sample's haplotagged BAM + the significant-CpG list
# into the inputs asm_readlevel.py needs, then run it.
#
#   1. modkit extract calls  -> per-read CpG methylation calls
#   2. join whatshap htlist  -> add hp (1/2) per read (H1->1, H2->2, drop none)
#   3. build merged-block regions from the significant CpGs
#   4. run asm_readlevel.py  -> region-level ASM (Atlas method)
#
# Usage: bash asm_prep_chr22.sh SAMPLE HAPLOTAGGED_BAM HTLIST SIGCPGS REF OUTDIR [MERGE_BP]
set -euo pipefail
SAMPLE=$1; BAM=$2; HTLIST=$3; SIGCPGS=$4; REF=$5; OUT=$6; MERGE=${7:-1000}
SC=/scratch/cy94/rs4477/cov_v2
MODKIT=/scratch/cy94/rs4477/micromamba/envs/mqtl_modkit/bin/modkit
PY=/scratch/cy94/rs4477/micromamba/envs/mqtl_tensorqtl/bin/python3
mkdir -p "$OUT"

echo "[1/4] $(date +%T) modkit extract calls"
"$MODKIT" extract calls --reference "$REF" -t 8 "$BAM" "$OUT/$SAMPLE.calls.raw.tsv"

echo "[2/4] $(date +%T) join haplotype tags"
"$PY" - "$OUT/$SAMPLE.calls.raw.tsv" "$HTLIST" "$OUT/$SAMPLE.calls.tsv" << 'PY'
import sys, pandas as pd
raw, htlist, out = sys.argv[1], sys.argv[2], sys.argv[3]
# haplotype map: readname -> 1/2  (H1->1, H2->2, drop none)
ht = pd.read_csv(htlist, sep="\t")
ht.columns = [c.lstrip("#") for c in ht.columns]
hpmap = {"H1": 1, "H2": 2, "1": 1, "2": 2}
ht["hp"] = ht["haplotype"].astype(str).map(hpmap)
ht = ht.dropna(subset=["hp"])[["readname", "hp"]]
ht["hp"] = ht["hp"].astype(int)

c = pd.read_csv(raw, sep="\t")
# modkit extract columns: read_id, ..., chrom, ref_position, ..., call_code, ...
rid = "read_id" if "read_id" in c.columns else c.columns[0]
c = c.rename(columns={rid: "readname"})
keep = {}
for want, opts in {"chrom": ["chrom", "chromosome"],
                   "ref_position": ["ref_position", "ref_pos", "position"],
                   "call": ["call_code", "call", "mod_code"]}.items():
    for o in opts:
        if o in c.columns:
            keep[o] = want; break
c = c.rename(columns=keep)[["readname", "chrom", "ref_position", "call"]]
c = c.merge(ht, on="readname", how="inner")   # keep only haplotyped reads
c = c.rename(columns={"readname": "read_id"})
c.to_csv(out, sep="\t", index=False)
print(f"  calls with haplotype: {len(c):,}  (hp1={int((c.hp==1).sum()):,}  hp2={int((c.hp==2).sum()):,})")
PY

echo "[3/4] $(date +%T) build merged-block regions (merge within ${MERGE}bp)"
"$PY" - "$SIGCPGS" "$OUT/regions.bed" "$MERGE" << 'PY'
import sys, pandas as pd
sig, out, merge = sys.argv[1], sys.argv[2], int(sys.argv[3])
pos = sorted(int(l.strip().split("_")[-1]) for l in open(sig) if l.strip())
blocks, s, e = [], pos[0], pos[0]
for p in pos[1:]:
    if p - e <= merge:
        e = p
    else:
        blocks.append((s, e)); s = e = p
blocks.append((s, e))
with open(out, "w") as f:
    for i, (a, b) in enumerate(blocks):
        f.write(f"chr22\t{a-250}\t{b+250}\tblock{i}_{a}_{b}\n")   # pad 250bp each side
print(f"  {len(pos)} sig CpGs -> {len(blocks)} merged regions")
PY

echo "[4/4] $(date +%T) run ASM caller"
"$PY" "$SC/asm_readlevel.py" --calls "$OUT/$SAMPLE.calls.tsv" \
  --regions "$OUT/regions.bed" --out "$OUT/$SAMPLE.asm_regions.tsv"
echo "done -> $OUT/$SAMPLE.asm_regions.tsv"
