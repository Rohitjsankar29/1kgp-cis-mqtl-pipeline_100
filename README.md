# chr22 methylation extraction — 95 new samples

Reuses your proven prototype (`scripts/prototype/01_extract_chr_modkit.prototype.sh`)
which streams chr22 from S3 into modkit. Nothing about the methylation recipe
changes — this just builds the manifest and runs the prototype as an array.

## 0. Get the pack onto Gadi
scp this zip to scratch and unzip (same as before), or drop the two files into
your repo. The wrapper expects the prototype at:
  /scratch/cy94/rs4477/1kgp-cis-mqtl-pipeline/scripts/prototype/01_extract_chr_modkit.prototype.sh

## 1. Build the manifest (LOGIN NODE — needs internet)
```bash
/scratch/cy94/rs4477/micromamba/envs/mqtl_tensorqtl/bin/python3 \
  extraction/00_make_chr22_manifest.py \
  --out /scratch/cy94/rs4477/1kgp-cis-mqtl/config/chr22_manifest_95.tsv \
  --n 95
```
This lists the bucket, skips your 5 done samples, and writes 95 rows. Check the
R9/R10 split it prints — you want a mix (it becomes the basecaller covariate).

## 2. TEST ON ONE SAMPLE FIRST
```bash
mkdir -p /scratch/cy94/rs4477/logs/extract
qsub -J 1-1 extraction/run_extract_chr22.pbs
```
Watch it: `qstat -u rs4477`. When done, check the output:
```bash
SID=$(sed -n 2p /scratch/cy94/rs4477/1kgp-cis-mqtl/config/chr22_manifest_95.tsv | cut -f1)
ls -la /scratch/cy94/rs4477/1kgp-cis-mqtl/methylation/chr22/$SID/
zcat /scratch/cy94/rs4477/1kgp-cis-mqtl/methylation/chr22/$SID/$SID.chr22.cov5.bedmethyl.gz | head
cat  /scratch/cy94/rs4477/1kgp-cis-mqtl/methylation/chr22/$SID/$SID.chr22.modkit.log
```
Confirm: the bedMethyl has rows, chrom is `chr22`, and the CpG count is in the
same ballpark as your existing 5 samples. If the job errored at `command -v
bgzip/tabix`, load htslib too (add `module load htslib` to the wrapper). If copyq
rejected `ncpus=4`, drop it to 1 (and `--threads 1`).

## 3. Launch the rest
```bash
qsub -J 2-95 extraction/run_extract_chr22.pbs
```
copyq runs a few array tasks at a time, so all 95 trickle through over hours.
Re-running is safe — the prototype skips samples whose output already exists.

## Notes
- BAMs are never stored: ~2 GB of chr22 is streamed per sample, not the ~100 GB
  whole-genome BAM.
- Output per sample: methylation/chr22/<SID>/<SID>.chr22.cov5.bedmethyl.gz (+ .tbi)
- Your existing 5 samples live under /g/data/xl04/.../methylation/chr22/ in a
  different layout — we reconcile both when building the 100-sample matrix next.
