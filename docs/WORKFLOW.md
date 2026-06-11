# WORKFLOW — step-by-step run order

All paths assume the Gadi layout: working area `/scratch/cy94/rs4477`,
read-only inputs under `/g/data/xl04/rs4477`. Logs go to
`/scratch/cy94/rs4477/logs/`.

## Stage 0 — inputs already in place
- ONT modBAM URLs: 1000g-ont S3 bucket (public, https).
- hg38 reference: `/scratch/cy94/rs4477/reference/hg38.fa`
- 1000G chr22 panel:
  `/g/data/xl04/.../genotypes/chr22/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz`

## Stage 1 — methylation extraction
```bash
# 1.1 build the manifest of new samples (LOGIN NODE — needs internet)
python3 scripts/01_extraction/00_make_chr22_manifest.py \
  --out /scratch/cy94/rs4477/1kgp-cis-mqtl/config/chr22_manifest_95.tsv --n 95

# 1.2 stream chr22 from S3  (copyq; one looping job, re-runnable)
qsub scripts/01_extraction/run_stream_chr22_loop.pbs

# 1.3 modkit pileup once BAMs are present  (normal queue, 8 threads)
qsub scripts/01_extraction/run_modkit_chr22_loop.pbs
```
Output: `methylation/chr22/<SID>/<SID>.chr22.cov5.bedmethyl.gz` (+ raw + tabix).
Each row is a CpG; col 4 = mod code (`m`=5mC, `h`=5hmC), col 10 = coverage,
col 11 = percent modified.

## Stage 2 — phenotype matrix
```bash
qsub scripts/02_matrix/run_build_matrix.pbs
```
Reads all per-sample bedMethyl (5mC rows, coverage ≥5), keeps CpGs covered in
≥90% of samples, imputes the few gaps with the per-CpG mean, M-value transforms,
and writes a sorted/bgzipped/tabixed tensorQTL BED:
`matrix/chr22_100samples/methylation_Mval.bed.gz`.

## Stage 3 — genotypes
```bash
qsub scripts/03_genotypes/build_genotypes.pbs
```
Subsets the 1000G panel to the 100 samples (GM→NA keep-list), renames NA→GM so
IDs match the methylation columns, keeps biallelic SNVs, forces `chr22`, drops
duplicates, MAF≥0.05 →
`genotypes/chr22_100samples/chr22.100samples.{bed,bim,fam}`.
The log's `in BOTH geno+pheno: 100` line confirms no samples are dropped.

## Stage 4 — covariates
```bash
qsub scripts/04_covariates/run_covariates.pbs
```
LD-prunes genotypes → 5 genotype PCs (population structure); PCA of the
methylation matrix → 10 methylation PCs (hidden technical / platform structure).
Output: `covariates/chr22_100samples/covariates.tsv` (15 × 100, covariates as
rows).

## Stage 5 — cis-mQTL
```bash
# 5a permutation (GPU; one empirical p-value per CpG)
qsub scripts/05_qtl/run_perm_gpu.pbs

# 5b FDR + significant CpGs
python3 scripts/05_qtl/05b_fdr_significant.py \
  --permutation .../tensorqtl/chr22_100s/permutation/chr22.perm.txt.gz \
  --out-dir .../tensorqtl/chr22_100s/fdr --label chr22 --fdr 0.05

# 5c nominal mapping for FDR-significant CpGs only
python3 scripts/05_qtl/05c_run_nominal_significant.py \
  --plink-prefix .../genotypes/chr22_100samples/chr22.100samples \
  --phenotype-bed .../matrix/chr22_100samples/methylation_Mval.bed.gz \
  --covariates .../covariates/chr22_100samples/covariates.tsv \
  --significant-cpgs .../tensorqtl/chr22_100s/fdr/chr22.significant_cpgs.txt \
  --out-dir .../tensorqtl/chr22_100s/nominal --prefix chr22.sig
```

## Gadi gotchas learned the hard way
- Array jobs need `#PBS -r y`; single-element ranges (`-J 1-1`) are illegal.
- `copyq` allows only **1 CPU** and limits array size → use a single looping job.
- Only `copyq` / data-mover nodes have external network; everything else is
  offline. `bgzip`/`tabix` come from the modkit env, not the samtools module.
- The V100 needs a torch built with `sm_70` (use `torch==2.4.1+cu121`).
