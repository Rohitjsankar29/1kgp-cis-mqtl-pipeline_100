# cis-mQTL prioritisation pipeline (1000 Genomes ONT)

A reproducible, CpG-centric cis-methylation QTL framework built on long-read (Oxford Nanopore) data from the 1000 Genomes Project, developed on chromosome 22 as a pilot before scaling genome-wide and applying to the PROPHECY Aboriginal cohort.

**Author:** Rohit Jaya Sankar — MSc, University of Western Australia · **Supervisor:** Dr Sam Buckberry · **Compute:** NCI Gadi (project cy94)

## Pipeline overview

```
ONT modBAM (S3)
   |  stream chr22 only (~2 GB/sample, not ~100 GB)          [01_extraction]
   v
modkit pileup  ->  per-sample CpG bedMethyl (5mC)            [01_extraction]
   |
   v
methylation matrix  ->  tensorQTL phenotype BED (M-values)   [02_matrix]
   |
1000G high-cov VCF  ->  PLINK genotypes (chr22, MAF>=0.05)   [03_genotypes]
   |
   v
covariates = sex + platform + 5 genotype PCs + 5 residual meth PCs  [04_covariates]
   |
   v
CpG-SNP filter  ->  drop CpGs whose C/G overlaps an in-sample SNP  [04b_filter]  <- NEW
   |
   v
tensorQTL  ->  permutation -> FDR -> nominal-for-significant  [05_qtl]
   |
   v
downstream: SuSiE fine-mapping - SV integration - annotation - prioritisation  [06-09]
   |
   v
ASM corroboration  ->  allele-specific methylation (Atlas method)  [10_asm]  <- NEW
```

## Repository layout

```
scripts/
  01_extraction/   make manifest, stream chr22 from S3, modkit methylation
  02_matrix/       collapse bedMethyl -> tensorQTL phenotype BED
  03_genotypes/    subset 1000G panel -> PLINK (with GM->NA ID mapping)
  04_covariates/   sex + platform + genotype PCs + residual methylation PCs
  05_tensorqtl/    04b_filter_cpg_snps.py (CpG-SNP artefact filter) + tensorQTL jobs
  05_qtl/          tensorQTL permutation / FDR / nominal / permuted-null + GPU jobs
  06_finemap/      SuSiE fine-mapping (PIP + 95% credible sets)
  07_sv_integration/  structural-variant association + LD tagging (novel)
  08_functional/   nearest gene/TSS, CpG islands, ENCODE cCREs
  09_prioritise/   multi-evidence weighted + elastic-net variant ranking
  10_asm/          allele-specific methylation module (Atlas method)  <- NEW
  plots/           make_figures.py - reproducible figure generation
config/            manifests, sample sheets, keep/rename lists (generated)
env/               environment + tool setup notes
docs/              WORKFLOW - DOWNSTREAM - METHODS_AND_CODE_WALKTHROUGH - RESULTS
results/           prioritised shortlists, ASM classification, readable tables
figures/           generated figures
```

See `docs/WORKFLOW.md` for the exact run order and `env/ENVIRONMENT.md` for the toolchain.

## Status (chr22 pilot, n = 100 - post CpG-SNP filter)

| Stage | Output | Status |
|---|---|---|
| Methylation extraction | 100 samples, chr22, 5mC | done |
| Phenotype matrix | 601,535 CpGs x 100 (M-values, >=90% samples) | done |
| Genotypes | 111,853 chr22 SNVs x 100 (MAF>=0.05) | done |
| Covariates | 12 (sex + platform_R9 + 5 genotype PCs + 5 residual methylation PCs) | done |
| **CpG-SNP filter (NEW)** | **17,685 CpGs removed (2.9%); artefacts removed pre-scan** | done |
| cis-mQTL permutation | lambda = 1.45; **6,539** significant CpGs at FDR < 0.05 | done |
| Permuted-null calibration | lambda_null = 1.10 - test calibrated, inflation is real signal | done |
| Nominal (significant CpGs) | per-pair stats for the 6,539 significant CpGs | done |
| SuSiE fine-mapping | credible sets across 6,356 CpGs; **75,529** candidate variants | done |
| SV integration (novel) | CpG-SV pairs flagged (~13% of hits) | done |
| Functional annotation | all 6,539 CpGs (gene/TSS, CpG islands, ENCODE cCREs) | done |
| Prioritisation | 75,529 ranked pairs; elastic-net R2 = 0.086; convergence rho = 0.559 | done |
| **ASM corroboration (NEW)** | **67 SD-ASM, 7 imprint (56 samples); 67/67 sign-test concordance** | done |

Pilot complete and validated - see `docs/RESULTS.md` for the full breakdown and `docs/METHODS_AND_CODE_WALKTHROUGH.md` for the stage-by-stage method explanation. Next: R10-only scale-up (single chemistry), then genome-wide on the full 1KGP cohort (unrelated subset), then PROPHECY.

## Headline results

**Top corroborated hit:** TBC1D22A (`chr22:46773373:C:T`) - fine-mapping uncertain (PIP 0.49) but allele-specific methylation consistent in 27/27 heterozygous carriers. ASM resolves what the statistics could not.

- **Sign test:** among the 67 SD-ASM pairs, ASM and mQTL effect directions agree 67/67 (p = 1.4e-20), oriented independently of the association - positive evidence the mQTLs are genuinely cis-causal.
- **Complementarity:** SD-ASM variants have median PIP 0.51 (where fine-mapping is uncertain) vs 0.998 for the top-50 by score; ASM corroborates the local, low-PIP variants fine-mapping cannot resolve (median SNP-CpG distance 187 bp).
- **Fine-mapping sanity check:** high-PIP is not an LD-sparsity artefact - PIP decreases with SNP-CpG distance (Spearman rho = -0.24).

## Key design decisions

**chr22-only streaming.** BAMs are read remotely from S3 with `samtools view <url> chr22`, pulling ~2 GB per sample instead of the ~100 GB whole-genome BAM - a ~50x reduction in data movement.

**5mC as the phenotype.** modkit is run with `--modified-bases 5mC` so the phenotype is CpG 5-methylation. (R9 basecallers compress 5mC dynamic range vs R10; the R9/R10 batch effect is controlled via the `platform_R9` covariate.)

**CpG-SNP filter before the scan (NEW).** A variant on the C or G of a CpG mechanically forces the methylation readout (genotype -> apparent methylation) - an artefact, not regulation. These were dominating the top of the ranking. Stage 04b removes any CpG whose C/G overlaps an in-sample variant (filtering against the full panel over-filters at 32%; in-sample gives the correct 2.9%). The former top-ranked variant (a CpG-SNP artefact) drops out.

**Chromosome convention is chr22 end-to-end.** Methylation and genotypes are both forced to chr22; a mismatch silently yields zero cis variants.

**Sample IDs.** ONT samples use GM##### / HG#####; the 1000G genotype panel labels the same individuals NA##### / HG#####. Stage 03 maps GM->NA to subset, then renames back so genotype IDs match the methylation columns exactly.

**In-sample LD.** SuSiE computes LD from the study genotypes - no external reference panel. This avoids out-of-sample mismatch false positives and is the only correct choice for transfer to ancestrally underrepresented cohorts (PROPHECY).

**GTEx-style 3-stage QTL.** permutation (1 p-value/CpG) -> FDR -> nominal mapping only for FDR-significant CpGs. This avoids the all-pairs memory blow-up of genome-wide nominal mapping in CpG-dense regions.

**ASM is within-individual (NEW).** The ASM layer (Atlas method) compares haplotypes within heterozygous individuals - a within-molecule readout orthogonal to the between-individual mQTL, and immune to the LD-transfer problem. Oriented by phased genotype only (never the mQTL slope), so the sign-test concordance is non-circular.

**GPU.** tensorQTL runs on a Tesla V100 (compute capability 7.0, Volta), so the torch build must include sm_70 - `torch==2.4.1+cu121` works; a cu130 build drops Volta and fails with "no kernel image available".

## Quick start

```bash
# 1. methylation  (copyq stream, then normal-queue modkit)
qsub scripts/01_extraction/run_stream_chr22_loop.pbs
qsub scripts/01_extraction/run_modkit_chr22_loop.pbs
# 2. phenotype matrix
qsub scripts/02_matrix/run_build_matrix.pbs
# 3. genotypes
qsub scripts/03_genotypes/build_genotypes.pbs
# 4. covariates
qsub scripts/04_covariates/run_covariates.pbs
# 4b. CpG-SNP filter (NEW - before the scan)
qsub scripts/05_tensorqtl/run_filter_cpgsnps.pbs
# 5. cis-mQTL permutation (GPU)
qsub scripts/05_qtl/run_perm_gpu.pbs
```

Paths in the PBS scripts are specific to the Gadi setup (`/scratch/cy94/rs4477`, `/g/data/xl04/rs4477`); adjust the variables at the top of each for a different environment.

## Downstream framework (stages 06-10)

```bash
qsub scripts/05_qtl/run_nominal.pbs      # 05c nominal per-pair stats (sig CpGs)
qsub scripts/05_qtl/run_null_gpu.pbs     # 05d permuted-null calibration
qsub scripts/06_finemap/run_finemap.pbs  # SuSiE PIP + 95% credible sets
# 07 SV integration, 08 functional annotation, 09 prioritisation (run as documented)
# 10 ASM corroboration (Atlas method, 56 samples)
qsub scripts/10_asm/run_asm_serial.pbs
python3 scripts/10_asm/combine_asm_oriented.py --asm-glob ... --out ...
python3 scripts/10_asm/join_asm_to_enriched.py --enriched ... --asm ... --out ...
```

## Figures

Regenerate all figures from the result tables (numbers computed live, nothing hardcoded):

```bash
python3 scripts/plots/make_figures.py \
  --enriched results/chr22.prioritised_variant_level.enriched_asm.txt.gz \
  --nominal  <downstream>/nominal/chr22.sig.cis_qtl_pairs.chr22.parquet \
  --outdir   figures/
```

| Figure | Shows |
|---|---|
| `figures/asm_composition_signtest.png` | ASM classification + 67/67 sign test |
| `figures/asm_distance.png` | ASM tests short-range pairs (187 bp median) |
| `figures/pip_vs_distance.png` | PIP falls with distance (no LD artefact) |
| `figures/asm_finemap_complementarity.png` | ASM vs fine-mapping complementarity |

## Documentation & results

- `docs/METHODS_AND_CODE_WALKTHROUGH.md` - stage-by-stage annotated explanation of every script (purpose, rationale, key code, parameters) plus a viva concept glossary.
- `docs/RESULTS.md` - the chr22 / n=100 post-filter pilot results: lambda=1.45, 6,539 significant CpGs, null lambda=1.10, credible-set resolution, SV integration, prioritisation convergence (rho=0.559), and the ASM validation (67 SD-ASM, 67/67 sign test).
- `docs/WORKFLOW.md`, `docs/DOWNSTREAM.md` - workflow and downstream notes.

## Honest scope

At n=100 the framework produces a **prioritised shortlist with within-molecule ASM corroboration**, not proof of causality for any single variant. ASM corroborates only the short-range subset (SNP and CpG within read span). Full causal claims await the larger cohort. LD is computed in-sample (no external reference panel) - the correct choice for transfer to ancestrally underrepresented cohorts (PROPHECY).
