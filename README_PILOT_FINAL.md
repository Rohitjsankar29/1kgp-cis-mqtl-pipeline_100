# cis-mQTL prioritisation pipeline (1000 Genomes ONT)

A reproducible, CpG-centric **cis-methylation QTL** framework built on long-read
(Oxford Nanopore) data from the 1000 Genomes Project, developed on chromosome 22
as a pilot before scaling genome-wide and applying to the PROPHECY Aboriginal
cohort.

**Author:** Rohit Jaya Sankar — MSc, University of Western Australia
**Supervisor:** Dr Sam Buckberry
**Compute:** NCI Gadi (project `cy94`)

---

## Pipeline overview

```
ONT modBAM (S3)
   │  stream chr22 only (~2 GB/sample, not ~100 GB)         [01_extraction]
   ▼
modkit pileup  →  per-sample CpG bedMethyl (5mC + 5hmC)     [01_extraction]
   │
   ▼
methylation matrix  →  tensorQTL phenotype BED (M-values)   [02_matrix]
   │
1000G high-cov VCF  →  PLINK genotypes (chr22, MAF≥0.05)    [03_genotypes]
   │
   ▼
covariates = sex + platform + 5 genotype PCs + 5 residual meth PCs  [04_covariates]
   │
   ▼
tensorQTL  →  permutation → FDR → nominal-for-significant    [05_qtl]
   │
   ▼
downstream: SuSiE fine-mapping · SV integration · annotation · prioritisation  [06–09] ✅
```

## Repository layout

```
scripts/
  01_extraction/   make manifest, stream chr22 from S3, modkit methylation
  02_matrix/       collapse bedMethyl → tensorQTL phenotype BED
  03_genotypes/    subset 1000G panel → PLINK (with GM→NA ID mapping)
  04_covariates/   sex + platform + genotype PCs + residual methylation PCs
  05_qtl/          tensorQTL permutation / FDR / nominal / permuted-null + GPU jobs
  06_finemap/      SuSiE fine-mapping (PIP + 95% credible sets)
  07_sv_integration/  structural-variant association + LD tagging (novel)
  08_functional/   nearest gene/TSS, CpG islands, ENCODE cCREs
  09_prioritise/   multi-evidence weighted + elastic-net variant ranking
config/            manifests, sample sheets, keep/rename lists (generated)
env/               environment + tool setup notes
docs/              WORKFLOW · DOWNSTREAM · METHODS_AND_CODE_WALKTHROUGH · PILOT_RESULTS
```

See **docs/WORKFLOW.md** for the exact run order and **env/ENVIRONMENT.md** for
the toolchain.

## Status (chr22 pilot, n = 100)

| Stage | Output | Status |
|-------|--------|--------|
| Methylation extraction | 100 samples, chr22, 5mC + 5hmC | ✅ |
| Phenotype matrix | 601,535 CpGs × 100 (M-values, covered in ≥90% samples) | ✅ |
| Genotypes | 111,853 chr22 SNVs × 100 (MAF≥0.05) | ✅ |
| Covariates | 12 (sex + platform + 5 genotype PCs + 5 residual methylation PCs) | ✅ |
| cis-mQTL permutation | λ = 1.446; **6,995** significant CpGs at FDR < 0.05 | ✅ |
| Permuted-null calibration | λ_null = 1.104 — test calibrated, inflation is real signal | ✅ |
| Nominal (significant CpGs) | ~10.1 M cis variant–CpG pairs | ✅ |
| SuSiE fine-mapping | 6,797 / 6,995 CpGs (97%) with a 95% credible set; 80,711 variants | ✅ |
| SV integration (novel) | 1,054 CpG–SV pairs flagged across 885 CpGs (~13%) | ✅ |
| Functional annotation | all 6,995 CpGs (gene/TSS, CpG islands, ENCODE cCREs) | ✅ |
| Prioritisation | 80,711 ranked pairs; elastic-net R² = 0.17; convergence ρ = 0.63 | ✅ |

**Pilot complete and validated** — see **`docs/PILOT_RESULTS.md`** for the full
result breakdown and **`docs/METHODS_AND_CODE_WALKTHROUGH.md`** for the
stage-by-stage method explanation. Next: nperm = 10000 final pass, then
genome-wide on the full 1KGP cohort (unrelated subset), then PROPHECY.

## Key design decisions

- **chr22-only streaming.** BAMs are read remotely from S3 with `samtools view
  <url> chr22`, pulling ~2 GB per sample instead of the ~100 GB whole-genome
  BAM — a ~50× reduction in data movement.
- **5mC and 5hmC both retained.** modkit is run without restricting
  `--modified-bases`, so every sample yields 5mC (`m`) and, where the basecaller
  supports it, 5hmC (`h`). 5mC is the primary phenotype; 5hmC is a secondary
  analysis. (R9 basecallers do not separate the two — handled via covariates.)
- **Chromosome convention is `chr22` end-to-end.** Methylation and genotypes are
  both forced to `chr22`; a mismatch silently yields zero cis variants.
- **Sample IDs.** ONT samples use `GM#####` / `HG#####`; the 1000G genotype panel
  labels the same individuals `NA#####` / `HG#####`. Stage 03 maps `GM→NA` to
  subset, then renames back so genotype IDs match the methylation columns exactly.
- **GTEx-style 3-stage QTL.** permutation (1 p-value/CpG) → FDR → nominal mapping
  **only** for FDR-significant CpGs. This avoids the all-pairs memory blow-up of
  genome-wide nominal mapping in CpG-dense regions.
- **GPU.** tensorQTL runs on a Tesla V100. The V100 is compute capability 7.0
  (Volta), so the torch build must include `sm_70` — `torch==2.4.1+cu121` works;
  a `cu130` build does **not** (it dropped Volta and fails with "no kernel
  image available").

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
# 5. cis-mQTL permutation (GPU)
qsub scripts/05_qtl/run_perm_gpu.pbs
```

Paths in the PBS scripts are specific to the Gadi setup (`/scratch/cy94/rs4477`,
`/g/data/xl04/rs4477`); adjust the variables at the top of each for a different
environment.

## Downstream framework (stages 06–09)

After the cis-mQTL scan, the framework fine-maps, integrates structural variants
and functional annotation, and produces a ranked candidate-causal-variant list:

```bash
qsub scripts/05_qtl/run_nominal.pbs      # 05c nominal per-pair stats (sig CpGs)
qsub scripts/05_qtl/run_null_gpu.pbs     # 05d permuted-null calibration (λ≈1?)
qsub scripts/06_finemap/run_finemap.pbs  # SuSiE PIP + 95% credible sets
# 07 SV integration, 08 functional annotation, 09 prioritisation
python3 scripts/09_prioritise/prioritise_framework.py --susie ... --nominal ... --out ...
```

## Documentation & results

- **`docs/METHODS_AND_CODE_WALKTHROUGH.md`** — stage-by-stage annotated explanation
  of every script (purpose, rationale, key code, parameters) plus a viva concept
  glossary.
- **`docs/PILOT_RESULTS.md`** — the chr22 / n=100 pilot results: λ=1.446,
  6,995 significant CpGs, null λ=1.104, 97% credible-set resolution, SV
  integration, and prioritisation convergence (ρ=0.63).
- **`docs/WORKFLOW.md`**, **`docs/DOWNSTREAM.md`** — workflow and downstream notes.
