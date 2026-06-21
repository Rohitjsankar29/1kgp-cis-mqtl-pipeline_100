# Pilot results — chr22, n = 100 (1000 Genomes ONT)

This records the outcome of the complete pipeline run on the chromosome-22 /
100-sample pilot. It is a *test sample*: the framework is validated here before
being run across all of 1000 Genomes and then applied to the PROPHECY Aboriginal
cohort. All numbers below come from the scripts in `scripts/` run on NCI Gadi
(project `cy94`, V100 GPU).

## Cohort and inputs

| Input | Value |
|---|---|
| Samples | 100 (1000G ONT, chr22) |
| SNV genotypes | 111,853 variants × 100 |
| Structural variants | 1,462 chr22 SVs × 100 (from the 1KGP SNV_INDEL_SV panel) |
| CpG phenotypes (matrix) | 601,535 CpGs (M-values); 597,579 tested in cis |
| Covariates | 12 — sex + platform_R10 + 5 genotype PCs + 5 residual methylation PCs (residual df ≈ 86) |

## Stage 05 — cis-mQTL mapping

| Metric | Value |
|---|---|
| Permutation (nperm = 1000) genomic inflation | **λ = 1.446** |
| Significant CpGs at FDR < 0.05 | **6,995** (of 597,579 tested) |
| Nominal pass (significant CpGs) | ~10.1 M cis variant–CpG pairs |
| **Permuted-null calibration** (genotype labels shuffled, 19,879 CpGs) | **λ_null = 1.104** |

**Calibration verdict.** Shuffling the genotype sample labels collapses the
inflation from 1.446 toward 1 (λ_null = 1.10), so the bulk of the observed
inflation depends on the genuine genotype–phenotype correspondence — it is real
signal, not a test artefact. The small residual (λ_null > 1) is consistent with
correlated cis tests, the finite permutation count, and mean-imputation.
Separating true cis-mQTL signal from cryptic relatedness requires the
unrelated-subset / kinship analysis planned for the full-1KGP run.

## Stage 06 — SuSiE fine-mapping

| Metric | Value |
|---|---|
| CpGs with ≥ 1 95% credible set | **6,797 / 6,995 (97%)** |
| Credible-set variant rows | **80,711** |
| Max causal signals per CpG (L) | 10 |

## Stage 07 — structural-variant integration (novel)

| Metric | Value |
|---|---|
| Flagged CpG–SV pairs | **1,054** across **885 CpGs** (~13% of significant CpGs SV-implicated) |
| SV-type breakdown (flagged) | 506 DEL · 279 INS · 267 DUP · 2 INV |
| Evidence type | 935 direct-association only · 6 LD-tag only · 113 both |
| Standout hit | chr22:48,698,960 INS, r² ≈ 1.0 with lead SNV, p ≈ 7e-24 |

**Caveat.** At n = 100 an SV in r² ≈ 1 with a SNV is statistically
indistinguishable from it; these are *nominations* of candidate SV-driven mQTLs,
not proof the SV (rather than the linked SNV) is causal.

## Stage 08 — functional annotation

All 6,995 significant CpGs annotated with nearest gene/TSS + distance, CpG-island
overlap, and ENCODE cCRE overlap (GENCODE v44, UCSC cpgIslandExt, ENCODE SCREEN
GRCh38 cCREs).

## Stage 09 — prioritisation

| Metric | Value |
|---|---|
| Ranked variant–CpG pairs | **80,711** (credible-set variants) |
| Feature vector | pip, maf, prox_cpg, sv_flag, var_ccre, var_cpg_island, prox_tss, neg_log10_p, abs_effect, neg_log10_fdr |
| Elastic-net (predicting PIP) | α = 0.0019, l1_ratio = 0.1, **R² = 0.17** |
| Dominant learned coefficients | neg_log10_p = +0.194 · neg_log10_fdr = −0.176 · prox_cpg = +0.036 · prox_tss = +0.022 |
| **Convergence** (Spearman, weighted vs elastic-net) | **ρ = 0.629** |

**Interpretation.** PIP is driven mostly by association strength (per-pair p), with
the negative fdr coefficient reflecting that strongly-significant CpGs have broader
credible sets (PIP spreads across more variants). Functional/SV flags add little
*marginal* power to predict PIP at this scale — expected at n = 100 — but remain
central to biological interpretation. The two scores converge (ρ = 0.63) without
being redundant, which is the framework's robustness criterion in the absence of
causal ground truth.

## What the pilot supports

1. **Calibrated** — covariates + ancestry PCs + residual factors control
   confounding; the permuted null confirms the inflation isn't an artefact.
2. **Localised and integrated** — SuSiE credible sets + SV + functional +
   association evidence combine into one transparent ranking.
3. **Reproducible and scalable** — every stage scripted, every intermediate on
   disk, heavy reads filtered, ready for full 1KGP and PROPHECY.

It does **not** claim causality of any single variant: the output is a prioritised
shortlist for follow-up.

## Next steps

- Final permutation at nperm = 10000 (smoother tails) and 2–3 extra null seeds
  (bulletproof λ_null).
- Genome-wide run across all autosomes on the full 1KGP cohort, using the 2,504
  unrelated subset (or a kinship model) to separate signal from relatedness.
- Apply to PROPHECY (age/sex/clinical covariates available; ancestry/LD-transfer
  and Indigenous data-governance considerations to address with the supervisor).
