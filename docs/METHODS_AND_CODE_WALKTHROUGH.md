# cis-mQTL Variant-Prioritisation Framework — Annotated Code & Methods Walkthrough

A stage-by-stage explanation of every script in the pipeline: what it does, the
scientific rationale, its inputs and outputs, the key code logic, and the
parameter choices you should be able to defend. Written for the chr22 / n=100
pilot, but the logic is identical at full 1KGP and PROPHECY scale.

---

## The pipeline at a glance

```
ONT modBAM (per sample)
   │  01  extraction        → per-sample methylation calls (modkit bedMethyl)
   ▼
CpG × sample matrix          02  build_matrix.py        (phenotypes)
SNV genotypes (plink)        03  build_genotypes.pbs    (genotypes)
covariates                   04  metadata + covariates  (confounders)
   │
   ▼  05  cis-mQTL mapping (tensorQTL on GPU)
       05a permutation   → one empirical p per CpG
       05b FDR           → significant CpGs
       05c nominal       → per-pair effect sizes / p-values
       05d permuted null → calibration check (λ≈1?)
   │
   ▼  DOWNSTREAM (the framework)
       06  SuSiE fine-mapping     → PIP + 95% credible sets
       07  SV integration (novel) → does a structural variant explain the signal?
       08  functional annotation  → genes / cCREs / CpG islands
       09  prioritisation         → ranked candidate causal variants
```

Each stage produces a file the next stage consumes. Nothing is hidden in state;
every intermediate is on disk and inspectable.

---

## The core idea (one paragraph to anchor everything)

DNA methylation at a CpG site is the **phenotype**. Nearby genetic variants are
the **predictors**. A *cis-methylation QTL* (cis-mQTL) is a variant whose genotype
is statistically associated with methylation level at a CpG within ±1 Mb. Finding
the association is the easy part; the hard part — and the contribution of this
framework — is going from "this CpG has a signal somewhere in a 2 Mb window of
correlated variants" to a **ranked shortlist of the variants most likely to be
causal**, integrating fine-mapping, structural variants, and functional
annotation. Because there is no causal ground truth, the framework's claim is
*prioritisation and method convergence*, not proof of causality.

---

# Stage 01 — Extraction: `00_make_chr22_manifest.py`

**Purpose.** Build the work list of ONT samples to process: which BAM to stream,
from where, and what modification calls it carries.

**Why it exists.** The 1000 Genomes ONT data lives in a public S3 bucket with
hundreds of BAMs. You need a reproducible, deterministic selection of samples
(excluding ones already done) plus per-sample metadata that downstream steps key
off — most importantly the **pore chemistry (R9 vs R10)**, which is the dominant
technical batch effect and later becomes a covariate.

**Key logic.**
- `list_s3()` pages through the bucket via the S3 REST API (`list-type=2`,
  continuation tokens) using only the Python standard library — so it runs on a
  login node with no extra packages.
- `parse_name()` regex-parses each BAM filename to recover `sample_id`, pore
  (`R10` if `-R10-` in the name else `R9`), modifications (`5mC+5hmC` if `5hmC`
  appears, else `5mC`), and basecaller.
- Already-processed samples (`DONE_DEFAULT`) are skipped so re-runs never collide.

**Output.** A TSV manifest (`sample_id`, `bam_url`, pore, modifications, …). The
extraction PBS jobs (`run_stream_chr22_loop.pbs`, `run_modkit_chr22_loop.pbs`)
read it to stream each chr22 BAM and run `modkit pileup`, producing one
**bedMethyl** file per sample with per-CpG methylation counts.

**Defend it.** "Sample selection is scripted and deterministic, with technical
metadata (pore chemistry) captured at extraction time so it can be modelled as a
covariate rather than silently confounding the association."

---

# Stage 02 — Phenotypes: `build_matrix.py`

**Purpose.** Collapse the per-sample bedMethyl files into a single CpG × sample
methylation matrix in the BED format tensorQTL expects.

**Why it exists.** QTL mapping needs an aligned phenotype matrix: the same CpG
sites across all samples, with missingness handled. bedMethyl is per-sample and
ragged (different CpGs covered per sample), so it must be intersected and cleaned.

**Key logic.**
- `load_sample()` reads modkit's 5mC rows (`mod_code == 'm'`) and keeps only CpGs
  with `valid_coverage ≥ min_cov` (default **5×**). Methylation is taken as
  `percent_modified / 100` → a **beta value** in [0,1].
- The per-sample series are concatenated into a matrix; CpGs are kept only if
  covered in **≥ 90 %** of samples (`--min-frac 0.90`). This is the
  quality/coverage filter — it removes CpGs too sparsely observed to test.
- Remaining gaps are imputed with the **per-CpG (row) mean** — a minimal,
  unbiased fill so tensorQTL sees a complete matrix.
- Optional **M-value transform**: `log2((β+ε)/(1−β+ε))`. Beta values are bounded
  and heteroscedastic near 0 and 1; the logit (M-value) is closer to homoscedastic
  and better behaved for linear models.
- Writes a sorted BED: `#chr start end phenotype_id <samples…>`, with
  `phenotype_id = chr22_<pos>`.

**Parameters to defend.** `min-cov 5` (standard ONT methylation floor),
`min-frac 0.90` (CpG must be broadly observed), mean-imputation (minimal and only
for the ≤10 % residual gaps), M-value (variance stabilisation).

---

# Stage 03 — Genotypes: `build_genotypes.pbs`

**Purpose.** Produce the plink genotype set (`.bed/.bim/.fam`) for the same
samples — chr22 SNVs for the 100 samples.

**Why it exists.** tensorQTL reads genotypes through a plink reader. This step
subsets the 1KGP chr22 callset to your cohort and to SNVs, yielding
111,853 variants × 100 samples. (SVs are handled separately in Stage 07, because
they need different parsing.)

**Defend it.** Genotypes and phenotypes must be on an identical, identically
ordered sample set; the `.fam` file is the canonical sample list the metadata and
SV steps derive from.

---

# Stage 04 — Covariates: `build_sample_metadata.py` + `build_covariates.py`

This is the most important methodological stage for controlling false positives,
and the one your supervisor specifically asked to strengthen.

## `build_sample_metadata.py`
**Purpose.** Assemble per-sample known covariates: **sex** and **platform**.

- **Sex** comes from the 1000G panel file, not the genotypes — chr22-only data
  has no X chromosome, so sex can't be inferred from the calls. The code handles
  the `GM#####` (Coriell) ↔ `NA#####` (1000G) ID aliasing.
- **Platform** (R9/R10) comes from the extraction manifest — the dominant batch.
- **Age** is deliberately omitted: 1000G is de-identified and doesn't release age.
  The script explicitly flags this and leaves an `age` slot for PROPHECY, where
  age *is* available. (This is why the covariate builder already accepts age — the
  pipeline is forward-compatible with the real cohort.)

## `build_covariates.py` (v2 — the model you're running)
**Purpose.** Build the covariate matrix that gets regressed out of methylation
before testing for genetic association.

The covariates are three layers:
1. **Known covariates** — sex (0/1), platform dummy (R10 vs R9 reference), age if
   present.
2. **Genotype PCs** — the top 5 principal components from `plink2 --pca`. These
   capture **ancestry / population structure**, the classic confounder in any
   genetic association study (allele frequencies and methylation both vary by
   ancestry, producing spurious associations if uncontrolled).
3. **Residual methylation PCs** — and this is the subtle, defensible part. Rather
   than taking raw PCs of the methylation matrix (which would re-capture sex,
   platform, and ancestry you've *already* modelled, double-counting them), the
   code first **regresses out the known covariates + genotype PCs**, then takes the
   top 5 PCs of the *residual*. So the methylation PCs capture only the
   **remaining hidden structure** (unknown batch, cell-composition, etc.) —
   PEER-style latent factors without the redundancy.

**Key code.** `residual_methylation_pcs()`:
- mean-centres the samples × CpGs matrix,
- if known covariates exist, builds `Cc = [1 | C]` and removes their fit via
  `lstsq` (`X − Cc·β`) — i.e. residualises,
- forms the samples × samples Gram matrix `G = XXᵀ` and eigendecomposes it
  (`eigh`), taking the top-k eigenvectors as the PCs. (Working on the Gram matrix
  is the efficient route when CpGs ≫ samples.)

**Output.** Covariates as **rows**, samples as **columns** (tensorQTL reads it
transposed). For the pilot: 12 covariates = sex + platform_R10 + 5 genoPC +
5 methPC, leaving residual df ≈ 86 at n=100.

**Defend it.** "Known technical and biological covariates and ancestry are modelled
explicitly; latent confounders are captured as principal components of the
*residual* methylation, so no source of structure is double-counted. The number of
covariates is kept well below n to preserve residual degrees of freedom."

---

# Stage 05 — cis-mQTL mapping (tensorQTL, GPU)

All four sub-stages share the same engine, `tensorqtl`, run on a V100 GPU. Two
transforms recur and are worth stating once:

- **Inverse-normal transform** (`inverse_normal_transform`): each CpG's values are
  rank-transformed to a standard normal (`norm.ppf((rank−0.5)/n)`). This makes the
  phenotype Gaussian per-site, so the linear model's p-values are well-calibrated
  and robust to outliers/skew. Applied consistently across **05a, 05c, 05d** so
  effect sizes and p-values are comparable across stages.
- **±1 Mb cis-window, MAF ≥ 0.05**: only variants within 1 Mb of the CpG are tested
  (cis), and only common variants (rare variants are underpowered at n=100).

## `05a_run_permutation.py` — the permutation pass
**Purpose.** One empirical (permutation) p-value per CpG, asking "does *any*
variant in this CpG's cis-window associate with its methylation, beyond chance?"

**Why permutation.** Each CpG is tested against thousands of correlated cis
variants. A naïve minimum p-value would be wildly anti-conservative. tensorQTL
permutes the sample labels many times to build the null distribution of the
*best* association per CpG, then fits a **beta distribution** to it
(`pval_beta`) — an accurate, smooth per-CpG empirical p-value.

**Key code.** `cis.map_cis(..., nperm=, window=, maf_threshold=, seed=)`. A
**chrom-mismatch guard** aborts if genotype and phenotype chromosome names are
disjoint (the #1 cause of "0 variants in window"). Pilot used `nperm=1000`; the
final run uses `nperm=10000` for smoother tails.

**Output.** One row per CpG with `pval_beta` (and the lead variant). Result:
**6,995 significant CpGs** survive FDR (next stage), genomic inflation **λ=1.446**.

## `05b_fdr_significant.py` — multiple-testing control
**Purpose.** Convert per-CpG p-values to q-values and select significant CpGs at
FDR < 0.05.

**Key logic.** Tries tensorQTL's **Storey q-values** (`post.calculate_qvalues`,
needs R's `qvalue`); if unavailable, falls back to a dependency-free
**Benjamini–Hochberg** implementation in numpy (`bh_qvalues`). Both control the
false discovery rate; Storey additionally estimates π₀ (the null proportion) for
slightly more power. Writes the FDR table and the **significant-CpG list** that
every downstream stage consumes.

**Defend it.** "We control FDR, not family-wise error, because the goal is a
discovery set for follow-up, not a single confirmatory test; BH is the
conservative fallback if the Storey estimator isn't available."

## `05c_run_nominal_significant.py` — per-pair statistics
**Purpose.** For the significant CpGs only, compute the **full per-variant**
statistics: `pval_nominal`, `slope` (effect size), `slope_se`, `af`.

**Why.** The permutation pass gives one p per CpG; fine-mapping and prioritisation
need the *individual* variant effect sizes and p-values. Running `map_nominal`
only on the significant CpGs (not all 600k) keeps it cheap.

**Output.** A parquet of all cis pairs for the 6,995 CpGs (~10 M rows). This feeds
the `neg_log10_p` and `abs_effect` features in Stage 09.

## `05d_permuted_null.py` — the calibration negative control
**Purpose.** Decide whether the observed inflation (λ=1.446) is a real
genotype–phenotype signal or an artefact of a mis-calibrated test.

**Key logic.** It **permutes the genotype sample labels** — each sample keeps its
own methylation and covariates but is paired with a *random* sample's genotypes.
This destroys any true genotype–phenotype association. Re-running the permutation
scan on a 20k-CpG random subset under this null should give **λ ≈ 1** if the test
is well-calibrated.

**Interpretation (state this carefully).** A flat null (λ≈1) confirms the test
itself is calibrated, so the real λ=1.446 reflects genuine genotype–phenotype
structure — *not* a software/df artefact. It does **not**, by itself, separate
true cis-mQTL signal from relatedness, because shuffling destroys relatedness-driven
association too. Distinguishing those is the job of an unrelated-subset or
kinship-model rerun, which matters at full 1KGP.

---

# Stage 06 — Fine-mapping: `finemap_susie.py`

**Purpose.** For each significant CpG, move from "there's a signal in this window"
to "here are the specific variants likely responsible, with probabilities."

**Why SuSiE.** Within a cis-window many variants are in tight linkage
disequilibrium (LD) and all show similar association — you can't tell which is
causal from marginal p-values. **SuSiE** (Sum of Single Effects) is a Bayesian
fine-mapping method that decomposes the signal into up to **L** independent causal
effects and returns:
- a **PIP** (posterior inclusion probability) per variant — P(this variant is
  causal), and
- **95 % credible sets** — minimal variant groups that collectively contain a
  causal variant with 95 % probability.

**Key code.** Uses tensorQTL's **native torch SuSiE**:
`susie.map(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df,
L=10, maf_threshold=0.05, window=, coverage=0.95, summary_only=True)`.
`L=10` allows up to 10 independent signals per CpG; `coverage=0.95` defines the
credible-set threshold. The parser is defensive (handles both the newer DataFrame
return and older dict return) — the function name and return shape changed across
tensorQTL releases, which is a real gotcha worth a sentence in the methods.

**Output.** `[phenotype_id, variant_id, pip, af, cs_id]`. Pilot result: **6,797 of
6,995 CpGs (97 %)** resolved at least one 95 % credible set; **80,711**
credible-set variant rows. These credible-set variants are the candidate set the
prioritiser ranks.

**Defend it.** "Marginal association can't localise a causal variant inside an LD
block; SuSiE provides per-variant posterior probabilities and credible sets,
which is the principled basis for prioritisation. At n=100 credible sets are wide,
so we nominate rather than claim causality."

---

# Stage 07 — Structural-variant integration (novel): `sv_annotate.py`

**Purpose.** Ask, for each significant CpG, whether a **structural variant** (SV —
deletion, insertion, duplication, inversion) could explain the mQTL — something a
SNV-only scan would miss entirely. This is the novel contribution.

**Why it matters.** SNP arrays and standard QTL pipelines ignore SVs, yet SVs
(especially deletions/insertions overlapping regulatory DNA) are strong candidate
causal variants for methylation changes. ONT long reads make SVs callable, so the
framework can test them.

**Key logic — two complementary tests per CpG:**
1. **Direct association.** `residualise()` removes the covariates from the CpG's
   methylation, then `ols_p()` regresses the residual on each cis SV's **dosage**
   (0/1/2 alt alleles, via `gt_to_dosage`) → an SV-mQTL **β and p-value** by a
   proper t-test.
2. **LD tagging.** `r²` (squared correlation) between the CpG's **lead SNV** dosage
   and each cis SV dosage. A high r² means the SV is a plausible causal candidate
   tagged by the lead SNP.

A CpG-SV pair is **flagged** if the SV is directly associated (`p < 1e-3`) **or**
in high LD with the lead SNV (`r² ≥ 0.8`).

**Inputs.** Methylation matrix, covariates, the permutation table (for the lead
SNV per CpG), the significant-CpG list, the SNV plink set (lead dosage + LD), and
the SV dosage matrix from `extract_sv_genotypes.pbs` (which streams chr22 SVs from
the 1KGP `SNV_INDEL_SV` panel, filtering symbolic ALT alleles `<DEL>`, `<INS>`…).

**Output.** `[phenotype_id, …, sv_id, svtype, sv_pos, distance, sv_beta, sv_p,
r2_lead, sv_implicated]`. Pilot: **1,054 flagged pairs across 885 CpGs (~13 % of
significant CpGs SV-implicated)**; standout cluster at chr22:49.94 Mb.

**The honest caveat (say this).** At n=100, an SV in r²≈1 with a SNV is
**statistically indistinguishable** from it — the framework *nominates* the SV as a
candidate; it does not prove the SV (rather than the linked SNV) is causal. The
code comment flags the next step: recode SVs to biallelic dosages, merge with the
SNV set, and re-run SuSiE so SVs compete with SNVs for PIP.

---

# Stage 08 — Functional annotation: `functional_annotate.py`

**Purpose.** Attach regulatory context to each significant CpG (and, in Stage 09,
to each candidate variant): nearest gene/TSS and distance, CpG-island overlap,
ENCODE candidate cis-regulatory element (cCRE) overlap.

**Why.** A variant/CpG sitting in a promoter, enhancer, or CpG island is more
plausibly functional. This is orthogonal evidence to the statistics — it doesn't
depend on the association at all.

**Key logic (pure pandas, no bedtools — fine at chr22 scale).**
- `nearest gene/TSS`: TSS starts are sorted; `np.searchsorted` finds the insertion
  point, then the nearer of the two flanking TSSs gives `nearest_gene` and signed
  `dist_to_tss`.
- `overlaps()`: interval containment test for CpG islands and cCREs, returning the
  element label (e.g. a cCRE class like `CA-H3K4me3`).

**Inputs.** Annotation tracks downloaded by `get_annotations.sh`: GENCODE v44 TSS,
UCSC `cpgIslandExt`, ENCODE SCREEN GRCh38 cCREs — all subset to chr22.

**Output.** Per-CpG: `nearest_gene, dist_to_tss, cpg_island, ccre`. Pilot: all
6,995 CpGs annotated.

---

# Stage 09 — Prioritisation

Two scripts: a simple transparent scorer and the full framework matching your
methods §3.1.6.

## `prioritise.py` (simple, transparent)
A single explicit linear score per credible-set variant:
`score = w_pip·PIP + w_func·functional + w_sv·SV`, with the functional weight set
by cCRE membership, TSS distance < 2 kb, or CpG-island overlap. Deliberately
simple so the ranking is fully auditable. Useful as a sanity baseline.

## `prioritise_framework.py` (the dissertation framework)
**Purpose.** Build a standardised multi-evidence feature vector per
credible-set variant and produce **two** scores plus a robustness check.

**Candidate set.** Variants in a SuSiE 95 % credible set (`--cs-only`, default on).

**Features assembled (the §3.1.6 vector):**
- `pip` — fine-mapping posterior (from SuSiE).
- `maf` — minor allele frequency (`min(af, 1−af)`).
- `prox_cpg` — `−log10(|var_pos − cpg_pos| + 1)`; closeness of the variant to the
  CpG it regulates.
- `sv_flag` — is the CpG SV-implicated (from Stage 07)?
- `var_ccre`, `var_cpg_island` — does the **variant** overlap a cCRE / CpG island?
  (computed with `overlap_flag` + merged intervals, at variant position).
- `prox_tss` — `−log10(distance to nearest TSS + 1)` at the variant.
- `neg_log10_p`, `abs_effect` — association strength and effect size (from the
  Stage 05c nominal parquet; **filtered read** pulls only credible-set variants and
  4 columns so memory stays small and it scales genome-wide).
- `neg_log10_fdr` — the CpG-level significance (from Stage 05b).

**Two scores:**
1. **Weighted additive** (Ionita-Laza style). Every feature is **z-scored**
   (`zscore`) so they're on a common scale, then combined with explicit
   `DEFAULT_WEIGHTS` (e.g. `pip=2.0`). Transparent and reproducible — anyone can
   read the weights.
2. **Elastic-net** (Gagliano style). An `ElasticNetCV` (cross-validated, mixing
   L1/L2) **regresses PIP on the other features**, learning data-driven weights
   instead of hand-set ones. Because there's no causal label, PIP is used as the
   self-supervised target: the model learns which annotations *predict* fine-mapping
   confidence. Reports `alpha`, `l1_ratio`, `R²`, and the learned coefficients.

**Robustness check.** **Spearman ρ** between the weighted-additive ranking and the
elastic-net ranking. High convergence ⇒ the ranking isn't an artefact of arbitrary
weights, and isn't redundant either. Pilot result: **ρ = 0.629**, elastic-net
**R² = 0.17**, with `neg_log10_p` (+) and `neg_log10_fdr` (−) the dominant learned
drivers (the negative fdr coefficient reflects strongly-significant CpGs having
broader credible sets, so PIP spreads across more variants).

**Output.** `chr22.prioritised_framework.txt.gz` — 80,711 ranked variant-CpG pairs
with both scores and the z-scored features (`z_` columns) for transparency.

**Defend it.** "With no causal ground truth we can't train a supervised classifier,
so we use two complementary weighting philosophies — fixed and learned — and treat
their *convergence* as the evaluation. The framework nominates and ranks; it does
not assert causality."

---

# Supporting tools

- **`qq_fdr.py`** — computes the **genomic inflation factor λ** (median observed
  χ² ÷ median expected χ²; λ≈1 is calibrated, λ>1 is inflation), applies BH-FDR,
  and draws the QQ plot. This is how λ=1.446 (real) and the null λ (control) are
  measured identically.
- **PBS wrappers** (`run_*_gpu.pbs`, `*.pbs`) — schedule the jobs on Gadi. The
  recurring pattern: request `gpuvolta` + 1 V100, set the micromamba env paths,
  redirect HOME-bound caches to `/scratch` (the env quota fix), and `qsub`. They
  contain no science — they're reproducible run records.

---

# Concepts you should be ready to define in a viva

| Concept | One-line definition you can give |
|---|---|
| **beta vs M-value** | beta = fraction methylated [0,1]; M-value = logit of beta, variance-stabilised for linear models. |
| **cis-window** | the ±1 Mb region around a CpG within which variants are tested (local regulation). |
| **inverse-normal transform** | rank-to-normal per CpG, so p-values are calibrated and outlier-robust. |
| **genotype PCs** | principal components of genotypes capturing ancestry/population structure (a confounder). |
| **residual methylation PCs** | latent factors from methylation *after* removing known covariates — hidden batch/cell-composition without double-counting. |
| **permutation p-value / pval_beta** | empirical per-CpG significance from permuting labels and fitting a beta distribution to the null of the best hit. |
| **FDR (BH / Storey)** | expected proportion of false positives among discoveries; controlled instead of family-wise error for a discovery set. |
| **genomic inflation λ** | median test statistic ÷ its null expectation; ≈1 = calibrated, >1 = inflation (signal, relatedness, or confounding). |
| **PIP** | posterior probability a variant is causal, from SuSiE fine-mapping. |
| **95 % credible set** | smallest variant group containing a causal variant with 95 % probability. |
| **LD / r²** | correlation² between two variants' dosages; high r² ⇒ statistically indistinguishable. |
| **SV dosage** | 0/1/2 alternate-allele count for a structural variant, treated like a SNP for association. |
| **cCRE** | ENCODE candidate cis-regulatory element (promoter/enhancer-like region). |
| **z-scoring** | rescaling each feature to mean 0, sd 1 so heterogeneous features combine fairly. |
| **elastic-net** | regularised regression mixing L1 (sparsity) and L2 (shrinkage); here learns annotation weights predicting PIP. |
| **Spearman convergence** | rank-correlation between the two scores; the robustness criterion in lieu of ground truth. |

---

# The three claims this framework can honestly support

1. **It is calibrated.** Covariates + ancestry PCs + residual factors control
   confounding; the permuted-null check tests that the inflation isn't a software
   artefact.
2. **It localises and integrates.** SuSiE narrows each signal to credible sets;
   SV, functional, and association evidence are combined into a single transparent
   ranking.
3. **It is reproducible and scalable.** Every stage is scripted, every intermediate
   is on disk, and the heavy reads are filtered — so the same code runs from a
   100-sample chr22 pilot to full 1KGP and on to PROPHECY.

What it does **not** claim: causality of any individual variant. At n=100, wide
credible sets and r²≈1 ties mean the output is a prioritised shortlist for
follow-up — which is exactly what §3.1.6 sets out to deliver.
