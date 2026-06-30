# Prioritisation note — variant-level collapse and CpG-SNP screen

A correction and refinement to the stage-09 prioritisation output, applied after
review.

## The issue

The prioritised table (`chr22.prioritised_framework.txt.gz`) ranks variant–CpG
**pairs**. A variant that regulates a cluster of nearby, co-methylated CpGs is the
credible-set lead for each of them, so it appears once per CpG and floods the top
of the ranking. In the top 30 pair-rows there were only **13 distinct variants**;
two of them dominated — `chr22:47918714:G:A` (10 rows) and `chr22:23924767:A:C`
(5 rows). The shortlist should be a list of distinct candidate *variants*, so the
unit was wrong (pairs, not variants).

## The fix

`scripts/09_prioritise/collapse_to_variant_level.py` collapses the table to one
row per distinct `variant_id`, keeping the variant's best-scoring CpG as
representative and adding:

- **`n_cpgs`** — how many distinct CpGs the variant is the credible-set lead for.
  This turns the duplication into information: the two dominant variants tag
  **83** and **143** CpGs respectively — i.e. they are candidate *broad regulatory
  variants* spanning large co-methylated blocks, not repeated single hits.
- **`cpg_snp`** — flags variants lying on the CpG (within 1 bp) with a
  CpG-altering change (C↔T / G↔A).

Output: `chr22.prioritised_variant_level.txt.gz` (see
`results/chr22.top15_variant_level.md` for the top 15).

## CpG-SNP screen

Three top variants are CpG-SNP-flagged, including the highest-ranked,
`chr22:44360552:G:A`. A SNP that abolishes the CpG (e.g. CG→CA) means reads
carrying the alt allele have no CpG to methylate, so genotype *mechanically*
determines the readout — an ultra-strong but **artefactual** association rather
than regulation. Consistently, all three flagged variants have `n_cpgs = 1` (the
local-artefact signature, versus the broad-regulator signature of the 83/143
variants).

**Handling (flag and separate, not silent deletion).** CpG-SNPs are reported but
clearly marked; the primary candidate shortlist is read *excluding* them. They are
not deleted, because a CpG-SNP can occasionally be the genuine causal variant —
flagging leaves the decision explicit and auditable.

**Caveat.** The flag is a heuristic (±1 bp window, C↔T / G↔A). Whether a given
allele change destroys the CpG depends on the strand convention of the phenotype
position (C vs G of the CpG); confirming that convention is what turns "probable"
into "confirmed" artefact.
