# Results — chr22 pilot (post CpG-SNP filter)

Prioritised cis-mQTL variants for chr22 (n=100), after removing CpG-SNP artefacts
(2.9% of CpGs whose C/G overlaps an in-sample SNP) and integrating the
allele-specific methylation (ASM) evidence layer.

**Summary:** 20,936 candidate credible-set variants across 6,539 FDR-significant CpGs ·
**67 SD-ASM corroborated** (genetic, allele-consistent) · **7 imprint** (flagged out).

Columns: **PIP** fine-mapping posterior · **Score** weighted multi-evidence score ·
**cCRE** ENCODE candidate regulatory element · **SV** structural-variant implicated ·
**#CpGs** CpGs tagged (regulatory footprint) · **ASM** SD-ASM/imprint/insufficient ·
**ASM het** heterozygous carriers showing allele-specific methylation.

## Top 15 prioritised variants

| Variant            | Gene            | cCRE       | SV   |   #CpGs |   PIP |   Score | ASM          |   ASM het |
|:-------------------|:----------------|:-----------|:-----|--------:|------:|--------:|:-------------|----------:|
| chr22:42500468:T:C | SERHL2          | dELS       | Y    |      10 | 0.899 |   17.85 | insufficient |         0 |
| chr22:37080103:C:T | ENSG00000231467 | PLS        | Y    |       5 | 0.981 |   17.82 | insufficient |         0 |
| chr22:45413743:C:A | RIBC2           | pELS       |      |      32 | 1     |   17.37 | insufficient |         0 |
| chr22:45163750:C:G | NUP50-DT        | none       |      |       1 | 0.998 |   17.12 | insufficient |         0 |
| chr22:36552356:G:C | ENSG00000229971 | none       | Y    |       2 | 0.979 |   17.1  | insufficient |         0 |
| chr22:38427131:T:C | ENSG00000228620 | dELS       |      |       2 | 0.997 |   17.01 | insufficient |         0 |
| chr22:43281649:C:T | SCUBE1-AS1      | CA         | Y    |       1 | 1     |   16.99 | insufficient |         0 |
| chr22:23763821:C:G | C22orf15        | pELS       | Y    |      27 | 0.998 |   16.89 | insufficient |         0 |
| chr22:39124391:A:G | COX5BP7         | dELS       | Y    |       6 | 1     |   16.88 | insufficient |         1 |
| chr22:23589994:G:A | ENSG00000272733 | none       | Y    |       1 | 1     |   16.71 | insufficient |         0 |
| chr22:45908328:C:T | ENSG00000235091 | CA-TF      | Y    |       1 | 1     |   16.59 | insufficient |         0 |
| chr22:45990729:A:G | WNT7B           | CA-H3K4me3 | Y    |       1 | 0.996 |   16.52 | insufficient |         0 |
| chr22:19291671:G:T | ENSG00000287146 | pELS       |      |       5 | 0.947 |   16.52 | insufficient |         2 |
| chr22:49347680:C:T | ENSG00000285722 | CA         | Y    |       1 | 1     |   16.46 | insufficient |         0 |
| chr22:17404078:T:G | ENSG00000229492 | dELS       | Y    |       6 | 0.996 |   16.44 | insufficient |         1 |

## SD-ASM corroborated cis-mQTLs (top 15 of 67)

These candidates are independently corroborated by sequence-dependent allele-specific
methylation — the methylation tracks the variant allele *within* heterozygous
individuals (allele_consistency = 1.0), a within-molecule line of evidence distinct
from the between-individual QTL association.

| Variant            | Gene            | cCRE   | SV   |   #CpGs |   PIP |   Score | ASM    |   ASM het |
|:-------------------|:----------------|:-------|:-----|--------:|------:|--------:|:-------|----------:|
| chr22:46773373:C:T | TBC1D22A        | pELS   |      |      10 | 0.489 |    7.14 | SD-ASM |        27 |
| chr22:42819686:A:G | ENSG00000274717 | none   | Y    |       8 | 0.124 |    4.77 | SD-ASM |        22 |
| chr22:42529125:A:G | RRP7A           | dELS   | Y    |      18 | 0.998 |   13.79 | SD-ASM |        21 |
| chr22:35299809:G:A | ENSG00000273176 | pELS   | Y    |       3 | 0.025 |    7.66 | SD-ASM |        15 |
| chr22:26924029:C:G | ENSG00000235271 | none   |      |       2 | 0.062 |    0.24 | SD-ASM |        15 |
| chr22:45247026:T:C | KIAA0930        | dELS   |      |      11 | 1     |   13.36 | SD-ASM |        15 |
| chr22:37245956:G:A | RAC2            | none   | Y    |       5 | 0.064 |    4.34 | SD-ASM |        14 |
| chr22:42501739:G:A | SERHL2          | dELS   | Y    |      12 | 0.752 |   10.32 | SD-ASM |        13 |
| chr22:23537508:C:T | PCAT14          | pELS   |      |      11 | 1     |   14.62 | SD-ASM |        11 |
| chr22:43154143:T:C | TSPO            | dELS   |      |       4 | 0.333 |    5.78 | SD-ASM |        10 |
| chr22:38761399:C:T | ENSG00000273076 | dELS   |      |       7 | 0.969 |   12.66 | SD-ASM |         9 |
| chr22:49658959:C:A | MIR3667HG       | none   |      |      10 | 0.506 |    7.12 | SD-ASM |         9 |
| chr22:44711245:A:G | PRR5-ARHGAP8    | dELS   |      |       6 | 0.998 |   13.27 | SD-ASM |         9 |
| chr22:43148817:C:T | TSPO            | none   |      |       1 | 0.463 |    6.88 | SD-ASM |         8 |
| chr22:29778240:A:G | ENSG00000287967 | none   |      |       8 | 0.571 |    6.81 | SD-ASM |         8 |

> **Note on ASM coverage.** ASM is short-range (SNP and CpG must sit on the same reads)
> and heterozygote-dependent, so it corroborates a *subset* of candidates; top-scored
> variants often show `insufficient` where too few het carriers have fragment coverage.
> The table above shows where the association and ASM layers converge.

## Full tables
- `results/chr22_prioritised_POSTfilter.html` — interactive, sortable (top 500)
- `results/chr22_top30_POSTfilter.md` — top 30 by score
- `results/chr22_SDASM_POSTfilter.md` — all 67 SD-ASM variants
- `results/chr22.prioritised_variant_level.enriched_asm.txt.gz` — full 20,936-variant table

## Method note — CpG-SNP filter
Variants that overlap the C or G of a CpG mechanically force the methylation readout
(genotype → apparent methylation), producing artefactual associations that dominated
the unfiltered ranking. The filter (`scripts/05_tensorqtl/04b_filter_cpg_snps.py`)
removes these before tensorQTL. Effect: significant CpGs 6,995 → 6,539; the former
top-ranked variant (a CpG-SNP artefact) drops out; elastic-net R² 0.17 → 0.09,
reflecting removal of the artificially-predictable artefacts.
