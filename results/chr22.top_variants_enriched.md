# Top candidate variants — enriched shortlist (chr22, n=100)

Variant-level shortlist with statistical rank **and** biological context in each
row. Built by `scripts/09_prioritise/collapse_to_variant_level.py` →
`scripts/09_prioritise/enrich_variant_table.py`. Full ranked + enriched table:
`results/chr22.prioritised_variant_level.enriched.txt.gz`.

Columns: `n_cpgs` = CpGs the variant leads (breadth); `cCRE` = ENCODE regulatory
class (dELS/pELS = distal/proximal enhancer-like, CA = chromatin-accessible);
`SV` = structural variant implicated at that CpG; `cpg_snp` = on-CpG
artefact flag; `PIP` = fine-mapping posterior; score = primary ranking.

| variant_id | nearest_gene | dist_to_tss | cCRE | n_cpgs | SV | cpg_snp | PIP | score |
| --- | --- | ---: | :---: | ---: | :---: | :---: | ---: | ---: |
| chr22:44360552:G:A | ENSG00000220702 | -4,999 | dELS | 1 | INS | **yes** | 1.000 | 27.70 |
| chr22:47918714:G:A | ENSG00000279712 | 2,065 | dELS | **83** | — | no | 1.000 | 27.23 |
| chr22:23924767:A:C | ENSG00000225282 | -2,007 | none | **143** | DUP | no | 0.983 | 26.14 |
| chr22:17404078:T:G | ENSG00000229492 | -14,583 | dELS | 6 | DEL,INS | no | 1.000 | 25.85 |
| chr22:43233182:C:T | ENSG00000234892 | 1,038 | pELS | 1 | — | no | 1.000 | 25.30 |
| chr22:27196780:G:T | ENSG00000238195 | -7,431 | dELS | 1 | — | no | 1.000 | 24.33 |
| chr22:49213760:G:A | RPL35P8 | -487 | none | 1 | DUP | **yes** | 1.000 | 23.94 |
| chr22:45100712:C:G | ENSG00000273243 | -55,540 | CA | 1 | — | no | 0.999 | 23.56 |
| chr22:37158799:C:G | ENSG00000235237 | -7,876 | dELS | 11 | — | no | 1.000 | 23.54 |
| chr22:23927818:T:C | ENSG00000225282 | -1,740 | none | 2 | DUP | no | 1.000 | 23.50 |
| chr22:23908608:C:G | ENSG00000225282 | -7,122 | dELS | 18 | DUP | no | 1.000 | 22.54 |
| chr22:27253743:A:G | ENSG00000233574 | -22,497 | dELS | 1 | — | **yes** | 1.000 | 21.69 |
| chr22:40728988:C:T | ENSG00000289292 | 43,849 | dELS | 14 | DEL | no | 0.999 | 21.35 |
| chr22:48238411:T:C | ENSG00000285707 | -20,708 | none | 2 | DEL | no | 1.000 | 21.26 |

*(plus further variants in the full enriched file)*

## What this shows

- **Regulatory enrichment:** most top variants overlap an ENCODE enhancer-like
  element (dELS/pELS) or accessible chromatin — the candidates sit in regulatory
  DNA, not at random, which independently supports them being functional.
- **Broad regulators corroborated:** the two highest-`n_cpgs` variants
  (143 and 83 CpGs) both sit near a TSS / in an enhancer — breadth and regulatory
  context agree.
- **Artefact screen:** three variants are CpG-SNP-flagged (incl. the top score);
  all three have `n_cpgs = 1`, the local-artefact signature. Handled as flag-and-
  separate (see `docs/PRIORITISATION_NOTE.md`).

## Caveat to verify

`cpg_island` was uniformly "no" across the shown rows; confirm the CpG-island
track loaded correctly (a healthy yes/no mix across all 6,995 CpGs) before relying
on that column.
