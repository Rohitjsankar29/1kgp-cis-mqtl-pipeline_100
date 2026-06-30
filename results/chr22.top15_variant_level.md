# Top 15 prioritised variants — variant-level (deduplicated)

One row per distinct candidate variant (collapsed from the pair-level table with
`scripts/09_prioritise/collapse_to_variant_level.py`). `n_cpgs` = number of CpGs
the variant is the credible-set lead for; `cpg_snp` = variant lies on the CpG with
a CpG-altering change (possible methylation-calling artefact — review/exclude).

| variant_id | best_cpg | n_cpgs | cpg_snp | pip | pval_nominal | slope | score_weighted | score_enet |
| --- | --- | ---: | :---: | ---: | --- | ---: | ---: | ---: |
| chr22:44360552:G:A | chr22_44360552 | 1 | **yes** | 1.000 | 1.8e-32 | -1.199 | 27.696 | 0.322 |
| chr22:47918714:G:A | chr22_47918793 | **83** | no | 1.000 | 2.2e-38 | -1.380 | 27.232 | 0.274 |
| chr22:23924767:A:C | chr22_23924893 | **143** | no | 0.983 | 2.1e-36 | 1.161 | 26.137 | 0.274 |
| chr22:17404078:T:G | chr22_17404114 | 6 | no | 1.000 | 1.2e-30 | -1.216 | 25.854 | 0.243 |
| chr22:43233182:C:T | chr22_43233179 | 1 | no | 1.000 | 2.9e-33 | -1.147 | 25.296 | 0.302 |
| chr22:27196780:G:T | chr22_27196780 | 1 | no | 1.000 | 4.3e-32 | -1.157 | 24.328 | 0.349 |
| chr22:49213760:G:A | chr22_49213760 | 1 | **yes** | 1.000 | 1.5e-27 | -1.105 | 23.938 | 0.307 |
| chr22:45100712:C:G | chr22_45100471 | 1 | no | 0.999 | 4.3e-32 | -1.430 | 23.561 | 0.118 |
| chr22:37158799:C:G | chr22_37158881 | 11 | no | 1.000 | 3.7e-32 | -1.174 | 23.544 | 0.246 |
| chr22:23927818:T:C | chr22_23925160 | 2 | no | 1.000 | 8.5e-25 | -1.349 | 23.505 | 0.246 |
| chr22:23908608:C:G | chr22_23919778 | 18 | no | 1.000 | 2.9e-24 | 1.148 | 22.540 | 0.130 |
| chr22:27253743:A:G | chr22_27253743 | 1 | **yes** | 1.000 | 1.0e-26 | 1.125 | 21.695 | 0.275 |
| chr22:40728988:C:T | chr22_40728788 | 14 | no | 0.999 | 6.2e-21 | 1.453 | 21.349 | 0.241 |
| chr22:48238411:T:C | chr22_48238408 | 2 | no | 1.000 | 2.8e-23 | 1.405 | 21.257 | 0.268 |
| chr22:38427131:T:C | chr22_38427128 | 2 | no | 0.997 | 4.1e-18 | 1.200 | 21.206 | 0.296 |

**Reading it.** Two variants tag large co-methylated blocks — `chr22:23924767:A:C`
(143 CpGs) and `chr22:47918714:G:A` (83 CpGs) — candidate broad regulatory
variants. Three variants are CpG-SNP-flagged (including the top-ranked one); all
three have `n_cpgs = 1`, the signature of a local calling artefact rather than a
regulatory domain. See `docs/PRIORITISATION_NOTE.md` for how these are handled.
