# chr22 Prioritised Variants + ASM (post CpG-SNP filter)

_20,936 candidate variants across the filtered significant CpGs. **67 SD-ASM** corroborated, 7 imprint. Top 50 by weighted score shown; full table in `chr22.enriched_asm.POSTfilter.tsv`._

**Columns:** PIP = fine-mapping posterior · Score = weighted multi-evidence score · cCRE = ENCODE regulatory element · SV = structural-variant implicated · #CpGs = CpGs tagged · TSS_dist = bp to nearest TSS · ASM = SD-ASM/imprint/insufficient · ASM_het = het carriers with ASM · ASM_delta = allele methylation difference · ASM_cons = allele consistency.

| Variant            | Gene            | cCRE       | SV   | Island   |   #CpGs |   TSS_dist |   PIP |   Score | ASM          |   ASM_het |   ASM_delta |   ASM_cons |
|:-------------------|:----------------|:-----------|:-----|:---------|--------:|-----------:|------:|--------:|:-------------|----------:|------------:|-----------:|
| chr22:42500468:T:C | SERHL2          | dELS       | Y    |          |      10 |      -3882 | 0.899 |   17.85 | insufficient |         0 |     nan     |        nan |
| chr22:37080103:C:T | ENSG00000231467 | PLS        | Y    |          |       5 |         69 | 0.981 |   17.82 | insufficient |         0 |     nan     |        nan |
| chr22:45413743:C:A | RIBC2           | pELS       |      | Y        |      32 |        354 | 1     |   17.37 | insufficient |         0 |     nan     |        nan |
| chr22:45163750:C:G | NUP50-DT        | none       |      |          |       1 |       -711 | 0.998 |   17.12 | insufficient |         0 |     nan     |        nan |
| chr22:36552356:G:C | ENSG00000229971 | none       | Y    |          |       2 |        114 | 0.979 |   17.1  | insufficient |         0 |     nan     |        nan |
| chr22:38427131:T:C | ENSG00000228620 | dELS       |      | Y        |       2 |       2840 | 0.997 |   17.01 | insufficient |         0 |     nan     |        nan |
| chr22:43281649:C:T | SCUBE1-AS1      | CA         | Y    |          |       1 |       5696 | 1     |   16.99 | insufficient |         0 |     nan     |        nan |
| chr22:23763821:C:G | C22orf15        | pELS       | Y    |          |      27 |        759 | 0.998 |   16.89 | insufficient |         0 |     nan     |        nan |
| chr22:39124391:A:G | COX5BP7         | dELS       | Y    |          |       6 |       3833 | 1     |   16.88 | insufficient |         1 |      -1     |          1 |
| chr22:23589994:G:A | ENSG00000272733 | none       | Y    |          |       1 |       6126 | 1     |   16.71 | insufficient |         0 |     nan     |        nan |
| chr22:45908328:C:T | ENSG00000235091 | CA-TF      | Y    |          |       1 |      16887 | 1     |   16.59 | insufficient |         0 |     nan     |        nan |
| chr22:45990729:A:G | WNT7B           | CA-H3K4me3 | Y    |          |       1 |      13562 | 0.996 |   16.52 | insufficient |         0 |     nan     |        nan |
| chr22:19291671:G:T | ENSG00000287146 | pELS       |      |          |       5 |        574 | 0.947 |   16.52 | insufficient |         2 |      -0.796 |          1 |
| chr22:49347680:C:T | ENSG00000285722 | CA         | Y    |          |       1 |     -25247 | 1     |   16.46 | insufficient |         0 |     nan     |        nan |
| chr22:17404078:T:G | ENSG00000229492 | dELS       | Y    |          |       6 |     -14612 | 0.996 |   16.44 | insufficient |         1 |      -1     |          1 |
| chr22:49392937:C:A | ENSG00000285722 | CA-CTCF    | Y    |          |       1 |      20007 | 1     |   16.42 | insufficient |         0 |     nan     |        nan |
| chr22:17219920:C:T | FAM32BP         | pELS       | Y    |          |       2 |       6438 | 0.993 |   16.35 | insufficient |         0 |     nan     |        nan |
| chr22:40685160:G:T | ENSG00000289292 | PLS        | Y    |          |       1 |        192 | 0.9   |   16.31 | insufficient |         0 |     nan     |        nan |
| chr22:23927818:T:C | ENSG00000225282 | pELS       | Y    |          |       1 |        -93 | 0.999 |   16.31 | insufficient |         0 |     nan     |        nan |
| chr22:43836136:T:C | EFCAB6-DT       | dELS       | Y    |          |       2 |      23688 | 0.999 |   16.3  | insufficient |         0 |     nan     |        nan |
| chr22:44740243:C:T | ARHGAP8         | dELS       | Y    |          |       1 |     -12322 | 0.977 |   16.26 | insufficient |         0 |     nan     |        nan |
| chr22:25405434:A:G | ENSG00000231466 | none       |      |          |       2 |       6827 | 0.998 |   16.23 | insufficient |         0 |     nan     |        nan |
| chr22:42440331:A:C | ENSG00000230107 | none       | Y    |          |       8 |       2731 | 0.999 |   16.21 | insufficient |         0 |     nan     |        nan |
| chr22:49357798:C:T | ENSG00000285722 | CA         | Y    |          |       1 |     -15153 | 1     |   16.21 | insufficient |         0 |     nan     |        nan |
| chr22:50178377:A:G | TRABD           | dELS       |      | Y        |       4 |      -7584 | 1     |   16.2  | insufficient |         0 |     nan     |        nan |
| chr22:47918714:G:A | ENSG00000279712 | none       | Y    |          |      83 |       1408 | 1     |   16.2  | insufficient |         1 |       0.718 |          1 |
| chr22:48603783:A:G | ENSG00000281732 | dELS       | Y    |          |       1 |      56388 | 1     |   16.04 | insufficient |         0 |     nan     |        nan |
| chr22:44858837:C:G | RPL6P28         | none       | Y    |          |       1 |     -42697 | 1     |   16.03 | insufficient |         0 |     nan     |        nan |
| chr22:40728988:C:T | ENSG00000289292 | dELS       | Y    |          |      14 |      44028 | 0.998 |   15.9  | insufficient |         0 |     nan     |        nan |
| chr22:44612034:A:G | LINC00229       | none       | Y    |          |       1 |     -13576 | 1     |   15.84 | insufficient |         0 |     nan     |        nan |
| chr22:36946958:C:A | CSF2RBP1        | dELS       | Y    |          |       3 |      -6237 | 1     |   15.83 | insufficient |         0 |     nan     |        nan |
| chr22:36070934:T:C | ENSG00000287269 | pELS       | Y    |          |       1 |      -1124 | 0.979 |   15.71 | insufficient |         0 |     nan     |        nan |
| chr22:23908608:C:G | ENSG00000225282 | dELS       | Y    |          |      18 |      -7122 | 1     |   15.65 | none         |         0 |     nan     |        nan |
| chr22:38570408:C:G | NPTXR           | dELS       |      |          |       1 |       9425 | 0.999 |   15.47 | insufficient |         0 |     nan     |        nan |
| chr22:20357407:A:G | ENSG00000287446 | none       |      | Y        |       1 |       3474 | 0.953 |   15.44 | insufficient |         0 |     nan     |        nan |
| chr22:18954277:C:G | ENSG00000280418 | none       | Y    |          |       1 |        110 | 0.999 |   15.41 | insufficient |         0 |     nan     |        nan |
| chr22:45102024:T:C | ENSG00000273243 | dELS       | Y    |          |       8 |     -54140 | 0.997 |   15.41 | insufficient |         0 |     nan     |        nan |
| chr22:24004423:C:T | ENSG00000235689 | CA-TF      | Y    |          |       4 |      -1148 | 0.939 |   15.34 | insufficient |         0 |     nan     |        nan |
| chr22:35004295:A:C | LINC02885       | none       | Y    |          |       1 |       1434 | 0.998 |   15.31 | insufficient |         0 |     nan     |        nan |
| chr22:30480905:T:C | SDC4P           | none       | Y    |          |       1 |      -1009 | 1     |   15.3  | insufficient |         0 |     nan     |        nan |
| chr22:37051815:G:A | TRIOBP          | none       |      |          |       1 |      40560 | 0.974 |   15.29 | insufficient |         0 |     nan     |        nan |
| chr22:48490892:G:A | ENSG00000281732 | none       |      |          |       1 |     -12974 | 0.994 |   15.22 | insufficient |         0 |     nan     |        nan |
| chr22:30490899:G:A | SDC4P           | CA-CTCF    | Y    |          |       2 |       8979 | 0.88  |   15.14 | insufficient |         1 |      -0.333 |          1 |
| chr22:37200293:A:G | C1QTNF6         | dELS       | Y    |          |       1 |        937 | 0.833 |   15.06 | insufficient |         0 |     nan     |        nan |
| chr22:21309286:C:G | SCARF2          | pELS       | Y    | Y        |       1 |        536 | 0.995 |   15.05 | insufficient |         0 |     nan     |        nan |
| chr22:21944102:G:A | PPM1F           | dELS       | Y    |          |       1 |      -4552 | 0.962 |   14.98 | insufficient |         0 |     nan     |        nan |
| chr22:31620173:G:A | LINC01521       | pELS       | Y    | Y        |       1 |        863 | 0.998 |   14.98 | insufficient |         0 |     nan     |        nan |
| chr22:23388183:T:C | ENSG00000290449 | none       |      | Y        |       1 |        617 | 0.846 |   14.96 | insufficient |         0 |     nan     |        nan |
| chr22:23975074:G:C | ENSG00000290199 | CA-CTCF    | Y    |          |       5 |        671 | 0.998 |   14.94 | insufficient |         1 |       1     |          1 |
| chr22:49725075:C:A | ENSG00000273192 | CA-H3K4me3 |      |          |       1 |     -17596 | 0.972 |   14.94 | insufficient |         0 |     nan     |        nan |