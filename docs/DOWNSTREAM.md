# DOWNSTREAM — fine-mapping, SV integration, functional annotation, prioritisation

Everything here takes the **FDR-significant CpGs** from stage 05 and turns each
cis-mQTL into a ranked, annotated, candidate-causal-variant call. The structural
(SV) integration is the novel part of the framework — most mQTL pipelines stop
at SNVs.

```
significant CpGs (05b)  +  nominal stats (05c)
        │
        ▼
06  SuSiE fine-mapping        per-CpG credible sets + PIP per variant
        │
        ▼
07  SV integration  ***NOVEL***
        │  do structural variants explain the signal?
        │  (a) direct test: methylation ~ SV dosage
        │  (b) LD: is the lead SNV tagging a cis SV?
        ▼
08  functional annotation     CpG + variant context: gene/promoter/CpG-island,
        │                      ENCODE cCRE / chromHMM, distance to TSS
        ▼
09  prioritisation score      PIP × functional weight × SV evidence  → ranked list
```

## 06 — fine-mapping (SuSiE)
SuSiE assigns each cis variant a posterior inclusion probability (PIP) and groups
likely-causal variants into 95% credible sets, resolving the lead-SNP-vs-LD-block
ambiguity that a single top hit can't.

Dependency: `susieR` + `rpy2` in the tensorQTL env (tensorQTL calls SuSiE through R):
```bash
micromamba install -p /scratch/cy94/rs4477/micromamba/envs/mqtl_tensorqtl \
  -c conda-forge r-base r-susier rpy2
```
Run: `scripts/06_finemap/run_finemap.pbs` → `chr22.susie.txt.gz`
(columns: phenotype_id, variant_id, pip, cs_id). **Verify the SuSiE return-object
keys** in `finemap_susie.py` against your installed tensorQTL version — the parser
is defensive but the field names have changed across releases.

## 07 — SV integration (novel)
The 1000G high-coverage panel carries SVs in the same VCF (`...SNV_INDEL_SV...`),
but stage 03 deliberately kept only biallelic SNVs. Here we bring the SVs back in
**only around significant CpGs** and ask whether an SV — not a SNV — is driving
the methylation signal.

1. `extract_sv_genotypes.pbs` — pull chr22 SV records for the 100 samples from the
   panel (symbolic-ALT / `SVTYPE` records), reheader NA→GM, dump a dosage matrix.
2. `sv_annotate.py` — for each significant CpG:
   - **direct association**: regress covariate-residualised methylation on each cis
     SV's dosage → SV-mQTL effect and p-value;
   - **LD tagging**: r² between the lead SNV and each cis SV;
   - **flag** CpGs where an SV is directly associated *or* in high LD (r²≥0.8) with
     the lead SNV → candidate SV-driven mQTLs.

   Output: `chr22.sv_annotation.txt.gz` (CpG, SV, svtype, distance, sv_p, r2_lead, flag).

This is the defensible first version. The fuller extension — putting SVs *into* the
SuSiE genotype matrix so they compete with SNVs for PIP — is noted at the bottom of
`sv_annotate.py` as the next step (it needs SVs recoded to a plink-friendly form).

## 08 — functional annotation
Annotate each significant CpG and its lead/credible variants with regulatory context
using `bedtools intersect` against chr22 subsets of standard tracks:
- GENCODE genes/promoters (nearest gene, distance to TSS),
- UCSC CpG islands,
- ENCODE cCREs (candidate cis-regulatory elements) / chromHMM states.

`get_annotations.sh` downloads + subsets the tracks to chr22; `functional_annotate.py`
produces `chr22.functional.txt.gz`. **Decide the track set with Sam** — for the
PROPHECY application you may want blood/immune-cell-specific regulatory maps.

## 09 — prioritisation score
Combine the three evidence streams into one ranked table of candidate causal variants:

```
score = w_pip * PIP
      + w_func * functional_weight     # e.g. in a cCRE/promoter near the CpG
      + w_sv   * sv_evidence           # direct SV assoc or high-LD SV tag
```

Weights are explicit and tunable (`prioritise.py --w-pip --w-func --w-sv`); the point
is a transparent, reproducible ranking, not a black box. Output:
`chr22.prioritised.txt.gz`, sorted by score.

## Status / decisions for Sam
- Fine-mapping needs the susieR/rpy2 install above.
- SV integration thresholds (LD r², SV-assoc p) are placeholders — set with Sam.
- Functional tracks: generic ENCODE now; tissue-specific for PROPHECY.
- Prioritisation weights: start equal, calibrate against known mQTLs.
- All of this is chr22-only; the same scripts scale per-chromosome genome-wide.
