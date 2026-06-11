# ENVIRONMENT — toolchain on Gadi

Tools are provided via Gadi modules and three micromamba environments under
`/scratch/cy94/rs4477/micromamba/envs/`.

## Gadi modules
- `samtools/1.22`  — samtools (streaming chr22)
- `bcftools/1.21`  — VCF subset + sample rename

## micromamba environments
| env | key packages | used for |
|-----|--------------|----------|
| `mqtl_modkit`     | modkit, bgzip, tabix         | methylation pileup, bgzip/tabix |
| `mqtl_tensorqtl`  | tensorqtl 1.0.10, torch, pandas, numpy, scipy | QTL mapping, matrix, covariates |
| `plink`           | plink2                       | genotype PLINK conversion, PCA |

Call binaries by absolute path (no activation needed), e.g.
`/scratch/cy94/rs4477/micromamba/envs/mqtl_tensorqtl/bin/python3`.

## GPU / torch — IMPORTANT
The permutation runs on a Tesla **V100-SXM2-32GB**, which is **compute
capability 7.0 (Volta)**. The torch build *must* include `sm_70` kernels:

```bash
# working build for the V100 (runs fine on the CUDA 13.0 driver — drivers are
# backward compatible):
/scratch/cy94/rs4477/micromamba/envs/mqtl_tensorqtl/bin/python3 \
  -m pip install "torch==2.4.1" --index-url https://download.pytorch.org/whl/cu121
```

A `+cu130` / torch ≥ 2.6 build was compiled only for CC ≥ 7.5 and fails on the
V100 with `CUDA error: no kernel image is available for execution on the device`,
even though `torch.cuda.is_available()` returns True.

## plink2 install (if missing)
```bash
micromamba create -y -p /scratch/cy94/rs4477/micromamba/envs/plink \
  -c bioconda -c conda-forge plink2
```

## Reference / data
- hg38 fasta: `/scratch/cy94/rs4477/reference/hg38.fa`
- 1000G chr22 high-coverage phased panel under `/g/data/xl04/rs4477/...`
- ONT modBAMs: 1000g-ont S3 bucket (public https), streamed per-chromosome.
