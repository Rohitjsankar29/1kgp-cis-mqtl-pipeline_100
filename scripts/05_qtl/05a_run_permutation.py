#!/usr/bin/env python3
"""
05a_run_permutation.py — cis-mQTL permutation pass (tensorQTL map_cis).
One empirical p-value per CpG. GPU-aware: uses CUDA automatically when torch
sees a supported device. The V100 needs a torch built with sm_70 (e.g. cu121).
"""
import argparse, sys, time
import numpy as np
import pandas as pd


def log(msg):
    print(time.strftime("%H:%M:%S"), msg, flush=True)


def inverse_normal_transform(df):
    """Rank-based inverse-normal transform per phenotype (row); samples in cols."""
    from scipy.stats import norm
    ranks = df.rank(axis=1, method="average")
    n = df.shape[1]
    return pd.DataFrame(norm.ppf((ranks - 0.5) / n),
                        index=df.index, columns=df.columns)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plink-prefix", required=True)
    ap.add_argument("--phenotype-bed", required=True)
    ap.add_argument("--covariates", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--cis-window", type=int, default=1000000)
    ap.add_argument("--maf-threshold", type=float, default=0.05)
    ap.add_argument("--n-permutations", type=int, default=10000)
    ap.add_argument("--inverse-normal", action="store_true")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()

    import torch
    torch.set_num_threads(a.threads)
    import tensorqtl
    from tensorqtl import genotypeio, cis
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    log(f"tensorqtl {tensorqtl.__version__} | torch {torch.__version__} | device={dev}")
    if dev == "cpu":
        log("Running on CPU — permutations are slow. Use the gpuvolta job for a GPU.")

    phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(a.phenotype_bed)
    log(f"Phenotypes: {phenotype_df.shape[0]:,} CpGs x {phenotype_df.shape[1]} samples")
    if a.inverse_normal:
        phenotype_df = inverse_normal_transform(phenotype_df)
        log("Applied inverse-normal transform to phenotypes")

    covariates_df = pd.read_csv(a.covariates, sep="\t", index_col=0).T
    log(f"Covariates: {covariates_df.shape[1]} covariates x "
        f"{covariates_df.shape[0]} samples ({list(covariates_df.columns)})")

    pr = genotypeio.PlinkReader(a.plink_prefix)
    genotype_df = pr.load_genotypes()
    variant_df = pr.bim.set_index("snp")[["chrom", "pos"]]
    log(f"Genotypes: {genotype_df.shape[0]:,} variants x {genotype_df.shape[1]} samples")

    # chrom-mismatch guard (the #1 cause of "0 variants in cis-window")
    gchr = set(variant_df["chrom"].astype(str).unique())
    pchr = set(phenotype_pos_df["chr"].astype(str).unique())
    if gchr.isdisjoint(pchr):
        log(f"Chrom names disagree! genotype={sorted(gchr)} phenotype={sorted(pchr)}")
        sys.exit(2)

    n = phenotype_df.shape[1]
    log(f"n={n} samples | {covariates_df.shape[1]} covariates | "
        f"residual df ~ {n - covariates_df.shape[1] - 2}")
    if n < 30:
        log("n<30: permutation p-values are weak. Treat results as VALIDATION.")
    log(f"map_cis: window=+/-{a.cis_window:,} bp | maf>={a.maf_threshold} | "
        f"nperm={a.n_permutations} | seed={a.seed}")

    cis_df = cis.map_cis(
        genotype_df, variant_df, phenotype_df, phenotype_pos_df,
        covariates_df=covariates_df, maf_threshold=a.maf_threshold,
        nperm=a.n_permutations, window=a.cis_window, seed=a.seed,
    )
    log(f"map_cis done: {cis_df.shape[0]} CpGs with a tested cis-window")
    cis_df.to_csv(a.out, sep="\t")
    log(f"Written: {a.out}")


if __name__ == "__main__":
    main()
