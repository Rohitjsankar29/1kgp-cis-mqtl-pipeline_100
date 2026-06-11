#!/usr/bin/env python3
"""
build_covariates.py — assemble tensorQTL covariates from data only:
  - genotype PCs   (from plink2 --pca .eigenvec)  -> population structure
  - methylation PCs (PCA of the phenotype matrix)  -> hidden technical structure
Output: covariates.tsv in tensorQTL layout (covariates as rows, samples as cols).
"""
import argparse, sys
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", required=True, help="methylation_Mval.bed.gz")
    ap.add_argument("--geno-pca", required=True, help="plink2 .eigenvec")
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-meth-pc", type=int, default=10)
    a = ap.parse_args()

    # --- genotype PCs ---
    gp = pd.read_csv(a.geno_pca, sep=r"\s+")
    iid = "IID" if "IID" in gp.columns else gp.columns[1]
    gp = gp.set_index(iid)
    pcs = [c for c in gp.columns if c.upper().startswith("PC")]
    geno = gp[pcs].copy()
    geno.columns = [f"genoPC{i+1}" for i in range(len(pcs))]
    print(f"geno PCs: {geno.shape}", file=sys.stderr)

    # --- methylation PCs (eigU of the 100x100 sample Gram matrix) ---
    df = pd.read_csv(a.matrix, sep="\t")
    samples = list(df.columns[4:])
    X = df[samples].to_numpy(dtype=float).T          # samples x CpGs
    X = (X - X.mean(0)) / (X.std(0) + 1e-8)          # standardise CpGs
    G = X @ X.T                                       # samples x samples
    val, vec = np.linalg.eigh(G)
    order = np.argsort(val)[::-1]
    k = min(a.n_meth_pc, vec.shape[1])
    meth = pd.DataFrame(vec[:, order[:k]], index=samples,
                        columns=[f"methPC{i+1}" for i in range(k)])
    print(f"meth PCs: {meth.shape}", file=sys.stderr)

    # --- merge, keep matrix sample order ---
    cov = meth.join(geno, how="left").loc[samples]
    if cov.isna().any().any():
        miss = cov.index[cov.isna().any(axis=1)].tolist()
        print(f"WARN: missing geno PCs for {miss} -> filled 0", file=sys.stderr)
        cov = cov.fillna(0.0)

    out = cov.T                                       # covariates x samples
    out.index.name = "id"
    out.to_csv(a.out, sep="\t")
    print(f"covariates: {out.shape[0]} x {out.shape[1]} -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
