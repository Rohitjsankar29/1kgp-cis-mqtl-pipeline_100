#!/usr/bin/env python3
"""
05b_fdr_significant.py — FDR control on permutation p-values.
Tries tensorQTL's Storey q-values (needs R+qvalue); otherwise falls back to a
dependency-free Benjamini-Hochberg in numpy. Writes the FDR table and the list
of significant CpGs (input to stage 5c).
"""
import argparse, os, time
import numpy as np
import pandas as pd


def log(m):
    print(time.strftime("%H:%M:%S"), m, flush=True)


def bh_qvalues(p):
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--permutation", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--label", default="chr22")
    ap.add_argument("--fdr", type=float, default=0.05)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    df = pd.read_csv(a.permutation, sep="\t", index_col=0)
    pcol = "pval_beta" if "pval_beta" in df.columns else "pval_perm"
    valid = df[df[pcol].notna()].copy()
    log(f"{valid.shape[0]} CpGs with valid {pcol} (of {df.shape[0]})")

    try:
        from tensorqtl import post
        post.calculate_qvalues(valid, fdr=a.fdr)   # adds 'qval' in place
        used = "tensorqtl Storey q-values"
    except Exception as e:
        log(f"tensorqtl q-values unavailable ({e}); using BH fallback")
        valid["qval"] = bh_qvalues(valid[pcol].values)
        used = "Benjamini-Hochberg (numpy)"

    sig = valid[valid["qval"] < a.fdr]
    log(f"FDR method: {used} | significant CpGs at FDR<{a.fdr}: {sig.shape[0]}")

    fdr_path = f"{a.out_dir}/{a.label}.cis_qtl.fdr.txt.gz"
    sig_path = f"{a.out_dir}/{a.label}.significant_cpgs.txt"
    valid.to_csv(fdr_path, sep="\t")
    sig.index.to_series().to_csv(sig_path, index=False, header=False)
    log(f"Written: {fdr_path}")
    log(f"Written: {sig_path} ({sig.shape[0]} CpGs)")
    if sig.shape[0] == 0:
        log("No significant CpGs (expected at small n). Stage 5c will be a no-op.")


if __name__ == "__main__":
    main()
