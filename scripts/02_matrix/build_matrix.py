#!/usr/bin/env python3
"""
build_matrix.py — collapse per-sample modkit bedMethyl into a tensorQTL
phenotype BED for cis-mQTL mapping.

Input: a sample sheet TSV (no header):  <sample_id>\t<path_to_bedmethyl.gz>
For each sample it reads the 5mC ('m') rows at valid coverage >= --min-cov,
takes percent-modified (col 11) as the methylation fraction (beta), builds a
CpG x sample matrix, keeps CpGs covered in >= --min-frac of samples, imputes
the few remaining gaps with the per-CpG mean, optionally M-value transforms,
and writes a sorted tensorQTL BED (#chr start end phenotype_id <samples...>).

modkit bedMethyl columns used (0-based): 0=chrom 1=start 3=mod_code
10=valid_coverage 11=percent_modified
"""
import argparse, gzip, sys
import numpy as np
import pandas as pd


def load_sample(path, min_cov, mod="m"):
    """Return dict {start_pos: beta} for CpGs of one sample passing min_cov."""
    pos, beta = [], []
    with gzip.open(path, "rt") as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 12 or c[3] != mod:
                continue
            try:
                cov = int(c[9])
            except ValueError:
                continue
            if cov < min_cov:
                continue
            pos.append(int(c[1]))
            beta.append(float(c[10]) / 100.0)
    return pd.Series(beta, index=pos, dtype="float64")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sheet", required=True, help="TSV: sample_id <tab> bedmethyl_path")
    ap.add_argument("--out", required=True, help="output BED (uncompressed; bgzip after)")
    ap.add_argument("--chrom", default="chr22")
    ap.add_argument("--min-cov", type=int, default=5)
    ap.add_argument("--min-frac", type=float, default=0.90,
                    help="keep CpGs covered in >= this fraction of samples")
    ap.add_argument("--mod", default="m", help="modkit mod code: m=5mC, h=5hmC")
    ap.add_argument("--mvalue", action="store_true", help="logit (M-value) transform")
    a = ap.parse_args()

    sheet = pd.read_csv(a.sheet, sep="\t", header=None, names=["sid", "path"])
    cols = {}
    for _, r in sheet.iterrows():
        s = load_sample(r["path"], a.min_cov, a.mod)
        cols[r["sid"]] = s
        print(f"  {r['sid']}: {len(s):>7} CpGs", file=sys.stderr)

    mat = pd.DataFrame(cols)                      # index=pos, cols=sample
    n = mat.shape[1]
    keep = mat.notna().sum(axis=1) >= a.min_frac * n
    mat = mat.loc[keep].sort_index()
    print(f"CpGs covered in >= {a.min_frac:.0%} of {n} samples: "
          f"{mat.shape[0]} (of {len(keep)} union)", file=sys.stderr)
    if mat.shape[0] == 0:
        sys.exit("No CpGs passed the coverage threshold — lower --min-frac")

    # impute remaining gaps with per-CpG (row) mean, vectorised
    vals = mat.to_numpy(dtype="float64", copy=True)
    rmean = np.nanmean(vals, axis=1)
    nan_i, nan_j = np.where(np.isnan(vals))
    vals[nan_i, nan_j] = rmean[nan_i]

    if a.mvalue:
        eps = 1e-2
        vals = np.log2((vals + eps) / (1.0 - vals + eps))

    starts = mat.index.astype(int).to_numpy()
    bed = pd.DataFrame(vals, columns=mat.columns)
    bed.insert(0, "phenotype_id", [f"{a.chrom}_{p}" for p in starts])
    bed.insert(0, "end", starts + 1)
    bed.insert(0, "start", starts)
    bed.insert(0, "#chr", a.chrom)
    bed = bed.sort_values(["#chr", "start"]).reset_index(drop=True)
    bed.to_csv(a.out, sep="\t", index=False, float_format="%.4f")
    print(f"Wrote {bed.shape[0]} CpGs x {n} samples -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
