#!/usr/bin/env python3
"""asm_readlevel.py — fragment-level, region-based allele-specific methylation,
following the Atlas method (Rosenski et al., Nat Commun 2025).

For each candidate region:
  1. take every read (fragment) covering >= MIN_CPG CpGs in the region
  2. classify the fragment by its average methylation:
        U (hypo)  if mean <= 0.35
        M (hyper)  if mean >= 0.65
        X (mixed)  otherwise
  3. split fragments by allele (from the HP tag at a het SNP: HP1 vs HP2)
  4. build the allele x {U,M} contingency table and run Fisher's exact test
     (X fragments are excluded from the 2x2, as in the Atlas)
  5. BH-FDR across regions; flag ASM at q < 0.01

Input --calls: a read-level methylation table (e.g. from `modkit extract calls`)
with columns: read_id, chrom, ref_position, call, hp
  - `call` is methylated if in --meth-codes (default: m,C+m,1,methylated)
  - `hp` is 1 or 2 (the haplotype = allele at the het SNP)
Input --regions: BED (chrom, start, end, name) of candidate regions.
"""
import argparse
import gzip

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


def opener(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", required=True)
    ap.add_argument("--regions", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-cpg", type=int, default=3)      # >=3 CpGs per fragment
    ap.add_argument("--u-max", type=float, default=0.35)   # U threshold
    ap.add_argument("--m-min", type=float, default=0.65)   # M threshold
    ap.add_argument("--min-frag", type=int, default=5)     # >=5 frags per allele
    ap.add_argument("--meth-codes", default="m,C+m,1,methylated")
    a = ap.parse_args()
    meth = set(a.meth_codes.split(","))

    calls = pd.read_csv(a.calls, sep="\t",
                        dtype={"read_id": str, "chrom": str})
    calls["is_meth"] = calls["call"].astype(str).isin(meth).astype(int)
    calls["ref_position"] = pd.to_numeric(calls["ref_position"], errors="coerce")
    calls["hp"] = pd.to_numeric(calls["hp"], errors="coerce")

    regions = pd.read_csv(a.regions, sep="\t", header=None,
                          names=["chrom", "start", "end", "name"],
                          dtype={"chrom": str})

    rows = []
    for r in regions.itertuples():
        sub = calls[(calls["chrom"] == r.chrom) &
                    (calls["ref_position"] >= r.start) &
                    (calls["ref_position"] < r.end) &
                    (calls["hp"].isin([1, 2]))]
        if sub.empty:
            continue
        # per-fragment average methylation over its CpGs in the region
        frag = sub.groupby("read_id").agg(hp=("hp", "first"),
                                          ncpg=("is_meth", "size"),
                                          meth=("is_meth", "mean"))
        frag = frag[frag["ncpg"] >= a.min_cpg]
        if frag.empty:
            continue
        frag["cls"] = np.where(frag["meth"] <= a.u_max, "U",
                       np.where(frag["meth"] >= a.m_min, "M", "X"))
        # contingency: allele (HP1/HP2) x {U, M}
        u1 = int(((frag.hp == 1) & (frag.cls == "U")).sum())
        m1 = int(((frag.hp == 1) & (frag.cls == "M")).sum())
        u2 = int(((frag.hp == 2) & (frag.cls == "U")).sum())
        m2 = int(((frag.hp == 2) & (frag.cls == "M")).sum())
        n1, n2 = u1 + m1, u2 + m2
        if n1 < a.min_frag or n2 < a.min_frag:
            continue
        _, p = fisher_exact([[u1, m1], [u2, m2]])
        # direction: which haplotype is more methylated
        meth1 = m1 / n1 if n1 else np.nan
        meth2 = m2 / n2 if n2 else np.nan
        rows.append((r.chrom, r.start, r.end, r.name, u1, m1, u2, m2,
                     round(meth1, 3), round(meth2, 3),
                     round(meth1 - meth2, 3), p))

    df = pd.DataFrame(rows, columns=["chrom", "start", "end", "name",
                                     "U_hp1", "M_hp1", "U_hp2", "M_hp2",
                                     "meth_hp1", "meth_hp2", "delta", "p"])
    if len(df):
        m = len(df)
        rank = df["p"].rank(method="first")
        df["q"] = (df["p"] * m / rank).clip(upper=1.0)
        df = df.sort_values("p").reset_index(drop=True)
    df.to_csv(a.out, sep="\t", index=False)
    n_asm = int((df["q"] < 0.01).sum()) if len(df) else 0
    print(f"regions tested: {len(df):,}   ASM (q<0.01): {n_asm:,}   -> {a.out}")
    if len(df):
        print(df.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
