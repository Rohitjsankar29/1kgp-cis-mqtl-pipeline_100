#!/usr/bin/env python3
"""Stage 09: integrate fine-mapping + functional + SV evidence into a single
ranked candidate-causal-variant list.

    score = w_pip * PIP  +  w_func * functional_weight  +  w_sv * sv_evidence

Weights are explicit and tunable so the ranking is transparent and reproducible.
"""
import argparse, os, time
import pandas as pd
import numpy as np


def log(m):
    print(time.strftime("%H:%M:%S"), m, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--susie", required=True, help="finemap output (phenotype_id, variant_id, pip, cs_id)")
    ap.add_argument("--functional", help="functional_annotate output")
    ap.add_argument("--sv", help="sv_annotate output")
    ap.add_argument("--out", required=True)
    ap.add_argument("--w-pip", type=float, default=1.0)
    ap.add_argument("--w-func", type=float, default=0.5)
    ap.add_argument("--w-sv", type=float, default=0.5)
    ap.add_argument("--cs-only", action="store_true", help="keep only variants in a credible set")
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)

    su = pd.read_csv(a.susie, sep="\t")
    if su.empty:
        log("empty fine-mapping input -- nothing to prioritise.")
        su.to_csv(a.out, sep="\t", index=False)
        return
    if a.cs_only and "cs_id" in su.columns:
        su = su[su["cs_id"].notna()]

    func_w = {}
    if a.functional and os.path.exists(a.functional):
        fa = pd.read_csv(a.functional, sep="\t")
        for _, r in fa.iterrows():
            w = 0.0
            if str(r.get("ccre", "none")) not in ("none", "", "nan"):
                w = max(w, 1.0)
            try:
                if abs(float(r.get("dist_to_tss", 1e9))) < 2000:
                    w = max(w, 1.0)
            except (TypeError, ValueError):
                pass
            if str(r.get("cpg_island", "no")) == "yes":
                w = max(w, 0.5)
            func_w[r["cpg"]] = w

    sv_w = {}
    if a.sv and os.path.exists(a.sv):
        sv = pd.read_csv(a.sv, sep="\t")
        if not sv.empty and "sv_implicated" in sv.columns:
            sv_w = sv.groupby("phenotype_id")["sv_implicated"].any().astype(float).to_dict()

    su["pip"] = pd.to_numeric(su["pip"], errors="coerce").fillna(0.0)
    su["func"] = su["phenotype_id"].map(func_w).fillna(0.0)
    su["sv"] = su["phenotype_id"].map(sv_w).fillna(0.0)
    su["score"] = a.w_pip * su["pip"] + a.w_func * su["func"] + a.w_sv * su["sv"]
    su = su.sort_values("score", ascending=False)
    su.to_csv(a.out, sep="\t", index=False)
    log(f"Wrote {a.out}: {su.shape[0]} variants ranked"
        + (f"; top score {su['score'].max():.3f}" if not su.empty else ""))


if __name__ == "__main__":
    main()
