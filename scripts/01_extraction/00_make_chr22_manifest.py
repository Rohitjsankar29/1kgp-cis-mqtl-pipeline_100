#!/usr/bin/env python3
"""
00_make_chr22_manifest.py — manifest of NEW 1KGP ONT samples for chr22 extraction.

Lists the 1000g-ont S3 bucket, EXCLUDES samples already done, takes the next N,
and writes a TSV that scripts/prototype/01_extract_chr_modkit.prototype.sh can
read directly. The prototype uses:
    col1 = sample_id
    col2 = bam_url     (must be https:// — confirmed working with samtools)
    col9 = modifications  ("5mC" or "5mC+5hmC" -> drives the per-sample mod flag)

Run on the LOGIN NODE (needs internet). Stdlib only — any python3 works.

Usage:
    python3 00_make_chr22_manifest.py \
        --out /scratch/cy94/$USER/1kgp-cis-mqtl/config/chr22_manifest_95.tsv \
        --n 95
"""
import argparse, re, sys, urllib.parse, urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

S3_HTTP = "https://s3.amazonaws.com/1000g-ont"
BAM_PREFIX = "PROCESSED_DATA/ALIGNED_TO_HG38/MINIMAP2_ALIGNED_BAMS/"
NS = "http://s3.amazonaws.com/doc/2006-03-01/"

# The five samples you have already processed — never re-pick these.
DONE_DEFAULT = "HG02470,HG02479,HG03027,HG03045,HG03079"


def list_s3(prefix):
    out, token = [], None
    while True:
        url = f"{S3_HTTP}?list-type=2&prefix={urllib.parse.quote(prefix)}&max-keys=1000"
        if token:
            url += f"&continuation-token={urllib.parse.quote(token)}"
        with urllib.request.urlopen(url) as r:
            root = ET.parse(r).getroot()
        for c in root.findall(f"{{{NS}}}Contents"):
            key = c.find(f"{{{NS}}}Key").text
            size = int(c.find(f"{{{NS}}}Size").text)
            out.append((key, size))
        t = root.find(f"{{{NS}}}IsTruncated")
        if t is not None and t.text == "true":
            token = root.find(f"{{{NS}}}NextContinuationToken").text
        else:
            return out


def parse_name(fname):
    m = re.match(r"^([A-Z]+\d+)-ONT-", fname)
    sid = m.group(1) if m else fname.split("-ONT-")[0]
    pore = "R10" if "-R10-" in fname else "R9"
    mods = "5mC+5hmC" if ("5hmC" in fname or "5hmCG" in fname) else "5mC"
    if "guppy657" in fname:
        bc = "guppy657"
    elif "guppy" in fname:
        bc = "guppy"
    else:
        dm = re.search(r"dorado(\d+)", fname)
        bc = "dorado" + dm.group(1) if dm else "unknown"
    return sid, pore, mods, bc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--n", type=int, default=95)
    ap.add_argument("--exclude", default=DONE_DEFAULT,
                    help="comma-separated sample IDs to skip (already done)")
    args = ap.parse_args()
    done = set(x.strip() for x in args.exclude.split(",") if x.strip())

    objs = list_s3(BAM_PREFIX)
    bams = [(k, s) for k, s in objs if k.endswith(".phased.bam")]
    rows, seen = [], set()
    for key, size in sorted(bams):
        fn = Path(key).name
        sid, pore, mods, bc = parse_name(fn)
        if sid in done or sid in seen:
            continue
        seen.add(sid)
        rows.append([
            sid,                          # 1 sample_id
            f"{S3_HTTP}/{key}",           # 2 bam_url (https)  <- samtools reads this
            f"{S3_HTTP}/{key}.bai",       # 3 bai_url
            fn,                           # 4 bam_filename
            round(size / 1e9, 2),         # 5 bam_size_gb
            pore,                         # 6 pore
            1 if pore == "R10" else 0,    # 7 basecaller_model (R9=0,R10=1)
            bc,                           # 8 basecaller
            mods,                         # 9 modifications  <- prototype mod flag
            str("5hmC" in mods).upper(),  # 10 has_5hmc
        ])
        if len(rows) >= args.n:
            break

    if not rows:
        sys.exit("No samples selected — check exclude list / bucket")

    hdr = ["sample_id", "bam_url", "bai_url", "bam_filename", "bam_size_gb",
           "pore", "basecaller_model", "basecaller", "modifications", "has_5hmc"]
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as f:
        f.write("\t".join(hdr) + "\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")

    r9 = sum(1 for r in rows if r[5] == "R9")
    r10 = len(rows) - r9
    hmc = sum(1 for r in rows if r[9] == "TRUE")
    print(f"Wrote {len(rows)} samples -> {args.out}", file=sys.stderr)
    print(f"Pore: R9={r9} R10={r10} | with 5hmC: {hmc}", file=sys.stderr)
    print(f"First 3: {[r[0] for r in rows[:3]]}", file=sys.stderr)


if __name__ == "__main__":
    main()
