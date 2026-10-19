#!/usr/bin/env python3
"""Drop leaves whose share of gaps in the concatenated alignment is above a threshold.

The gap fraction is computed the way horizontal.py builds the concatenation: a genome
absent from a block contributes that block's full length as gaps. Dropped genomes are
removed from every block, so the blocks still add up to the tree alignment.
"""
import argparse
import os
import sys


def read_fasta(path):
    seqs, name = {}, None
    with open(path) as fin:
        for line in fin:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name is not None:
                seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def gap_fractions(blocks):
    """blocks: {path: {name: seq}} -> {name: fraction of gap columns in the concatenation}"""
    lengths = {p: len(next(iter(s.values()))) if s else 0 for p, s in blocks.items()}
    total = sum(lengths.values())
    names = sorted({n for s in blocks.values() for n in s})
    frac = {}
    for n in names:
        gaps = 0
        for p, seqs in blocks.items():
            seq = seqs.get(n)
            gaps += lengths[p] if seq is None else seq.count("-") + seq.count(".")
        frac[n] = gaps / total if total else 1.0
    return frac


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--clade", required=True)
    ap.add_argument("--max-gap-fraction", type=float, default=None,
                    help="drop leaves above this fraction; unset = report only")
    ap.add_argument("--outdir", default="filtered")
    ap.add_argument("alignments", nargs="+")
    args = ap.parse_args()

    blocks = {p: read_fasta(p) for p in sorted(args.alignments)}
    frac = gap_fractions(blocks)
    thr = args.max_gap_fraction
    drop = {n for n, f in frac.items() if thr is not None and f > thr}

    os.makedirs(args.outdir, exist_ok=True)
    for path, seqs in blocks.items():
        kept = {n: s for n, s in seqs.items() if n not in drop}
        if seqs and not kept:
            sys.exit(f"dropping {len(drop)} leaves would empty {path}; raise --max-gap-fraction")
        with open(os.path.join(args.outdir, os.path.basename(path)), "w") as fout:
            for n, s in kept.items():
                fout.write(f">{n}\n{s}\n")

    with open(f"leaf_gaps_{args.clade}.tsv", "w") as fout:
        fout.write("accession\tgap_fraction\texcluded\n")
        for n in sorted(frac):
            fout.write(f"{n}\t{frac[n]:.4f}\t{int(n in drop)}\n")
    print(f"{args.clade}: {len(frac)} leaves, {len(drop)} above {thr}", file=sys.stderr)


if __name__ == "__main__":
    main()
