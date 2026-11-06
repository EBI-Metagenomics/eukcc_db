#!/usr/bin/env python3
"""Drop leaves whose share of gaps in the concatenated alignment is above a threshold.

The gap fraction is computed the way horizontal.py builds the concatenation: a genome
absent from a block contributes that block's full length as gaps. Dropped genomes are
removed from every block, so the blocks still add up to the tree alignment.
"""
import argparse
import os
import sys

from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment


def read_block(path):
    """One alignment block: {name: record}; MultipleSeqAlignment checks the equal lengths"""
    alignment = MultipleSeqAlignment(SeqIO.parse(path, "fasta"))
    return {record.id: record for record in alignment}


def gap_fractions(blocks):
    """blocks: {path: {name: record}} -> {name: fraction of gap columns in the concatenation}"""
    lengths = {p: len(next(iter(s.values()))) if s else 0 for p, s in blocks.items()}
    total = sum(lengths.values())
    names = sorted({n for s in blocks.values() for n in s})
    frac = {}
    for n in names:
        gaps = 0
        for p, records in blocks.items():
            record = records.get(n)
            gaps += lengths[p] if record is None else record.seq.count("-") + record.seq.count(".")
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

    blocks = {p: read_block(p) for p in sorted(args.alignments)}
    frac = gap_fractions(blocks)
    thr = args.max_gap_fraction
    drop = {n for n, f in frac.items() if thr is not None and f > thr}

    os.makedirs(args.outdir, exist_ok=True)
    for path, records in blocks.items():
        kept = [r for n, r in records.items() if n not in drop]
        if records and not kept:
            sys.exit(f"dropping {len(drop)} leaves would empty {path}; raise --max-gap-fraction")
        for r in kept:
            r.description = ""  # header = the name only
        # fasta-2line: one line per sequence, as the blocks were written before
        SeqIO.write(kept, os.path.join(args.outdir, os.path.basename(path)), "fasta-2line")

    with open(f"leaf_gaps_{args.clade}.tsv", "w") as fout:
        fout.write("accession\tgap_fraction\texcluded\n")
        for n in sorted(frac):
            fout.write(f"{n}\t{frac[n]:.4f}\t{int(n in drop)}\n")
    print(f"{args.clade}: {len(frac)} leaves, {len(drop)} above {thr}", file=sys.stderr)


if __name__ == "__main__":
    main()
