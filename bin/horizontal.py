#!/usr/bin/env python3
import os
import logging
from operator import itemgetter
import argparse
from collections import defaultdict

# now we need to concatenate horizontally


def readFasta(file):
    seqs = {}
    with open(file) as f:
        for line in f:
            if line.startswith(">"):
                name = line.strip()
                seqs[name] = ""
            else:
                seqs[name] += line.strip()
    return seqs


def horizontalConcat(files, profiles, output):
    seqs = {}
    allnames = set()
    for profile, f in zip(profiles, files):
        # print("Reding {} {}".format(profile,f))
        seqs[profile] = readFasta(f)
        for name, seq in seqs[profile].items():
            allnames.add(name)

    # length matter
    # as I will need to pad single seqs
    lengths = defaultdict(int)
    for profile, seq in seqs.items():
        for k, s in seq.items():
            lengths[profile] = len(s)
            break

    # add missing seqences as strings of spaces
    for profile, seq in seqs.items():
        for name in allnames:
            if name not in seqs[profile].keys():
                seqs[profile][name] = "".join(lengths[profile] * ["-"])

    mergedSeqs = {}
    # sorted: a set's order changes between runs (string hashing is randomised)
    for name in sorted(allnames):
        s = ""
        for p in profiles:
            s += seqs[p][name]
        mergedSeqs[name] = s

    # write output:
    with open(output, "w") as f:
        for name, seq in mergedSeqs.items():
            seq = seq.replace(".", "-")
            f.write("{}\n{}\n".format(name, seq))

    return output


parser = argparse.ArgumentParser(description="Get the top entry from hmmsearch results")
parser.add_argument("outfile", type=str, help="output file")
parser.add_argument("faa", type=str, nargs="+", help="fasta file")
args = parser.parse_args()

# define logging
logLevel = logging.INFO
logging.basicConfig(format="%(asctime)s %(message)s", datefmt="%m/%d/%Y %H:%M:%S: ", level=logLevel)


files = args.faa
profiles = [os.path.basename(x) for x in args.faa]
# make sure it is orted here, as we will need to preserve the order
profiles.sort()
horizontalConcat(files, profiles, args.outfile)
