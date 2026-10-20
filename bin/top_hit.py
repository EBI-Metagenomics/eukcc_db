#!/usr/bin/env python3
import os
import logging
from pyfaidx import Fasta
from operator import itemgetter
import argparse


def read_hmmer(p, profile):
    if not os.path.exists(p):
        logging.error("File {} does not exist".format(p))
        exit(1)
    results = []
    with open(p) as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            fields = line.strip().split()
            if fields[2].startswith(profile):
                key = fields[0]
                bitscore = float(fields[5])
                results.append([key, bitscore])
    if len(results) > 0:
        results = sorted(results, key=itemgetter(1), reverse=True)

    return results


def return_top(hits, faa_p):
    if len(hits) == 0:
        return None
    faa = Fasta(faa_p)
    key = hits[0][0]
    if key in faa.keys():
        return faa[key]


parser = argparse.ArgumentParser(description="Get the top entry from hmmsearch results")
parser.add_argument("--faa", type=str, help="fasta file")
parser.add_argument("--hits", type=str, help="hits file")
parser.add_argument("--profile", type=str, help="profile")
parser.add_argument("--name", type=str, help="name file", default="NAME")
parser.add_argument("outfile", type=str, help="output file")
args = parser.parse_args()

# define logging
logLevel = logging.INFO
logging.basicConfig(format="%(asctime)s %(message)s", datefmt="%m/%d/%Y %H:%M:%S: ", level=logLevel)


hmmresults = read_hmmer(args.hits, args.profile)
top = return_top(hmmresults, args.faa)


with open(args.outfile, "w") as fout:
    if top is not None:
        fout.write(">{}\n".format(args.name))
        fout.write("{}\n".format(top))
