#!/usr/bin/env python
# coding: utf-8

import os
import logging
import argparse


def read_tsv(path, header=True, sep="\t", comment="#"):
    """
    simple csv reader
    """
    head = None
    content = []

    with open(path) as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            line = line.strip()
            if header and head is None:
                if sep is None:
                    head = line.split()
                else:
                    head = line.split(sep)
                continue
            else:
                if sep is None:
                    elements = line.split(sep)
                else:
                    elements = line.split()
            if head is not None and header:
                content.append(dict(zip(head, elements)))
            else:
                content.append(elements)
    return content


def add_GA(profile, basedir, score, fout):
    inputfile = os.path.join(basedir, profile, "hmmer.hmm")
    with open(inputfile) as fin:
        for line in fin:
            if line.startswith("DATE"):
                fout.write("GA    {} {}\n".format(score, score))
            fout.write(line)
    return


def main(**kwargs):
    data = read_tsv(kwargs["info"])
    with open(kwargs["outfile"], "w") as fout:
        for row in data:
            add_GA(row["profile"], kwargs["panther"], row["cutoff"], fout)

    return


if __name__ == "__main__":
    # define arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("outfile", help="scores table", type=str)
    parser.add_argument("--panther", help="panther dir", type=str)
    parser.add_argument("--info", help="scores table", type=str)
    parser.add_argument("--quiet", help="run quietly", action="store_true", default=False)
    parser.add_argument("--debug", help="run quietly", action="store_true", default=False)
    args = parser.parse_args()

    # define logging
    logLevel = logging.INFO
    if args.quiet:
        logLevel = logging.WARNING
    elif args.debug:
        logLevel = logging.DEBUG
    logging.basicConfig(format="%(asctime)s %(message)s", datefmt="%d/%m/%Y %H:%M:%S: ", level=logLevel)

    main(**vars(args))
