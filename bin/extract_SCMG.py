#!/usr/bin/env python
# coding: utf-8

import os
import logging
from collections import defaultdict
import argparse
from random import seed
import gzip

seed(42)


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


def read_hmmer(path):
    data = read_tsv(path, sep=None, header=False)
    filename = os.path.basename(path)
    name = filename.replace("hmmerGA_", "").replace(".tbl", "")
    logging.debug("Reading file of acc {}".format(name))

    # reduce to what we need:
    perProfile = defaultdict(list)
    for row in data:
        profile = row[2].split(".", 1)[0]
        perProfile[profile].append(name)

    return perProfile


def main(**kwargs):
    # load all tables into memory
    logging.info("Reading in all the data")
    with gzip.open("scmgs.gz", "wb") as fout:
        for tbl in kwargs["tbls"]:
            for profile, acc in read_hmmer(tbl).items():
                if len(acc) == 1:
                    fout.write("{},{}\n".format(acc[0], profile).encode())
    logging.info("Done")
    return


if __name__ == "__main__":
    # define arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("tbls", help="run quietly", nargs="+")
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
