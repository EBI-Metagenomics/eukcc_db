#!/usr/bin/env python
# coding: utf-8

import logging
import argparse
from ete3 import Tree


def main(**kwargs):
    # load all tables into memory
    logging.info("Reading the tree")
    t = Tree(kwargs["tree"])
    R = t.get_midpoint_outgroup()
    t.set_outgroup(R)
    t.write(format=1, outfile=kwargs["outfile"])

    logging.info("Done")
    return


if __name__ == "__main__":
    # define arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("tree", help="run quietly", type=str)
    parser.add_argument("outfile", help="run quietly", type=str)
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
