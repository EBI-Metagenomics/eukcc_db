#!/usr/bin/env python
# coding: utf-8
import os
import logging
import json
import gzip
import argparse
import csv


def lcsv(f):
    with open(f) as fin:
        csvr = csv.DictReader(fin, delimiter="\t")
        for row in csvr:
            return row
            break


def read_eukcc(f):
    d = {
        "completeness": 0,
        "contamination": 0,
        "set_species": "NA",
        "number_markers": 0,
        "marker_prevalence": 98,
    }
    with open(f) as fin:
        c = csv.DictReader(fin, delimiter="\t")
        for row in c:
            d["completeness"] = row["completeness"]
            d["contamination"] = row["contamination"]
            # EukCC 2 writes fasta, completeness, contamination, ncbi_lng; ngenomes and n
            # only exist in EukCC 1 output
            d["set_species"] = row.get("ngenomes", "NA")
            d["number_markers"] = row.get("n", "NA")
            if "ncbi_lng" in row:
                d["estimator"] = "EukCC 2"
            break
    return d


def read_eukcc2(f):
    with gzip.GzipFile(f, "r") as fin:
        data = json.loads(fin.read().decode("utf-8"))
    # construct a dict that is concistent
    d = {
        "completeness": data["quality"]["completeness"],
        "contamination": data["quality"]["contamination"],
        "set_species": data["scmg_species"],
        "number_markers": len(data["profiles"]),
        "marker_prevalence": data["scmg_prevalence"],
    }
    return d


def main(**kwargs):
    # load all tables into memory
    estimator = "unkown"
    es = {}
    if kwargs["estimated"].endswith("tsv"):
        # this is eukcc1
        estimator = "EukCC"
        es = read_eukcc(kwargs["estimated"])
        estimator = es.pop("estimator", estimator)
    elif kwargs["estimated"].endswith("gz"):
        # this is eukcc2
        es = read_eukcc2(kwargs["estimated"])
        estimator = "EukCC 2"
    else:
        logging.error("Not yet implemeneted this ending")
        estimator = "BUSCO?"
    si = lcsv(kwargs["simulated"])
    # merged dict
    m = {**es, **si}
    m["estimator"] = estimator

    with open("summed_{}_{}.tsv".format(estimator, m["simulation_prefix"]), "w") as fout:
        cw = csv.DictWriter(fout, list(m.keys()), delimiter="\t")
        cw.writeheader()
        cw.writerow(m)
    logging.info("Done")
    return


if __name__ == "__main__":
    # define arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--estimated", "-e", help="estimated qualit")
    parser.add_argument("--simulated", "-s", help="simulated qualit")
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
