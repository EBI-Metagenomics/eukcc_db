#!/usr/bin/env python
# coding: utf-8
import os
import logging
from collections import defaultdict
import argparse
from collections import Counter
import operator
from statistics import median
from ete3 import NCBITaxa

# from glob import glob
# import shutil
import csv

ncbi = NCBITaxa()


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
    name = filename.replace("training_", "").replace(".tbl", "")
    logging.debug("Reading file of acc {}".format(name))

    # reduce to what we need:
    perProfile = defaultdict(list)
    for row in data:
        profile = row[2].split(".", 1)[0]
        perProfile[profile].append([name, float(row[5])])

    return perProfile


def median_score(scores):
    """
    function to define the cuttoff based on the median of a plateau
    """
    max_hit = max(x[1] for x in scores)
    cutoffs = [x[0] for x in scores if x[1] == max_hit]
    return [round(median(cutoffs), 2), max_hit]


def max_score(scores):
    """
    function to define the cuttoff based on the median of a plateau
    """
    max_hit = max(x[1] for x in scores)
    cutoffs = [x[0] for x in scores if x[1] == max_hit]
    return [round(max(cutoffs), 2), max_hit]


def test_median():
    scores = [[10, 5], [9, 5], [8, 5], [11, 1], [7, 3]]
    h = median_score(scores)
    assert h == [9, 5]


def score_scmgs(hits, taxinfo, win_for=0.7, max_gnms=100):
    # clades count 10
    hits = Counter(hits)
    score = 0
    if win_for > 1:
        logging.error("winf for xcant be bigger than 1")
        exit(1)

    for clade, genomes in taxinfo.items():
        singletons = sum([1 for genome in genomes if hits[genome] == 1])
        fr = singletons / len(genomes)
        if fr >= win_for:
            score += max_gnms
        else:
            score += fr * max_gnms * 0.8

    return round(score, 2)


def maximize(data, taxinfo):
    # collect all bitscores in this profile
    # we will evaluate the distrubution at each bitscore hit
    steps = [x[1] for x in data]
    steps.append(0)
    # making unique
    steps = list(set(steps))
    steps.sort()

    # actual cutoffs are between two steps always
    cutoffs = [20, 50, 100, 200]
    for i, step in enumerate(steps):
        if i + 1 == len(steps):
            # no need to get a step behond maximum
            break
        cutoffs.append(round((step + steps[i + 1]) / 2, 2))
    # enfore minimal cutoff of 20
    cutoffs = [x for x in cutoffs if x >= 20]

    # for each possible cutoff see how many scmgs are there
    scores = []
    max_gnms = max([len(clade) for key, clade in taxinfo.items()])
    logging.debug("Max gain set as {}".format(max_gnms))
    for cutoff in cutoffs:
        # get genome names for each hit
        hits = [x[0] for x in data if x[1] >= cutoff]
        oldscore = len(set([k for k, v in Counter(hits).items() if v == 1]))
        scmgs = score_scmgs(hits, taxinfo, max_gnms=max_gnms)
        scores.append([cutoff, scmgs, oldscore])

    # which cutoff produces the maximum number of SCMGs
    # sort by SCMGscore and bitscore
    sorted_scores = sorted(scores, key=operator.itemgetter(1, 0), reverse=True)

    return sorted_scores[0]


def main(**kwargs):

    # load all tables into memory
    logging.info("Reading in all the data")
    data = defaultdict(list)
    nspecies = len(kwargs["tbls"])
    genomes = set()
    for tbl in kwargs["tbls"]:
        perprofile = read_hmmer(tbl)
        for profile, d in perprofile.items():
            data[profile].extend(d)
            genomes.add(d[0][0])

    # ncbi info
    # table clade,taxid, acc
    taxinfo = defaultdict(set)
    # partition the data
    with open(kwargs["taxinfo"]) as fin:
        for line in fin:
            row = line.strip().split(",")
            # clade = row[0]
            taxid = row[1]
            acc = row[2]
            if acc not in genomes:
                continue
            lng = ncbi.get_lineage(taxid)
            lineage2ranks = ncbi.get_rank(lng)
            ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
            ranks = {"{}".format(rank): ranks2lineage.get(rank, "NA") for rank in ["clade", "phylum"]}
            groups = [lng[3], ranks["clade"]]
            for x in reversed(groups):
                if x != "NA":
                    taxinfo[x].add(acc)
                    break
    n_m = sum([len(s) for k, s in taxinfo.items()])
    logging.debug(
        "Loaded taxid, we have {} clades with {} members and {} genomes".format(len(taxinfo), n_m, len(genomes))
    )
    # for each profile
    # find bitscore that maximizes the number
    # of single copy occurenc for this profile
    # safe how distant the next change is
    logging.info("Computing cutoffs")
    with open("bitscore_cutoffs.csv", "w") as fout:
        header = ["profile", "cutoff", "score", "max_singletons", "species"]
        writer = csv.writer(fout, delimiter="\t")
        writer.writerow(header)
        # optimize per profile
        for profile, hits in data.items():
            logging.debug("Maximizing profile {}".format(profile))
            x = maximize(hits, taxinfo)
            row = [profile]
            row.extend(x)
            row.append(nspecies)
            writer.writerow(row)

    logging.info("Done writing all results")
    return


if __name__ == "__main__":
    # define arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--taxinfo", help="give taxinfo list")
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
