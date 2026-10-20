#!/usr/bin/env python
# coding: utf-8

import os
import logging
from collections import defaultdict
import argparse
from collections import Counter
from random import sample, seed
from operator import itemgetter

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
    name = filename.replace("hmmerGA_", "").rsplit(".", 1)[0]
    logging.debug("Reading file of acc {}".format(name))

    # reduce to what we need:
    perProfile = defaultdict(list)
    for row in data:
        profile = row[2].split(".", 1)[0]
        perProfile[profile].append(name)

    return perProfile


def score_profile(lst, ignore=None, lenient=False):
    """
    score a profile based on how many SCMG are contained
    if many score is high, if few, score is low
    we can also pass a number of accesions to ignore for the score

    The second score traces how many singletons in total would be covered
    This is to rank profiles even if the number of novel singleton coverage
    is the same, the hidden coverage (Of genomes that are already
    well covered) might be different
    """
    if ignore is None:
        ignore = set()
    # in lenient mode we allow also doubletons to be counted as a win
    if lenient:
        max_cov = 2
    else:
        max_cov = 1
    score = 0
    hidden_score = 0
    for key, value in Counter(lst).items():

        increase = 0
        if value <= max_cov:
            increase += 1
            if value == 1:
                hidden_score += 1
        elif value > 2:
            # doubletons get no negative penalty
            increase += -value

        # if we have enough markers for a given species,
        # any marker gets a score of 0
        if key in ignore and increase > 0:
            increase = 0
        score = score + increase
    return (score, hidden_score)


def greedy_marker_selection(data, profiles=None, n_target=10, maxiter=500, soft_target=5):

    iter_i = 0
    ignore = set()
    softignore = set()
    keys = list(data.keys())
    universe = set()
    for k, v in data.items():
        universe = universe | set(v)

    # at the end we might want to be lenient and also allow doubletns
    lenient = False

    # allow the exclusion of already selected profiles
    if profiles is None:
        choosen_profiles = set()
    else:
        choosen_profiles = profiles.copy()
    keys = [key for key in keys if key not in choosen_profiles]
    # ignore any accession that has more than n_targets
    score = defaultdict(int)
    while iter_i < maxiter:
        logging.debug("iteration {}".format(iter_i))
        iter_i += 1

        scores = []
        for profile in keys:
            genomes = data[profile]
            scores.append(score_profile(genomes, ignore=ignore, lenient=lenient))
        # get max positions
        scores_sort = sorted(scores, key=itemgetter(0, 1), reverse=True)
        max_idx = [idx for idx, scores in enumerate(scores) if scores == scores_sort[0]]
        # allow ay that covers either max or 20
        # max_idx = [idx for idx, scores in enumerate(scores) if scores == scores_sort[0] or scores[0] >= 15]
        print(max_idx)
        idx = sample(max_idx, 1)[0]

        # switch to lenient if this matches
        if lenient is False and scores_sort[0][0] < 10 and scores_sort[0][1] < 50:
            logging.info("switching to lenient scan, top score: {}".format(scores_sort[0]))
            lenient = True

        if scores_sort[0][0] < 0.04 * (len(universe) - len(ignore)):
            logging.info("No profile covers more than 1% of the remaining universe, so we stop here")
            break

        if scores_sort[0][1] < 5:
            logging.warning("No profile covers enough (3) genomes anymore")
            i = 0
            for s in scores_sort:
                print(s)
                i += 1
                if i > 5:
                    break
            break

        logging.debug("Added profile {} with scores {} ".format(keys[idx], scores[idx]))
        choosen_profiles.add(keys[idx])
        new_profile = keys[idx]
        del keys[idx]

        for acc in set(data[new_profile]):
            score[acc] += 1

        for key, value in score.items():
            if value >= n_target:
                ignore.add(key)
            elif value >= soft_target:
                softignore.add(key)

        logging.debug("Found enough profiles for {} out of {}".format(len(ignore), len(universe)))
        if len(ignore | softignore) == len(universe):
            logging.info("Found enough profiles to cover everything")
            break

    stats = {}
    for acc, cov in score.items():
        # identify multitons
        n = 0
        for profile in choosen_profiles:
            if data[profile].count(acc) > 1:
                n += 1
        logging.debug("Covered {} with {} profiles, {} of which are multitons".format(acc, cov, n))
        stats[acc] = {"covered": cov, "multitons": n}

    return (choosen_profiles, stats)


def main(**kwargs):
    # load all tables into memory
    logging.info("Reading in all the data")
    data = defaultdict(list)
    for tbl in kwargs["tbls"]:
        perprofile = read_hmmer(tbl)
        for profile, d in perprofile.items():
            data[profile].extend(d)

    n_target = 8
    n_soft = 4
    logging.info("Greedy marker selection")
    profiles, stats = greedy_marker_selection(data, n_target=n_target, soft_target=n_soft, maxiter=40)
    # for the worst we launch again
    worst_data = defaultdict(list)
    offenders = set()
    for acc, stat in stats.items():
        if stat["covered"] < n_target:
            offenders.add(acc)
    logging.info("Found {} genomes with low coverage".format(len(offenders)))
    for profile, accs in data.items():
        worst_data[profile] = [x for x in accs if x in offenders]

    profiles_1, stats_new = greedy_marker_selection(worst_data, n_target=n_soft, soft_target=n_soft, maxiter=10)
    profiles, stats = greedy_marker_selection(
        data, profiles=profiles_1, n_target=n_target, soft_target=n_soft, maxiter=40
    )

    # profiles = greedy_marker_selection(data)
    clade = kwargs["clade"]
    with open("profiles.csv", "w") as fout:
        for p in profiles:
            fout.write("{},{}\n".format(clade, p))
    return


if __name__ == "__main__":
    # define arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--clade", help="which clade", type=str, default="base")
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
