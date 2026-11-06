#!/usr/bin/env python
# coding: utf-8

import logging
from collections import Counter, defaultdict
import argparse
import os
import csv
from random import seed, shuffle
from ete3 import NCBITaxa

ncbi = NCBITaxa()
seed(124367665)


# backup fasta handler, so we can use readonly directories
class fa_class:
    def __init__(self, seq, name, long_name):
        self.seq = seq
        self.name = name
        self.long_name = long_name

    def __str__(self):
        return self.seq

    def __len__(self):
        return len(self.seq)


def Fasta(path):
    """
    Iterator for fasta files
    """
    entry = False

    with open(path) as fin:
        for line in fin:
            if line.startswith(">"):
                if entry is not False:
                    entry.seq = "".join(entry.seq)
                    yield entry
                # define new entry
                long_name = line.strip()[1:]
                name = long_name.split()[0]
                entry = fa_class([], name, long_name)
            else:
                entry.seq.append(line.strip())
        # yield last one
        entry.seq = "".join(entry.seq)
        yield entry


def read_tsv(path, header=True, sep="\t"):
    """
    simple csv reader
    """
    head = None
    content = []

    with open(path) as fin:
        for line in fin:
            line = line.strip()
            if header and head is None:
                head = line.split(sep)
                continue
            else:
                elements = line.split(sep)
            if head is not None and header:
                content.append(dict(zip(head, elements)))
            else:
                content.append(elements)
    return content


def filter_genome(acc, at_least=500000):
    p = "genomes/{}.fna".format(acc)
    return os.path.exists(p)


def getRanks(
    taxid,
    desired_ranks=["domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"],
    fillNA=True,
):
    lng = ncbi.get_lineage(taxid)
    lineage2ranks = ncbi.get_rank(lng)
    ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
    ranks = {"{}".format(rank): ranks2lineage.get(rank, "NA") for rank in desired_ranks}
    # for all ranks with NA value, we insert a joined statement of the lineage below
    # the next non NA rank
    if fillNA:
        for rank in desired_ranks:
            if ranks[rank] == "NA":
                next_none_na = None
                for r in desired_ranks[desired_ranks.index(rank) + 1 :]:
                    if ranks[r] != "NA":
                        next_none_na = ranks[r]
                        break
                ranks[rank] = "-".join([str(x) for x in lng[: lng.index(next_none_na)]])

    return ranks


def sample_evenly(
    d, level="phylum", at_most=20, per_clade=2, ignore=None, only_refseq=False, min_bucket=1, protect=None
):
    if ignore is None:
        ignore = set()
    # protected ranks (for the hold-out: phylum and order): a species is sampled only if
    # its taxon at every protected rank keeps at least min_bucket - 1 species behind, so
    # with min_bucket=2 the sample never takes the last genome of a phylum or an order
    # out of the database, in any round. Taxa are counted over every species, before the
    # ignore list
    protect = list(protect or []) if min_bucket > 1 else []
    size = {rk: Counter(v["lng"][rk] for v in d.values()) for rk in protect}
    taken = {rk: Counter() for rk in protect}

    def eligible(sp):
        return all(size[rk][d[sp]["lng"][rk]] - taken[rk][d[sp]["lng"][rk]] > min_bucket - 1 for rk in protect)

    # add genbanks to ignore if we want refseq only
    if only_refseq:
        for k, v in d.items():
            if v["acc"].startswith("GCA"):
                ignore.add(k)

    # get a dict sorted by the balanced level
    sample_list = defaultdict(set)
    for k, v in d.items():
        # skip any species we want to ignore
        if k in ignore:
            continue
        taxid = v["lng"][level]
        sample_list[taxid].add(k)

    # shuffle all clades real good
    # sorted lists shuffled with the module seed: the same input gives the same sample
    # (a set's order would depend on PYTHONHASHSEED). key=str: buckets mix taxids and
    # the pseudo-phylum strings built by getRanks(fillNA=True).
    sample_list = {k: sorted(v) for k, v in sorted(sample_list.items(), key=lambda kv: str(kv[0]))}
    for l in sample_list.values():
        shuffle(l)
    avail = sum([len(x) for k, x in sample_list.items()])

    logging.debug("There are {} species to sample from".format(avail))
    sampled = set()
    while avail > 0 and len(sampled) < at_most:
        for clade, species in sample_list.items():
            if len(species) > 0:
                # the last eligible species of the shuffled bucket (without protected
                # ranks, simply the last one)
                pick = next((i for i in range(len(species) - 1, -1, -1) if eligible(species[i])), None)
                if pick is None:
                    avail = avail - len(species)
                    species.clear()
                    continue
                sp = species.pop(pick)
                sampled.add(sp)
                for rk in protect:
                    taken[rk][d[sp]["lng"][rk]] += 1
                avail = avail - 1
    logging.debug("Done sampling")
    return set(sampled)


def main(**kwargs):
    # read in the info csv
    table = read_tsv(kwargs["info"], header=False, sep=",")
    logging.debug("We loaded info for {} entries".format(len(table)))

    # for each table element we need to get the lineage information
    # so we work with a dictionary
    d = {}
    for row in table:
        acc = row[2]
        p = "genomes/{}.fna".format(acc)
        p2 = "genomes/{}.faa".format(acc)
        if not os.path.exists(p):
            if not os.path.exists(p2):
                logging.info("Genome {} does not exist".format(row[2]))
                continue

        # save genome info
        d[row[1]] = {"lng": getRanks(row[1]), "acc": row[2], "taxid": row[1], "group": row[0]}
    logging.debug("{} Entries ended up as valid entries".format(len(d)))

    # now that we have the ineage information
    # we can sample across the level of interest
    # groups = [{"name": "test_set", "size": 1, "taxids": set()}, {"name": "training_set", "size": 1, "taxids": set()}]
    # sumrequired = sum([i["size"] for i in groups])

    # sample testing and training set
    if kwargs["testing"] is True:
        testing_set = sample_evenly(
            d,
            level=kwargs["balance_level"],
            at_most=kwargs["n_testing"],
            per_clade=1,
            only_refseq=True,
            min_bucket=kwargs["testing_min_bucket"],
            protect=[r for r in kwargs["testing_protect"].split(",") if r],
        )
        logging.info("Sampled {} species as testing set".format(len(testing_set)))
    else:
        testing_set = set()

    training_set = sample_evenly(
        d, level=kwargs["balance_level"], at_most=kwargs["n_training"], per_clade=2, ignore=testing_set.copy()
    )
    logging.info("Sampled {} species as training set".format(len(training_set)))

    metaeuk_set = sample_evenly(
        d,
        level=kwargs["balance_level"],
        at_most=kwargs["n_metaeuk"],
        per_clade=1,
        ignore=testing_set.copy(),
        only_refseq=True,
    )
    logging.info("Sampled {} species as metaeuk set".format(len(metaeuk_set)))

    # remaining set
    remainer = set(d.keys()) - (testing_set | training_set)
    logging.info("{} species remain".format(len(remainer)))

    p = "species_classification.csv"
    logging.debug("Categorizing {} entries".format(len(d)))
    x = defaultdict(int)
    with open(p, "w") as fout:
        writer = csv.writer(fout, delimiter="\t")
        for k, v in d.items():
            classified = "NA"
            if k in testing_set:
                classified = "testing"
            if k in training_set:
                classified = "training"
            if k in remainer:
                classified = "remaining"
            x[classified] += 1
            r = [v["group"], v["acc"], classified, k, ";".join([str(x) for x in ncbi.get_lineage(k)])]
            writer.writerow(r)
            # metaeuk set
            if k in metaeuk_set:
                classified = "metaeuk"
                r = [v["group"], v["acc"], classified, k, ";".join([str(x) for x in ncbi.get_lineage(k)])]
                writer.writerow(r)

    print(x)
    return


if __name__ == "__main__":
    # define arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("prefix", help="prefix for the output", type=str)
    parser.add_argument("--info", help="taxonomic information", type=str)
    parser.add_argument("--testing", help="run quietly", action="store_true", default=False)
    parser.add_argument("--balance_level", help="level on which balancing is done", default="phylum", type=str)
    parser.add_argument("--minbp", default=1000000, help="Minimal baspairs required", type=int)
    # sample sizes. Sampling goes round-robin over the phylum buckets and checks the size
    # after each full round, so a sample can exceed these numbers
    parser.add_argument("--n_testing", default=30, help="hold-out size", type=int)
    parser.add_argument("--n_training", default=45, help="calibration set size", type=int)
    parser.add_argument("--n_metaeuk", default=25, help="MetaEuk reference size", type=int)
    parser.add_argument(
        "--testing_min_bucket",
        default=1,
        type=int,
        help="the hold-out leaves at least this many - 1 species of every protected taxon",
    )
    parser.add_argument(
        "--testing_protect",
        default="",
        type=str,
        help="comma-separated ranks the hold-out never empties, e.g. phylum,order",
    )
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
