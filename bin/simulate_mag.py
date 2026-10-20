#!/usr/bin/env python
import logging
import argparse
import random
import textwrap
import os


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


def genome_size(path):
    return sum([len(rec.seq) for rec in Fasta(path)])


def decision(probability):
    return random.random() < probability


def mergeOverlaps(intervals):
    """https://codereview.stackexchange.com/a/69249"""
    sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
    merged = []

    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1] = [lower[0], upper_bound]  # replace by merged interval
            else:
                merged.append(higher)

    return merged


def sample_from_genome(path, bp):
    """
    simple function to sample n bp from a genome
    """
    fragment_sizes = list(range(5000, 50000, 523))
    gs = genome_size(path)
    if bp == 0:
        fraction = 0
    else:
        fraction = float(bp) / float(gs)
    if bp > gs:
        logging.debug("Cant sample more than what we need, will add as much as we can")
        fraction = 1
        bp = gs
    logging.debug("Will retain {}bp which are {}% of the genome".format(bp, round(100 * fraction, 2)))

    seqs = []
    sum_len = 0
    sample_space = [rec for rec in Fasta(path)]
    while len(sample_space) > 0 and sum_len < bp and bp > 0:
        for rec in sample_space:
            n = 0
            max_n = len(rec.seq)
            frags = []
            # scan along the chromosome
            while n < max_n:
                fragment_sample_bp = random.choice(fragment_sizes)
                if decision(fraction):
                    if sum_len + fragment_sample_bp > bp:
                        # retain only the n first bp, as we have to much
                        fragment_sample_bp = int(bp - sum_len)
                    frags.append([n, n + fragment_sample_bp])

                # increase end point
                n = n + fragment_sample_bp
            # merge overlaps as to not overly fragment everything:
            frags = mergeOverlaps(frags)

            # now for each frag add it as a sequences
            for f in frags:
                seqs.append(rec.seq[f[0] : f[1]])

            # count how much we have
            sum_len = sum([len(x) for x in seqs])

            # remove rec from sample_space
            sample_space.remove(rec)

            if sum_len >= bp:
                break

    logging.debug("sampled {} bp".format(sum_len))
    return seqs


def break_line(s, lw=120):
    return textwrap.fill(s, width=lw)


def main(**kwargs):

    # parse fractions
    fractions = [float(x) / 100 for x in kwargs["fraction"].split("/")]
    genomes = [kwargs["genome"], kwargs["contaminate"]]

    # size of genome A
    logging.info("Measuring the input genome")
    size_genome = genome_size(kwargs["genome"])
    logging.debug("Calculating how many bp I have to keep and add")
    sizes = [int(round(f * size_genome, 0)) for f in fractions]

    mag_seqs = []
    logging.debug("Aggregating sequences")
    for genome, bp in zip(genomes, sizes):
        seqs = sample_from_genome(genome, bp)
        n = sum([len(x) for x in seqs])
        mag_seqs.append({"bp": n, "seqs": seqs, "percent": round(100 * n / size_genome, 2)})

    logging.debug("Writing MAG to file")
    with open("{}_mag.fasta".format(kwargs["prefix"]), "w") as fout:
        # write refernce
        i = 1
        logging.debug("Writing genome MAG to file")
        for seq in mag_seqs[0]["seqs"]:
            fout.write(">reference_{}\n".format(i))
            for line in textwrap.wrap(seq):
                fout.write("{}\n".format(line))
            i += 1
        # write contamiantion
        i = 1
        logging.debug("Writing cont MAG to file")
        for seq in mag_seqs[1]["seqs"]:
            fout.write(">contaminant_{}\n".format(i))
            for line in textwrap.wrap(seq):
                fout.write("{}\n".format(line))
            i += 1

    # write stats
    logging.debug("Writing stats to file")
    with open("{}_stats.tsv".format(kwargs["prefix"]), "w") as fout:
        compl = mag_seqs[0]["percent"]
        cont = mag_seqs[1]["percent"]
        fout.write(
            "simulation_genome\tsimulation_contaminant\tsimulated_completeness\tsimulated_contamination\tsimulation_prefix\n"
        )
        fout.write(
            "{}\t{}\t{}\t{}\t{}\n".format(
                os.path.basename(kwargs["genome"]),
                os.path.basename(kwargs["contaminate"]),
                compl,
                cont,
                kwargs["prefix"],
            )
        )


if __name__ == "__main__":
    # define arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome", "-g", help="Base genome", type=str, required=True)
    parser.add_argument("--contaminate", "-c", help="Base genome", type=str, required=True)
    parser.add_argument("--prefix", "-p", help="output prefix", type=str, default="simulated")
    parser.add_argument("--fraction", "-f", help="completeness/contamination", type=str, required=True)
    parser.add_argument("--quiet", help="run quietly", action="store_true", default=False)
    parser.add_argument("--debug", help="run quietly", action="store_true", default=False)
    # seeds the fragment sampling: the same call gives the same MAG
    parser.add_argument("--seed", help="random seed", type=int, default=None)
    args = parser.parse_args()
    random.seed(args.seed)

    # define logging
    logLevel = logging.INFO
    if args.quiet:
        logLevel = logging.WARNING
    elif args.debug:
        logLevel = logging.DEBUG
    logging.basicConfig(format="%(asctime)s %(message)s", datefmt="%d/%m/%Y %H:%M:%S: ", level=logLevel)

    main(**vars(args))
