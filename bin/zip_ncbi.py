#!/usr/bin/env python3

import argparse
from collections import defaultdict
import logging
from statistics import median
from ete3 import NCBITaxa

ncbi = NCBITaxa()

# Embryophyta (land plants) and Metazoa
ignorelist = set(["3193", "33208"])


def merge_ncbi(data, path, stats, sizes, names):
    with open(path) as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            else:
                row = line.strip().split("\t")
                stats["rows"] += 1
                # keep representative/reference genomes with a full assembly. NCBI retired
                # the "representative genome" label, current summaries only use "reference"
                if row[4] not in ("representative genome", "reference genome") or row[13] != "Full":
                    continue
                stats["reference_full"] += 1
                # get acc and taxid
                acc = row[0]
                taxid = row[6]
                # genome_size (column 26) exists in summaries since 2022
                if len(row) > 25 and row[25].isdigit():
                    sizes[acc] = int(row[25])
                names[acc] = row[7]
                if data[taxid].startswith("GCF"):
                    # we have a refseq entry for this taxid, so we continue
                    logging.debug("For taxid {} we have acc {}".format(taxid, data[taxid]))
                    continue
                elif data[taxid] == "":
                    # use whatever we get
                    data[taxid] = acc
                elif data[taxid].startswith("GCA") and acc.startswith("GCF"):
                    # replace genbank wirth refseq
                    data[taxid] = acc

    return data


def get_lineage(taxid):
    """[] for a taxid the taxonomy does not know (newer than the dump, or deleted)"""
    try:
        return ncbi.get_lineage(taxid)
    except ValueError:
        return []


def resolve_by_name(name):
    """taxid of the most specific taxon an organism name still resolves to: the full
    name, then the first two words, then the genus"""
    words = name.split()
    for n in (name, " ".join(words[:2]), words[0] if words else ""):
        if n:
            hits = ncbi.get_name_translator([n]).get(n)
            if hits:
                return str(hits[0])
    return None


def lineage_key(taxid):
    """phylum of a taxid; without one, the highest rank below it that is set"""
    lng = get_lineage(taxid)
    ranks = ncbi.get_rank(lng)
    for rank in ("phylum", "class", "order", "family", "genus"):
        for t in lng:
            if ranks.get(t) == rank:
                return t
    return taxid


def genbank_fill(data, stats):
    """keep a GenBank-only species only if no RefSeq species shares its phylum"""
    keys = {taxid: lineage_key(taxid) for taxid in data}
    refseq_keys = {keys[t] for t, acc in data.items() if acc.startswith("GCF")}
    kept = defaultdict(str)
    for taxid, acc in data.items():
        if acc.startswith("GCA") and keys[taxid] in refseq_keys:
            stats["genbank_fill_dropped"] += 1
            continue
        kept[taxid] = acc
    return kept


def relative_size(data, sizes, min_fraction, stats):
    """drop genomes smaller than min_fraction of the median of their genus (or, with
    fewer than 3 sized species in the genus, of their family): partial assemblies and
    organelle-only records filed as reference genomes"""
    parent = {}
    for taxid in data:
        lng = get_lineage(taxid)
        ranks = ncbi.get_rank(lng)
        parent[taxid] = {r: t for t in lng for r in ("genus", "family") if ranks.get(t) == r}
    members = defaultdict(list)
    for taxid, acc in data.items():
        if acc in sizes:
            for rank, t in parent[taxid].items():
                members[(rank, t)].append(sizes[acc])
    kept = defaultdict(str)
    for taxid, acc in data.items():
        ref = None
        for rank in ("genus", "family"):
            group = members.get((rank, parent[taxid].get(rank)), [])
            if len(group) >= 3:
                ref = median(group)
                break
        if acc in sizes and ref and sizes[acc] < min_fraction * ref:
            logging.info("{} ({} bp) is below {} of its relatives ({} bp)".format(acc, sizes[acc], min_fraction, ref))
            stats["relative_size_dropped"] += 1
            continue
        kept[taxid] = acc
    return kept


PILOT_RANKS = ("phylum", "class", "order", "family", "genus")


def taxon_path(taxid):
    """(phylum, class, order, family, genus) of a taxid; a rank that is not set takes the
    taxon above it, and a missing phylum the highest rank below it (lineage_key)"""
    lng = get_lineage(taxid)
    ranks = ncbi.get_rank(lng)
    at = {ranks.get(t): t for t in lng}
    path = []
    for rank in PILOT_RANKS:
        path.append(at.get(rank, path[-1] if path else lineage_key(taxid)))
    return path


def pilot(data, sizes, n, stats):
    """a small, deterministic subset for a pilot run that keeps the diversity of the
    division: n species, each pick the species of the highest taxon not covered yet (a new
    phylum, then a new class, order, family, genus), from the least sampled phylum, class...
    first; ties go to RefSeq and smaller genomes, so the hold-out and the MetaEuk reference
    still have RefSeq genomes to draw from"""
    paths = {taxid: taxon_path(taxid) for taxid in data}
    counts = defaultdict(int)
    left = dict(data)
    kept = defaultdict(str)

    def key(taxid):
        seen = [counts[(i, t)] for i, t in enumerate(paths[taxid])]
        new = next((i for i, c in enumerate(seen) if c == 0), len(seen))
        acc = left[taxid]
        return (new, seen, not acc.startswith("GCF"), sizes.get(acc, 0), acc)

    while left and len(kept) < n:
        taxid = min(left, key=key)
        kept[taxid] = left.pop(taxid)
        for i, t in enumerate(paths[taxid]):
            counts[(i, t)] += 1
    stats["pilot_dropped"] = len(data) - len(kept)
    return kept


def main(files, outputfile, label, statsfile, fill, min_relative_size, n_pilot=None):
    # simple dict storing taxid and accession for each entry
    data = defaultdict(str)
    stats = defaultdict(int)
    sizes = {}
    names = {}

    # iterate over the files given
    for path in files:
        data = merge_ncbi(data, path, stats, sizes, names)
    stats["species"] = len(data)

    # a taxid NCBI deleted (delnodes.dmp) that the summary still uses: fall back to the
    # taxon the organism name resolves to (species, else genus)
    for taxid in list(data):
        if get_lineage(taxid):
            continue
        acc = data.pop(taxid)
        new = resolve_by_name(names.get(acc, ""))
        if new and new not in data:
            logging.warning("taxid {} of {} ({}) is not in the taxonomy: using {}".format(taxid, acc, names.get(acc), new))
            stats["taxid_from_name"] += 1
            data[new] = acc
        else:
            logging.warning("taxid {} of {} ({}) is not in the taxonomy: dropped".format(taxid, acc, names.get(acc)))
            stats["unknown_taxid"] += 1

    if fill:
        data = genbank_fill(data, stats)

    # land plants and animals out before the relative size check, so they do not
    # weigh on the medians
    excluded = set()
    for taxid in data:
        if set(str(x) for x in get_lineage(taxid)) & ignorelist:
            stats["excluded_lineage"] += 1
            excluded.add(taxid)
    data = defaultdict(str, {t: a for t, a in data.items() if t not in excluded})
    if min_relative_size:
        data = relative_size(data, sizes, min_relative_size, stats)
    if n_pilot:
        data = pilot(data, sizes, n_pilot, stats)

    with open(outputfile, "w") as fout:
        for taxid, acc in data.items():
            stats["candidates"] += 1
            fout.write("{},{},{}\n".format(label, taxid, acc))

    if statsfile:
        with open(statsfile, "w") as fout:
            steps = ["rows", "reference_full", "species", "taxid_from_name", "unknown_taxid"]
            if fill:
                steps.append("genbank_fill_dropped")
            steps.append("excluded_lineage")
            if min_relative_size:
                steps.append("relative_size_dropped")
            if n_pilot:
                steps.append("pilot_dropped")
            for step in steps + ["candidates"]:
                fout.write("{}\t{}\t{}\n".format(label, step, stats[step]))


parser = argparse.ArgumentParser(
    description="Merge ncbi summary files, givin refseq higher priority, reducing to representative genome"
)
parser.add_argument("--label", help="label", default="None")
parser.add_argument("--stats", help="write the count after each filter to this file", default=None)
parser.add_argument(
    "--genbank_fill",
    action="store_true",
    default=False,
    help="GenBank genomes only for phyla without any RefSeq species",
)
parser.add_argument(
    "--min_relative_size",
    type=float,
    default=None,
    help="drop genomes below this fraction of the median size of their genus/family",
)
parser.add_argument("--pilot", type=int, default=None, help="keep only this many species (pilot run)")
parser.add_argument("output", type=str, help="path to output file")
parser.add_argument("summaries", nargs="+", help="path to all the sumamry files")
parser.add_argument("--quiet", action="store_true", help="supress information", default=False)
parser.add_argument("--debug", action="store_true", help="Make it more verbose", default=False)

args = parser.parse_args()

# define logging
logLevel = logging.INFO
if args.quiet:
    logLevel = logging.WARNING
elif args.debug:
    logLevel = logging.DEBUG
logging.basicConfig(format="%(asctime)s %(message)s", datefmt="%m/%d/%Y %H:%M:%S: ", level=logLevel)

args = parser.parse_args()
main(args.summaries, args.output, args.label, args.stats, args.genbank_fill, args.min_relative_size, args.pilot)
