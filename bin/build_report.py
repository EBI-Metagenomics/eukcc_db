#!/usr/bin/env python3
"""Traceability of a build: where every candidate genome went, and the counts behind
each sub-database.

exclusions.tsv     one row per genome (and sub-database) that did not become a leaf
build_summary.tsv  per sub-database: NCBI rows -> filters -> genomes kept -> profiles
                   calibrated -> SCMGs -> placement markers -> tree
coverage.tsv       per sub-database, group (fungi, plant, protozoa) and phylum: candidates, hold-out,
                   genomes in the sub-database, leaves, non-standard genetic codes
                   (with --taxdb)

Per sub-database files are read from --subdb-dir by name:
species_classification_<s>.csv, bitscore_cutoffs_<s>.csv, profiles_<s>.csv,
leaf_gaps_<s>.tsv, scmgs_<s>.gz and <s>.fasta (the tree alignment).
"""
import argparse
import csv
import glob
import gzip
import os
import sqlite3
import statistics
import sys
from collections import Counter, defaultdict

GROUPS = ("fungi", "plant", "protozoa")


def phyla(path, taxids):
    """{taxid: phylum} from an ete3 taxa.sqlite, read with plain SQL"""
    out = {t: "NA" for t in taxids}
    try:
        con = sqlite3.connect(path)
        track = {}
        for t in taxids:
            row = con.execute("select track from species where taxid = ?", (int(t),)).fetchone()
            if row:
                track[t] = [int(x) for x in row[0].split(",")]
        nodes = {x for lng in track.values() for x in lng}
        info = {}
        nodes = list(nodes)
        for i in range(0, len(nodes), 500):
            chunk = nodes[i:i + 500]
            q = "select taxid, spname, rank from species where taxid in ({})".format(",".join("?" * len(chunk)))
            for taxid, name, rank in con.execute(q, chunk):
                info[taxid] = (name, rank)
    except sqlite3.Error:
        return out
    for t, lng in track.items():
        out[t] = next((info[x][0] for x in lng if info.get(x, ("", ""))[1] == "phylum"), "no phylum")
    return out


def read_exclude(path):
    """{accession: reason} of an exclude_genomes.tsv (same reader as zip_ncbi.py)"""
    out = {}
    with open(path) as fin:
        for line in fin:
            if line.startswith("#") or not line.strip():
                continue
            acc, _, reason = line.rstrip("\n").partition("\t")
            if acc.strip() != "accession":
                out[acc.strip()] = reason.strip()
    return out


def lines(path):
    with open(path) as fin:
        return {x.strip() for x in fin if x.strip()}


def fasta(path):
    """{name: length}"""
    out, name = {}, None
    with open(path) as fin:
        for line in fin:
            if line.startswith(">"):
                name = line[1:].split()[0]
                out[name] = 0
            elif name is not None:
                out[name] += len(line.strip())
    return out


def per_subdb(directory, prefix, suffix):
    """{subdb: path} for <directory>/<prefix><subdb><suffix>"""
    out = {}
    for p in glob.glob(os.path.join(directory, prefix + "*" + suffix)):
        name = os.path.basename(p)[len(prefix):-len(suffix)]
        if name:
            out[name] = p
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--summaries", nargs="+", default=[], help="<catalogue>_<group>.txt as downloaded")
    ap.add_argument("--zip-stats", nargs="+", default=[], help="<group>_stats.tsv from zip_ncbi.py")
    ap.add_argument("--candidates", required=True, help="ncbi_table.txt: group,taxid,accession")
    ap.add_argument("--downloaded", required=True, help="accessions out of download_genome")
    ap.add_argument("--sized", required=True, help="accessions out of size_filter")
    ap.add_argument("--classification", required=True, help="species_classification.csv (global roles)")
    ap.add_argument("--predicted", required=True, help="accessions out of metaeuk_pred")
    ap.add_argument("--searched", required=True, help="'subdb<TAB>accession' out of hmmsearch_GA")
    ap.add_argument("--subdb-dir", required=True)
    ap.add_argument("--panther-families", type=int, default=None)
    ap.add_argument("--minbp", default="minbp")
    ap.add_argument("--taxdb", help="ete3 taxa.sqlite, for coverage.tsv")
    ap.add_argument("--genetic-codes", help="accession, taxid, nuclear genetic code")
    ap.add_argument("--exclude-genomes", help="exclude_genomes.tsv: accessions dropped by zip_ncbi.py --exclude")
    args = ap.parse_args()

    downloaded, sized, predicted = lines(args.downloaded), lines(args.sized), lines(args.predicted)
    roles = defaultdict(set)
    with open(args.classification) as fin:
        for row in csv.reader(fin, delimiter="\t"):
            roles[row[1]].add(row[2])
    searched = defaultdict(set)
    for row in lines(args.searched):
        sub, acc = row.split("\t")
        searched[sub].add(acc)

    d = args.subdb_dir
    alignments = {s: fasta(p) for s, p in per_subdb(d, "", ".fasta").items()}
    gappy, gap_tables = {}, per_subdb(d, "leaf_gaps_", ".tsv")
    for sub, p in gap_tables.items():
        with open(p) as fin:
            for r in csv.DictReader(fin, delimiter="\t"):
                if r["excluded"] == "1":
                    gappy[(sub, r["accession"])] = r["gap_fraction"]

    # ---- one row per genome that did not become a leaf --------------------------------
    exclusions = []
    group_of, taxid_of = {}, {}
    with open(args.candidates) as fin:
        for group, taxid, acc in csv.reader(fin):
            if acc in group_of:
                continue
            group_of[acc] = group
            taxid_of[acc] = taxid
            if acc not in downloaded:
                exclusions.append((acc, group, "-", "download_genome", "download failed"))
            elif acc not in sized:
                exclusions.append((acc, group, "-", "size_filter", "assembly below {} bp".format(args.minbp)))
            elif "testing" in roles[acc]:
                exclusions.append((acc, group, "-", "choose_representative", "hold-out (testing)"))
            elif not roles[acc]:
                exclusions.append((acc, group, "-", "choose_representative", "not classified"))
            elif acc not in predicted:
                exclusions.append((acc, group, "-", "metaeuk_pred", "no proteome predicted"))
            else:
                for sub in (group, "base"):
                    if sub not in searched:
                        continue
                    if acc not in searched[sub]:
                        exclusions.append((acc, group, sub, "hmmsearch_GA", "search failed"))
                    elif (sub, acc) in gappy:
                        exclusions.append((acc, group, sub, "filter_gappy_leaves",
                                           "gap fraction {}".format(gappy[(sub, acc)])))
                    elif sub in alignments and acc not in alignments[sub]:
                        exclusions.append((acc, group, sub, "horizontal", "no placement marker"))

    # the hand-made list: every summary row it matches (an entry without a version matches
    # every version), and a warning for an entry no summary has (typo, suppressed assembly)
    exclude = read_exclude(args.exclude_genomes) if args.exclude_genomes else {}
    found, seen = set(), set()
    for p in args.summaries if exclude else []:
        group = os.path.basename(p).rsplit(".", 1)[0].split("_", 1)[1]
        with open(p) as fin:
            for line in fin:
                acc = line.split("\t", 1)[0]
                key = acc if acc in exclude else acc.split(".")[0]
                if line.startswith("#") or key not in exclude or (acc, group) in seen:
                    continue
                seen.add((acc, group))
                found.add(key)
                exclusions.append((acc, group, "-", "zip_summary",
                                   "listed in exclude_genomes: {}".format(exclude[key])))
    for key in sorted(set(exclude) - found):
        print("WARNING: {} is listed in exclude_genomes but not in the NCBI assembly summaries".format(key), file=sys.stderr)
        exclusions.append((key, "-", "-", "zip_summary", "listed in exclude_genomes, not found in the NCBI assembly summaries"))

    with open("exclusions.tsv", "w") as fout:
        fout.write("accession\tgroup\tsubdb\tstage\treason\n")
        for row in sorted(exclusions):
            fout.write("\t".join(row) + "\n")

    # ---- counts per sub-database ----------------------------------------------------
    rows = []

    def add(sub, step, n, what):
        rows.append((sub, step, n, what))

    # genome funnel: per group, base is the sum of the three
    funnel = defaultdict(Counter)
    for p in args.summaries:
        cat, group = os.path.basename(p).rsplit(".", 1)[0].split("_", 1)
        with open(p) as fin:
            funnel[group]["ncbi_" + cat] += sum(1 for x in fin if not x.startswith("#"))
    for p in args.zip_stats:
        with open(p) as fin:
            for group, step, n in csv.reader(fin, delimiter="\t"):
                funnel[group][step] += int(n)
    for acc, group in group_of.items():
        f = funnel[group]
        f["downloaded"] += acc in downloaded
        f["size_ok"] += acc in sized
        f["holdout"] += "testing" in roles[acc]
        f["metaeuk_reference"] += "metaeuk" in roles[acc]
        f["proteome"] += acc in predicted and "testing" not in roles[acc]
    for g in GROUPS:
        funnel["base"].update(funnel[g])

    steps = [
        ("ncbi_refseq", "rows in the RefSeq assembly_summary.txt"),
        ("ncbi_genbank", "rows in the GenBank assembly_summary.txt"),
        ("rows", "rows after the release-date cut-off (= all rows without --cutoff_date)"),
        ("reference_full", "reference/representative genome with a Full assembly"),
        ("user_excluded", "dropped: listed in exclude_genomes"),
        ("species", "one genome per species (RefSeq preferred)"),
        ("taxid_from_name", "taxid deleted by NCBI, replaced by the taxon the organism name resolves to"),
        ("unknown_taxid", "dropped: taxid not in the taxonomy and organism name not found"),
        ("genbank_fill_dropped", "dropped: GenBank genome of a phylum that has RefSeq species"),
        ("excluded_lineage", "dropped: Embryophyta or Metazoa"),
        ("relative_size_dropped", "dropped: below min_relative_size of the median of its genus/family"),
        ("pilot_dropped", "dropped: not in the pilot subset (--pilot_species)"),
        ("candidates", "candidate genomes"),
        ("downloaded", "downloaded"),
        ("size_ok", "assembly >= {} bp".format(args.minbp)),
        ("holdout", "hold-out (testing), kept out of every database"),
        ("metaeuk_reference", "RefSeq proteomes in the MetaEuk reference (also in the tree)"),
        ("proteome", "genomes with a predicted proteome: the genomes of the sub-database"),
    ]
    cls = per_subdb(d, "species_classification_", ".csv")
    cutoffs = per_subdb(d, "bitscore_cutoffs_", ".csv")
    markers = per_subdb(d, "profiles_", ".csv")
    scmgs = per_subdb(d, "scmgs_", ".gz")

    for sub in ("base",) + GROUPS:
        if sub not in searched and sub not in funnel:
            continue
        for step, what in steps:
            if step in funnel[sub]:
                add(sub, step, funnel[sub][step], what)

        if sub in cls:
            with open(cls[sub]) as fin:
                n = sum(1 for r in csv.reader(fin, delimiter="\t") if r[2] == "training")
            add(sub, "calibration_genomes", n, "genomes that calibrate the bit score cut-offs (~calibration_size, round-robin by phylum)")
        add(sub, "searched", len(searched.get(sub, ())), "genomes searched with the calibrated cut-offs (--cut_ga)")
        if args.panther_families is not None:
            add(sub, "panther_families", args.panther_families, "PANTHER family HMMs")
        if sub in cutoffs:
            with open(cutoffs[sub]) as fin:
                n = sum(1 for _ in fin) - 1
            add(sub, "calibrated_profiles", n, "profiles with a calibrated cut-off (GA line)")
        if sub in scmgs:
            per_genome = Counter()
            profiles = set()
            with gzip.open(scmgs[sub], "rt") as fin:
                for line in fin:
                    acc, profile = line.strip().split(",")
                    per_genome[acc] += 1
                    profiles.add(profile)
            add(sub, "scmg_profiles", len(profiles),
                "profiles single-copy in at least one genome: the SCMG candidates EukCC draws from")
            add(sub, "scmg_rows", sum(per_genome.values()), "genome,profile pairs in the SCMG table")
            if per_genome:
                add(sub, "scmgs_per_genome_median", statistics.median(per_genome.values()),
                    "median number of single-copy profiles per genome")
        if sub in markers:
            with open(markers[sub]) as fin:
                n = sum(1 for x in fin if x.strip())
            add(sub, "placement_markers", n, "markers chosen by the greedy selection (tree and placement)")
        if sub in gap_tables:
            add(sub, "dropped_gappy", sum(1 for s, _ in gappy if s == sub),
                "leaves dropped for too many gaps in the concatenation")
        if sub in alignments:
            aln = alignments[sub]
            add(sub, "leaves", len(aln), "leaves of the reference tree")
            add(sub, "alignment_columns", next(iter(aln.values()), 0), "columns of the concatenated alignment")

    # ---- coverage per group and phylum ----------------------------------------------
    if args.taxdb:
        code = {}
        if args.genetic_codes:
            with open(args.genetic_codes) as fin:
                for acc, _taxid, gc in csv.reader(fin, delimiter="\t"):
                    code[acc] = gc
        phylum_of = phyla(args.taxdb, set(taxid_of.values()))
        cov = defaultdict(Counter)
        for acc, group in group_of.items():
            key = (group, phylum_of[taxid_of[acc]])
            in_db = acc in predicted and "testing" not in roles[acc]
            for sub in (group, "base"):
                c = cov[(sub,) + key]
                c["candidates"] += 1
                c["refseq"] += acc.startswith("GCF")
                c["holdout"] += "testing" in roles[acc]
                c["in_subdb"] += in_db
                c["leaves"] += acc in alignments.get(sub, {})
                c["nonstandard_code"] += code.get(acc, "1") != "1"
        cols = ["candidates", "refseq", "holdout", "in_subdb", "leaves", "nonstandard_code"]
        with open("coverage.tsv", "w") as fout:
            fout.write("subdb\tgroup\tphylum\t" + "\t".join(cols) + "\n")
            for key in sorted(cov, key=lambda k: (k[0] != "base", k[0], -cov[k]["candidates"], k[1], k[2])):
                fout.write("\t".join(key) + "\t" + "\t".join(str(cov[key][c]) for c in cols) + "\n")

    with open("build_summary.tsv", "w") as fout:
        fout.write("subdb\tstep\tcount\tdescription\n")
        for r in rows:
            fout.write("\t".join(str(x) for x in r) + "\n")
    for r in rows:
        print("{:<9} {:<24} {:>8}  {}".format(*r))
    print("{} candidates, {} exclusion rows".format(len(group_of), len(exclusions)), file=sys.stderr)


if __name__ == "__main__":
    main()
