#!/usr/bin/env python

"""iTOL annotation files for a reference tree, to check it against the NCBI taxonomy.

<prefix>_labels.txt      LABELS: leaf accession -> "Species name (accession)"
<prefix>_<level>.txt     DATASET_COLORSTRIP, one per --levels entry: an NCBI rank
                         (phylum, class, order, family, ...) or "subdb", the
                         sub-database of the genome (fungi, plant, protozoa)

Drag the tree and the files onto https://itol.embl.de (help: itol.embl.de/help.cgi#annot).
"""

import argparse
import colorsys
import csv
import os
from collections import Counter
from ete3 import NCBITaxa, Tree

# every colour comes from the hue family of a sub-database, so a category reads as the
# sub-database its genomes are in: fungi red to orange, plant greens, protozoa blue to purple
HUES = {"fungi": (0, 40), "plant": (80, 160), "protozoa": (190, 290)}
SUBDB_COLORS = {"fungi": "#e6550d", "plant": "#31a354", "protozoa": "#3182bd"}
MISSING = "#d9d9d9"


def colorstrip(path, label, assignments, colors):
    """assignments: {leaf: category}; colors: {category: hex}"""
    legend = [c for c in colors if c in set(assignments.values())]
    with open(path, "w") as fout:
        fout.write("DATASET_COLORSTRIP\nSEPARATOR TAB\n")
        fout.write("DATASET_LABEL\t{}\nCOLOR\t#000000\nSTRIP_WIDTH\t25\n".format(label))
        fout.write("LEGEND_TITLE\t{}\n".format(label))
        fout.write("LEGEND_SHAPES\t" + "\t".join("1" for _ in legend) + "\n")
        fout.write("LEGEND_COLORS\t" + "\t".join(colors[c] for c in legend) + "\n")
        fout.write("LEGEND_LABELS\t" + "\t".join(legend) + "\n")
        fout.write("DATA\n")
        for leaf in sorted(assignments):
            cat = assignments[leaf]
            fout.write("{}\t{}\t{}\n".format(leaf, colors[cat], cat))


def shades(n, hues):
    """n distinguishable colours inside a hue range: hue spread over the range, lightness
    cycling so neighbours in the range still differ"""
    lo, hi = hues
    out = []
    for i in range(n):
        h = (lo + (hi - lo) * i / (n - 1)) if n > 1 else (lo + hi) / 2
        r, g, b = colorsys.hls_to_rgb(h / 360.0, (0.42, 0.60, 0.30)[i % 3], 0.75)
        out.append("#{:02x}{:02x}{:02x}".format(int(round(r * 255)), int(round(g * 255)), int(round(b * 255))))
    return out


def read_legend(path):
    """{category: colour} from the legend of an existing colour strip"""
    labels, colors = [], []
    with open(path) as fin:
        for line in fin:
            f = line.rstrip("\n").split("\t")
            if f[0] == "LEGEND_LABELS":
                labels = f[1:]
            elif f[0] == "LEGEND_COLORS":
                colors = f[1:]
    return dict(zip(labels, colors))


def family(color):
    """the sub-database whose hue range holds a colour (None for grey)"""
    r, g, b = (int(color[i:i + 2], 16) / 255.0 for i in (1, 3, 5))
    h, _l, sat = colorsys.rgb_to_hls(r, g, b)
    if sat < 0.1:
        return None
    for sub, (lo, hi) in HUES.items():
        if lo - 1 <= h * 360 <= hi + 1:
            return sub
    return None


def distance(a, b):
    return sum((int(a[i:i + 2], 16) - int(b[i:i + 2], 16)) ** 2 for i in (1, 3, 5)) ** 0.5


def palette(assignments, missing, subdb_of, reference=None):
    """Each category gets a shade of the sub-database most of its leaves are in; within
    a sub-database, the most frequent categories come first.

    With a reference ({category: colour} from an earlier annotation), a category keeps its
    colour if it is still in the same sub-database. Categories that are new, or moved to
    another sub-database, get the shade of their family farthest from every colour
    already in the strip, so no other colour changes."""
    counts = Counter(v for v in assignments.values() if v != missing)
    home = {}
    for cat in counts:
        subs = Counter(subdb_of.get(leaf) for leaf, c in assignments.items() if c == cat)
        home[cat] = subs.most_common(1)[0][0]
    colors = {}
    for sub in sorted(set(home.values()), key=str):
        cats = sorted((c for c in counts if home[c] == sub), key=lambda c: (-counts[c], c))
        if sub not in HUES:
            colors.update({c: MISSING for c in cats})
        elif reference is None:
            colors.update(zip(cats, shades(len(cats), HUES[sub])))
        else:
            kept = {c: reference[c] for c in cats if c in reference and family(reference[c]) == sub}
            colors.update(kept)
            candidates = shades(60, HUES[sub])
            for cat in cats:
                if cat in kept:
                    continue
                used = [v for v in list(colors.values()) + list(reference.values()) if family(v) == sub]
                colors[cat] = max(candidates, key=lambda c: min((distance(c, u) for u in used), default=0))
    if missing in assignments.values():
        colors[missing] = MISSING
    return colors


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--taxinfo", required=True, help="taxid,accession of the leaves")
    ap.add_argument("--groups", required=True, help="group,taxid,accession (ncbi_table.txt)")
    ap.add_argument("--taxdb", required=True, help="ete3 taxa.sqlite")
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--levels", nargs="+", required=True)
    ap.add_argument(
        "--reference",
        help="directory with earlier annotation files (<prefix>_<level>.txt): categories keep their colours",
    )
    args = ap.parse_args()

    ncbi = NCBITaxa(dbfile=args.taxdb)
    leaves = Tree(args.tree, format=1).get_leaf_names()
    taxid = {}
    with open(args.taxinfo) as fin:
        for t, acc in csv.reader(fin):
            taxid[acc] = int(t)
    group = {}
    with open(args.groups) as fin:
        for g, _t, acc in csv.reader(fin):
            group[acc] = g

    lineage = {}
    for t in set(taxid.values()):
        try:
            lng = ncbi.get_lineage(t)
        except ValueError:
            lng = []
        ranks = ncbi.get_rank(lng)
        names = ncbi.get_taxid_translator(lng)
        lineage[t] = {ranks[x]: names[x] for x in lng}
        lineage[t]["_name"] = names.get(lng[-1], str(t)) if lng else str(t)

    with open("{}_labels.txt".format(args.prefix), "w") as fout:
        fout.write("LABELS\nSEPARATOR TAB\nDATA\n")
        for leaf in sorted(leaves):
            name = lineage[taxid[leaf]]["_name"] if leaf in taxid else "NA"
            fout.write("{}\t{} ({})\n".format(leaf, name, leaf))

    for level in args.levels:
        if level == "subdb":
            assignments = {leaf: group.get(leaf, "NA") for leaf in leaves}
            colors = dict(SUBDB_COLORS)
            colors.update({c: MISSING for c in set(assignments.values()) - set(colors)})
        else:
            missing = "no {}".format(level)
            assignments = {
                leaf: lineage[taxid[leaf]].get(level, missing) if leaf in taxid else missing for leaf in leaves
            }
            ref_file = os.path.join(args.reference, "{}_{}.txt".format(args.prefix, level)) if args.reference else None
            reference = read_legend(ref_file) if ref_file and os.path.exists(ref_file) else None
            colors = palette(assignments, missing, group, reference)
        colorstrip("{}_{}.txt".format(args.prefix, level), level, assignments, colors)
        print("{} {}: {} categories over {} leaves".format(args.prefix, level, len(set(assignments.values())), len(leaves)))


if __name__ == "__main__":
    main()
