#!/usr/bin/env python3
"""Structural checks on a packaged EukCC2 database. Exits 1 on failure.

Checks, per sub-database: every path in db_info.json exists; the placement blocks add up
to the tree alignment; every profile in <clade>_search.hmm carries a GA line; every leaf
is in taxinfo. Across sub-databases: base is the union of the clades, and no hold-out
genome made it into any taxinfo.
"""
import argparse
import glob
import json
import os
import sys

HMM_KEYS = {"hmm_db", "hmm_placement_search"}


def fasta_lengths(path):
    """{name: sequence length}"""
    out, name = {}, None
    with open(path) as fin:
        for line in fin:
            if line.startswith(">"):
                name = line[1:].split()[0]
                out[name] = 0
            elif name is not None:
                out[name] += len(line.strip())
    return out


def taxinfo(path):
    with open(path) as fin:
        return {line.strip().split(",")[1] for line in fin if line.strip()}


def check_subdb(root, sub, errors):
    info = json.load(open(os.path.join(root, sub, "db_info.json")))
    files = info["files"]
    base = os.path.join(root, sub)
    for key, rel in files.items():
        p = os.path.normpath(os.path.join(base, rel))
        ok = os.path.exists(p) or (key in HMM_KEYS and os.path.exists(p + ".h3f"))
        if not ok:
            errors.append(f"{sub}: {key} -> {rel} does not exist")
    if "etedb" not in files:
        errors.append(f"{sub}: db_info.json has no etedb")

    tree = fasta_lengths(os.path.join(base, files["backbone_alignment"]))
    width = next(iter(tree.values()), 0)
    blocks = sorted(glob.glob(os.path.join(base, files["hmm_placement_place"], "*.aln")))
    block_sum = sum(next(iter(fasta_lengths(b).values()), 0) for b in blocks)
    if block_sum != width:
        errors.append(f"{sub}: placement blocks add up to {block_sum}, tree alignment is {width}")

    with open(os.path.join(base, files["hmm_placement_search"])) as fin:
        lines = fin.read().splitlines()
    names = sum(1 for x in lines if x.startswith("NAME "))
    gas = sum(1 for x in lines if x.startswith("GA "))
    if names != gas or names != len(blocks):
        errors.append(f"{sub}: search.hmm has {names} profiles, {gas} GA lines, {len(blocks)} blocks")

    tax = taxinfo(os.path.join(base, files["taxinfo"]))
    missing = set(tree) - tax
    if missing:
        errors.append(f"{sub}: {len(missing)} leaves not in taxinfo, e.g. {sorted(missing)[0]}")
    print(f"{sub}: {len(tree)} leaves, {len(blocks)} markers, alignment {width} columns", file=sys.stderr)
    return tax


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("dbdir")
    ap.add_argument("--classification", help="species_classification.csv, to check the hold-out")
    args = ap.parse_args()

    errors = []
    subs = sorted(d for d in os.listdir(args.dbdir) if d.startswith("db_"))
    tax = {s: check_subdb(args.dbdir, s, errors) for s in subs}

    clades = [s for s in subs if s != "db_base"]
    if "db_base" in tax and clades:
        union = set().union(*(tax[c] for c in clades))
        if tax["db_base"] != union:
            errors.append(f"db_base has {len(tax['db_base'])} genomes, the clades together {len(union)}")

    if args.classification:
        with open(args.classification) as fin:
            holdout = {r.split("\t")[1] for r in fin if r.split("\t")[2] == "testing"}
        leaked = holdout & set().union(*tax.values())
        if leaked:
            errors.append(f"{len(leaked)} hold-out genomes are in the database, e.g. {sorted(leaked)[0]}")

    for e in errors:
        print(f"FAIL {e}", file=sys.stderr)
    if errors:
        sys.exit(1)
    print("all checks passed", file=sys.stderr)


if __name__ == "__main__":
    main()
