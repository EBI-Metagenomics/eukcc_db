#!/usr/bin/env python3
"""Compare EukCC estimates against the simulated truth, per MAG and per sub-database.

Reads the summed_*.tsv rows written by quality_summary.py (one MAG per file) and reports
the thesis error, |Δcompleteness| + 5·|Δcontamination| (thesis §4.3.1, eq. 4.3).
"""
import argparse
import csv
import statistics as st
import sys
from collections import defaultdict


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--classification", required=True, help="species_classification.csv, to label sub-databases")
    ap.add_argument("tables", nargs="+")
    args = ap.parse_args()

    group = {}
    with open(args.classification) as fin:
        for r in csv.reader(fin, delimiter="\t"):
            group[r[1]] = r[0]

    rows = []
    for path in args.tables:
        with open(path) as fin:
            for r in csv.DictReader(fin, delimiter="\t"):
                acc = r["simulation_genome"].rsplit(".", 1)[0]
                c, k = num(r.get("completeness")), num(r.get("contamination"))
                sc, sk = num(r["simulated_completeness"]), num(r["simulated_contamination"])
                ok = c is not None and k is not None
                rows.append({
                    "mag": r["simulation_prefix"], "genome": acc,
                    "contaminant": r["simulation_contaminant"].rsplit(".", 1)[0],
                    "subdb": group.get(acc, "NA"),
                    "sim_completeness": sc, "sim_contamination": sk,
                    "completeness": c, "contamination": k,
                    "error": abs(sc - c) + 5 * abs(sk - k) if ok else None,
                })

    with open("quality_per_mag.tsv", "w") as fout:
        w = csv.DictWriter(fout, list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: r["mag"]))

    by = defaultdict(list)
    for r in rows:
        by[r["subdb"]].append(r)
        by["all"].append(r)
    with open("quality_summary.tsv", "w") as fout:
        fout.write("subdb\tmags\tfailed\tmedian_abs_dcompl\tmedian_abs_dcont\tmedian_error\n")
        for subdb, rs in sorted(by.items()):
            done = [r for r in rs if r["error"] is not None]

            def med(f):
                return "{:.2f}".format(st.median(f(r) for r in done)) if done else "NA"

            fout.write("\t".join([
                subdb, str(len(rs)), str(len(rs) - len(done)),
                med(lambda r: abs(r["sim_completeness"] - r["completeness"])),
                med(lambda r: abs(r["sim_contamination"] - r["contamination"])),
                med(lambda r: r["error"]),
            ]) + "\n")
    print("{} estimates".format(len(rows)), file=sys.stderr)


if __name__ == "__main__":
    main()
