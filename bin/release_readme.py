#!/usr/bin/env python3
"""README.md of a database release: its YYYY_MM version, where every input came from,
the parameters of the build, its contents and how to rebuild it.

Reads only what the release directory already holds (provenance/params.json, run.txt,
taxonomy.tsv, exclude_genomes.tsv, versions.yml, build_summary.tsv, db_*/db_info.json), plus the NCBI
assembly_summary.txt files of the run, to record their md5.

    release_readme.py --db eukcc2_db_<name>/ --output eukcc2_db_<name>/README.md
        --summaries <catalogue>_<group>.txt ...
"""
import argparse
import csv
import glob
import hashlib
import json
import os
import sys

SUBDBS = ("base", "fungi", "plant", "protozoa")


def read_kv(path, sep):
    """key<sep>value lines -> dict; {} when the file is not there"""
    if not os.path.exists(path):
        return {}
    with open(path) as fin:
        return dict(line.rstrip("\n").split(sep, 1) for line in fin if sep in line)


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fin:
        for chunk in iter(lambda: fin.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def rows(path):
    with open(path) as fin:
        return sum(1 for line in fin if not line.startswith("#"))


def version(params):
    """YYYY_MM of the cut-off date, otherwise of the day the run started"""
    day = params.get("cutoff_date") or params["run_stamp"]
    return str(day)[:7].replace("-", "_")


def flag(value):
    """a parameter as Nextflow writes it: true/false, not True/False"""
    return json.dumps(value) if isinstance(value, bool) else value


def table(header, body):
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    out += ["| " + " | ".join(str(x) for x in r) + " |" for r in body]
    return "\n".join(out)


def readme(db, summaries):
    prov = os.path.join(db, "provenance")
    params = json.load(open(os.path.join(prov, "params.json")))
    run = read_kv(os.path.join(prov, "run.txt"), ": ")
    tax = read_kv(os.path.join(prov, "taxonomy.tsv"), "\t")
    tools = read_kv(os.path.join(prov, "versions.yml"), ": ")
    counts = {}
    with open(os.path.join(prov, "build_summary.tsv")) as fin:
        for r in csv.DictReader(fin, delimiter="\t"):
            counts[(r["subdb"], r["step"])] = r["count"]
    info = {}
    for f in glob.glob(os.path.join(db, "db_*", "db_info.json")):
        i = json.load(open(f))
        info[i["db_clade"]] = i
    base = info.get("base") or next(iter(info.values()))
    ver = version(params)
    dirname = os.path.basename(os.path.normpath(db))
    run_day = params["run_stamp"][:10]
    cutoff = params.get("cutoff_date")

    def c(sub, step):
        return f"{int(counts[(sub, step)]):,}" if (sub, step) in counts else "-"

    contents = table(
        ["sub-database", "leaves", "placement markers", "alignment columns", "calibrated profiles", "SCMG profiles"],
        [[s, c(s, "leaves"), c(s, "placement_markers"), c(s, "alignment_columns"),
          c(s, "calibrated_profiles"), c(s, "scmg_profiles")] for s in SUBDBS if s in info],
    )

    if params.get("summaries"):
        ncbi_src = f"`{params['summaries']}/<catalogue>_<group>.txt` (pinned with `--summaries`)"
    else:
        ncbi_src = (f"`{params['ncbi_url']}/genomes/<catalogue>/<group>/assembly_summary.txt`, "
                    f"downloaded by the run that started on {run_day}")
    ncbi = table(["file", "rows", "md5"],
                 [[f"`{os.path.basename(p)}`", f"{rows(p):,}", f"`{md5(p)}`"] for p in sorted(summaries)])

    if params.get("taxdb"):
        tax_src = f"a ready `taxa.sqlite` (`--taxdb {params['taxdb']}`)"
    else:
        tax_src = (f"`{tax.get('source', params.get('taxdump'))}` (`--taxdump`); nodes.dmp dated "
                   f"{tax.get('ncbi_date', '?')}, md5 `{tax.get('md5', '?')}`, loaded with ete3 {tax.get('ete3', '?')}")
    excl = os.path.join(prov, "exclude_genomes.tsv")
    if os.path.exists(excl):
        n = sum(1 for line in open(excl) if line.strip() and not line.startswith(("#", "accession\t")))
        excl_src = (f"`{params['exclude_genomes']}` (`--exclude_genomes`): {n} accessions, md5 `{md5(excl)}`, "
                    f"copied to `provenance/exclude_genomes.tsv`")
    else:
        excl_src = "none"
    panther_src = (f"local copy `{params['panther']}` (`--panther`) of `{params['panther_url']}`"
                   if params.get("panther") else f"`{params['panther_url']}`")

    holdout = c("base", "holdout")
    rebuild = [
        "nextflow run EBI-Metagenomics/eukcc_db -profile codon_slurm \\",
        f"    --name {ver} --outdir <outdir> \\",
        f"    --summaries <outdir of this build>/provenance/ncbi/summaries \\",
        f"    --taxdb {dirname}/etetoolkit/taxa.sqlite \\",
        (f"    --exclude_genomes {dirname}/provenance/exclude_genomes.tsv \\" if os.path.exists(excl)
         else "    --exclude_genomes false \\"),
        f"    --panther_url {params['panther_url']} --panther_version {params['panther_version']} \\",
        (f"    --cutoff_date {cutoff} " if cutoff else "    ")
        + f"--modelfinder {flag(params['modelfinder'])} --iqmodel {params['iqmodel']} --seed {params['seed']}",
    ]

    sel = [
        ("`--cutoff_date`", cutoff or "none", "assemblies released up to this date"),
        ("`--divisions`", ", ".join(f"{g}/{cat}" for g, cat in params["divisions"]), "NCBI divisions read"),
        ("`--genbank_fill_groups`", ", ".join(params["genbank_fill_groups"]) or "none",
         "GenBank genomes only for phyla without RefSeq species"),
        ("`--exclude_genomes`", f"{c('base', 'user_excluded')} dropped" if os.path.exists(excl) else "none",
         "genomes listed by hand that never enter the database"),
        ("`--minbp`", f"{params['minbp']:,}", "minimum assembly size (bp)"),
        ("`--min_relative_size`", params["min_relative_size"], "minimum size as a fraction of the genus median"),
        ("`--holdout_size`", params["holdout_size"], f"RefSeq species kept out of every sub-database ({holdout} taken)"),
        ("`--holdout_protect_ranks`", ", ".join(params["holdout_protect_ranks"]), "ranks the hold-out never empties"),
        ("`--metaeuk_size`", params["metaeuk_size"],
         f"RefSeq proteomes in the MetaEuk reference ({c('base', 'metaeuk_reference')} taken)"),
        ("`--calibration_size`", params["calibration_size"], "species per sub-database that calibrate the GA cut-offs"),
        ("`--evalue`", params["evalue"], "hmmsearch E-value on the calibration genomes"),
        ("`--max_leaf_gap_fraction`", params["max_leaf_gap_fraction"], "leaves with more gaps leave the tree"),
        ("`--iqtree` / `--modelfinder`", f"{flag(params['iqtree'])} / {flag(params['modelfinder'])}", "tree method, model selection"),
        ("`--iqmodel`", params["iqmodel"], "substitution model when ModelFinder is off"),
        ("`--seed`", params["seed"], "IQ-TREE and sampling seed"),
    ]

    return f"""# EukCC2 database {ver}

| | |
|---|---|
| version | **{ver}** (YYYY_MM of {"the cut-off date" if cutoff else "the build"}) |
| directory | `{dirname}/` |
| created | {base.get('db_creation_date', '?')} |
| NCBI cut-off | {cutoff or "none: every assembly in the NCBI summaries below"} |
| pipeline | {run.get('pipeline', '?')}, Nextflow {run.get('nextflow', '?')} |
| tree model | {base.get('iqtree_model', '?')} |
| PANTHER | {base.get('panther_version', '?')} |

Use it with `eukcc --db eukcc2_db_{ver}`.

## Contents

{contents}

`base` holds every genome; `fungi`, `plant` and `protozoa` follow the NCBI division of
each genome. Only `db_base/` carries the MetaEuk reference (`metaeukdb/`).

## Sources

| input | where it came from |
|---|---|
| NCBI genomes | {ncbi_src} |
| NCBI taxonomy | {tax_src} |
| PANTHER library | {panther_src}, version {params['panther_version']} |
| excluded genomes | {excl_src} |

NCBI assembly summaries of this build:

{ncbi}

## Parameters

{table(["parameter", "value", ""], sel)}

Every parameter of the run is in `provenance/params.json`.

## Tools

{table(["tool", "version"], sorted(tools.items()))}

## Provenance

| file | |
|---|---|
| `provenance/build_summary.tsv` | count after each step, per sub-database |
| `provenance/exclusions.tsv` | every candidate that did not become a leaf, with the stage and reason |
| `provenance/coverage.tsv` | taxa covered by each sub-database |
| `provenance/species_classification.csv` | role of every genome (hold-out, MetaEuk, calibration, remaining) |
| `provenance/taxonomy.tsv` | the taxonomy dump: source, md5, date |
| `provenance/exclude_genomes.tsv` | the genomes excluded by hand, with the reason |
| `provenance/params.json`, `run.txt`, `versions.yml` | parameters, run and tool versions |

## Rebuild

NCBI changes every day. To rebuild this release, pin the NCBI snapshot and the taxonomy of
this build (the assembly summaries are in `<outdir>/provenance/ncbi/summaries/` of the run):

```bash
{chr(10).join(rebuild)}
```
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", required=True, help="the release directory, eukcc2_db_<name>/")
    ap.add_argument("--summaries", nargs="+", default=[], help="<catalogue>_<group>.txt as downloaded")
    ap.add_argument("--output", required=True, help="the README to write, e.g. eukcc2_db_<name>/README.md")
    args = ap.parse_args()

    with open(args.output, "w") as fout:
        fout.write(readme(args.db, args.summaries))
    print(f"{args.output}: {len(args.summaries)} NCBI summaries", file=sys.stderr)


if __name__ == "__main__":
    main()
