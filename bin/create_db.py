#!/usr/bin/env python3
import shutil
import os
from datetime import datetime
import argparse
import logging
import json
import re
import tarfile


def main(**opt):
    logging.info("Creating DB")
    for key, value in opt.items():
        logging.info("{}: {}".format(key, value))

    ww = {}
    clade = opt["clade"]
    logging.info("Creating DB for clade {}".format(clade))
    db_dir = "db_{}".format(opt["clade"])
    os.makedirs(db_dir, exist_ok=True)
    os.makedirs(os.path.join(db_dir, "backbone", "placement", "place"), exist_ok=True)
    os.makedirs(os.path.join(db_dir, "scmgs"), exist_ok=True)

    if clade == "base":
        logging.info("Copy metaeuk")
        logging.info("we have a metaeuk db and this is the base db")
        ww["metaeukdb"] = "metaeukdb/microeuks_metaeuk.db"
        os.makedirs(os.path.join(db_dir, "metaeukdb"), exist_ok=True)
        for f in opt["metaeuk"]:
            if not f.endswith(".db"):
                # microeuks_metaeuk.db.index -> .index (split on "." would give .db.index)
                suffix = f.split(".db", 1)[1]
            else:
                suffix = ""
            p = os.path.join(db_dir, ww["metaeukdb"] + suffix)
            shutil.copyfile(f, p)

    logging.info("Copy taxinfo")
    # taxinfo
    ww["taxinfo"] = "backbone/{}_taxinfo.csv".format(clade)
    shutil.copyfile(opt["taxinfo"], os.path.join(db_dir, ww["taxinfo"]))

    # SCMGs
    logging.info("Copy SCMGs")
    ww["scmgs"] = "scmgs/{}_SCMGs.csv.gz".format(clade)
    shutil.copyfile(opt["scmgs"], os.path.join(db_dir, ww["scmgs"]))

    # tree
    logging.info("Copy tree")
    ww["backbone_tree"] = "backbone/{}_rooted.tre".format(clade)
    shutil.copyfile(opt["backbone_tree"], os.path.join(db_dir, ww["backbone_tree"]))

    # backbone alignement
    ww["backbone_alignment"] = "backbone/{}_tree.aln".format(clade)
    shutil.copyfile(opt["backbone_alignment"], os.path.join(db_dir, ww["backbone_alignment"]))

    # hmmdb with indices
    logging.info("Copy Hmm db")
    ww["hmm_db"] = "scmgs/{}_PANTHERv{}.hmm".format(clade, opt["panther"])
    for f in opt["hmm_db"]:
        if not f.endswith(".hmm"):
            suffix = f.rsplit(".", 1)[1]
            suffix = ".{}".format(suffix)
        else:
            suffix = ""
        p = os.path.join(db_dir, ww["hmm_db"] + suffix)
        shutil.copyfile(f, p)

    # hmmplacement with indices
    logging.info("Copy hmm placement search files")
    ww["hmm_placement_search"] = "backbone/placement/{}_search.hmm".format(clade)
    for f in opt["placement_hmm"]:
        if not f.endswith(".hmm"):
            suffix = f.rsplit(".", 1)[1]
            suffix = ".{}".format(suffix)
        else:
            suffix = ""
        p = os.path.join(db_dir, ww["hmm_placement_search"] + suffix)
        shutil.copyfile(f, p)

    # placement profiles aln and hmms
    logging.info("Copy many aln and hmms")
    ww["hmm_placement_place"] = "backbone/placement/place/"
    profile = re.compile(r"PTHR[0-9]+")
    for aln in opt["places_aln"]:
        m = profile.search(aln)
        name = m.group(0)
        logging.debug("Copy {}: {}".format(name, aln))
        p = os.path.join(db_dir, ww["hmm_placement_place"], "{}.aln".format(name))
        shutil.copyfile(aln, p)
    for aln in opt["places_hmm"]:
        m = profile.search(aln)
        name = m.group(0)
        logging.debug("Copy {}: {}".format(name, aln))
        p = os.path.join(db_dir, ww["hmm_placement_place"], "{}.hmm".format(name))
        shutil.copyfile(aln, p)

    # the taxonomy the database was built with. Without this key EukCC resolves taxids
    # with the user's own ~/.etetoolkit. package_db puts the file at the release root.
    ww["etedb"] = "../etetoolkit/taxa.sqlite"

    logging.info("creating DB json")
    db_json = {}
    db_json["files"] = ww.copy()
    db_json["name"] = "EukCC2 database"
    db_json["db_clade"] = clade
    db_json["db_version"] = 2
    db_json["db_creation_date"] = str(datetime.date(datetime.now()))
    db_json["iqtree_model"] = opt["iqtree"]
    db_json["panther_version"] = opt["panther"]
    print(db_json)
    with open(os.path.join(db_dir, "db_info.json"), "w") as fp:
        json.dump(db_json, fp)
    return db_dir


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--clade", help="taxonomic information", type=str)
    parser.add_argument("--backbone_tree", help="backbone tree", type=str)
    parser.add_argument("--backbone_alignment", help="backbone tree", type=str)
    parser.add_argument("--scmgs", help="taxonomic information", type=str)
    parser.add_argument("--hmm_db", help="taxonomic information", type=str, nargs="+")
    parser.add_argument("--placement_hmm", help="taxonomic information", type=str, nargs="+")
    parser.add_argument("--places_aln", help="taxonomic information", type=str, nargs="+")
    parser.add_argument("--places_hmm", help="taxonomic information", type=str, nargs="+")
    parser.add_argument("--taxinfo", help="taxonomic information", type=str)
    parser.add_argument("--iqtree", help="IQ tree model, needed for epa ng", type=str)
    parser.add_argument("--panther", help="panther version used", type=str)
    parser.add_argument("--metaeuk", help="taxonomic information", type=str, nargs="+")

    parser.add_argument("--tar", help="run quietly", action="store_true", default=False)
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

    d = main(**vars(args))

    if args.tar:
        logging.info("Creating tar archive")
        tar = tarfile.open("db_{}.tar.gz".format(args.clade), "w:gz")
        tar.add(d, arcname="eukcc2_db_{}".format(args.clade))
        tar.close()
