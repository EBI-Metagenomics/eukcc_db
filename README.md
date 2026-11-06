# eukcc-db

Builds the reference databases of [EukCC](https://github.com/EBI-Metagenomics/EukCC) 2
from NCBI genomes and the PANTHER library, following chapter 4 of [Saary (2022)](https://libcatalog.embl.de/cgi-bin/koha/opac-detail.pl?biblionumber=36464).

The code available here is [legacy code](https://drive.google.com/file/d/1cTVMqHifC5HRAQCzPFSFB-V68G-YFeBY/view?usp=sharing) shared by Paul Saary, with some adjustments to allow reproducibility. `docs/workflow.md` explains each step of the workflow, following the diagram in `assets`, and which steps (modules) were changed from the legacy code and why.

The output is a directory EukCC reads with `--db`: four sub-databases (`base`, `fungi`,
`plant`, `protozoa`), each with a reference tree, its placement markers and the
calibrated profile HMMs, plus the NCBI taxonomy they were built with and the provenance
of every step.

## How a database is built

- NCBI assembly_summary.txt per division (fungi, plant, protozoa; RefSeq + GenBank)
  - reference genomes with a Full assembly, minus the bad references of `assets/exclude_genomes.tsv`,
    one per species, RefSeq first
  - fungi: GenBank only for phyla without any RefSeq species
  - drop Embryophyta and Metazoa; drop assemblies < 2 Mb or < 10% of their genus median
  - hold-out: ~30 RefSeq species, one per phylum at a time, never the last genome of a phylum or an order,
     kept out of every database (it can be larger: 43 species with NCBI in October 2026)
  - MetaEuk reference: annotated RefSeq proteomes of ~25 species (default), one per phylum at a time
  - MetaEuk protein prediction for every remaining genome, with the standard genetic code (as EukCC)
  - per sub-database (base = all genomes; fungi, plant, protozoa = the NCBI division):
       calibration set (~45 species by phylum) vs PANTHER -> one bit score cut-off per
       profile (GA line) -> every genome searched with `--cut_ga`
       - SCMG table (profiles with exactly one hit per genome)
       - greedy selection of placement markers
       - best hit per genome and marker -> FAMSA -> trimAl -> gap filter -> supermatrix
       - ModelFinder (once, on base) -> IQ-TREE -> midpoint root

## Requirements

- Nextflow >= 25.04, Java 17.
- Docker or Singularity

## Running

```bash
# the NCBI taxonomy of a month (required)
TAXDUMP=https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2026-11-01.zip

# EBI Codon: SLURM + Singularity, resources per process in conf/codon.config
nextflow run . -profile codon_slurm --outdir results --taxdump $TAXDUMP

# local
nextflow run . -profile docker --outdir results --taxdump $TAXDUMP

# Apple silicon: the images are linux/amd64
nextflow run . -profile docker,arm --outdir results --taxdump $TAXDUMP

# Test
nextflow run . -profile test,docker -stub-run --outdir results_stub

# Evaluation on MAGs simulated from the hold-out genomes
nextflow run . -profile docker --mode testing --db results/eukcc2_db_<YYYY_MM> --outdir results
```

### A monthly release

The release `2026_11`: genomes released up to the cut-off date and the NCBI taxonomy of
that month, with the tree model fixed (`--iqmodel` is only used with `--modelfinder false`).

```bash
nextflow run . -profile codon_slurm --outdir results_2026_11 \
    --name 2026_11 \
    --cutoff_date 2026-11-01 \
    --taxdump https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2026-11-01.zip \
    --modelfinder false --iqmodel LG+R10
```

The database goes to `results_2026_11/eukcc2_db_2026_11/`, with a `README.md` that records
the version, the sources and the parameters of the build.

## Main parameters

| parameter | default | |
|---|---|---|
| `--mode` | `build` | `build` a database, or `testing` to evaluate one |
| `--outdir` | required | |
| `--name` | `YYYY_MM` of `--cutoff_date`, otherwise of the run | the database goes to `<outdir>/eukcc2_db_<name>/` |
| `--cachedir` | `cache` | genomes and PANTHER, kept between runs (`storeDir`) |
| `--taxdump` | required | the NCBI dump, e.g. `taxdump_archive/taxdmp_YYYY-MM-01.zip` (monthly archive); loaded into the ete3 `taxa.sqlite` that ships with the database; the genetic codes always come from it |
| `--taxdb` | | a ready `taxa.sqlite` instead |
| `--panther` / `--panther_url` | PANTHER 16.0 public tarball (3.9 GB) | a directory holding `books/` |
| `--summaries` | download | a directory with `<catalogue>_<group>.txt`, to pin the NCBI snapshot |
| `--pilot_species` | | pilot run: keep only N species per division, spread over the taxonomy (a new phylum first, then a new class, order, family, genus; RefSeq and small genomes first on ties); no files to prepare |
| `--cutoff_date` | | `YYYY-MM-DD`: only assemblies released up to that date |
| `--exclude_genomes` | `assets/exclude_genomes.tsv` | bad reference genomes that never enter the database (see below); `false` = none |
| `--minbp` | 2000000 | minimum assembly size (the smallest nuclear genomes, Microsporidia, are ~2.2 Mb) |
| `--min_relative_size` | 0.1 | drop genomes below this fraction of the median size of their genus (or family): partial assemblies; `null` = off |
| `--genbank_fill_groups` | `['fungi']` | divisions where GenBank genomes enter only for phyla without RefSeq species |
| `--holdout_min_species` | 2 | the hold-out never takes a protected taxon below this many species, in any round: 2 keeps at least one genome of each |
| `--holdout_protect_ranks` | `['phylum', 'order']` | ranks the hold-out never empties |
| `--holdout_size` / `--metaeuk_size` / `--calibration_size` | 30 / 25 / 45 | sample sizes; sampling takes one species per phylum at a time and checks the size only after going through all phyla, so a sample can be larger |
| `--max_leaf_gap_fraction` | 0.9 | leaves with more gaps in the concatenation leave the tree; `null` = report only |
| `--modelfinder` / `--iqmodel` | true / | select the model, or fix it: `--iqmodel` is required with `--modelfinder false` (thesis: `LG+R7`) |
| `--iqtree` | true | false = FastTree |
| `--db` | | `--mode testing`: the database to evaluate (`eukcc2_db_<name>/`) |
| `--contamination_source` | `self` | `--mode testing`: contaminant of each simulated MAG: the genome itself (`self`), another hold-out species (`other`) or `both` |

### Excluded genomes

Some NCBI reference genomes can have low quality and are kept out.
`--exclude_genomes` is a TSV with a header and two columns, `accession` and `reason`;
lines starting with `#` are comments.

## Outputs

```
<outdir>/
├── eukcc2_db_<name>/            the database (--db for EukCC)
│   ├── README.md                version (YYYY_MM), sources, parameters and contents
│   ├── db_base/
│   ├── db_fungi/
│   ├── db_plant/
│   ├── db_protozoa/
│   ├── etetoolkit/
│   │   └── taxa.sqlite          the taxonomy of the build
│   └── provenance/              build_summary.tsv, exclusions.tsv, versions.yml,
│                                params.json, genome roles, cut-offs, markers, greedy
│                                log, leaf gaps; not used by EukCC
├── provenance/                  the same, plus NCBI summaries, hmmsearch --cut_ga
│                                tables, ModelFinder report, taxonomy.tsv (which NCBI
│                                dump); not used by EukCC
├── itol/
│   └── <subdb>/                 tree + iTOL annotation files; not used by EukCC
├── metaeuk_db/                  the MetaEuk reference database (a copy ships in db_base/)
└── pipeline_info/               Nextflow reports
```

Here we include some outputs that are not present on the original workflow:

### Traceability

- `provenance/build_summary.tsv`: for each sub-database, the rows NCBI listed, the count
  after each filter (category and assembly level, one per species, lineage, download,
  size), the hold-out, the genomes of the sub-database, the calibration set, the PANTHER
  families, the calibrated profiles, the SCMG candidates (profiles single-copy in at least
  one genome) and SCMGs per genome, the placement markers, the leaves and the alignment
  width.
- `provenance/coverage.tsv`: per sub-database, group (fungi, plant, protozoa) and phylum, the
  candidates, hold-out, genomes in the sub-database, leaves and genomes with a
  non-standard genetic code.
- `provenance/exclusions.tsv`: every candidate that did not become a leaf, with the stage
  and the reason (listed in `--exclude_genomes`, download, size, hold-out, no proteome,
  search failed, too many gaps, no placement marker).
- `provenance/exclude_genomes.tsv`: the `--exclude_genomes` list of the build; its md5 is
  in the release `README.md`.
- `provenance/markers/find_tree_MG_<subdb>.log`: every iteration of the greedy marker selection.

### iTOL

`itol/<subdb>/` holds the rooted tree and annotation files to drop onto
[iTOL](https://itol.embl.de): `<subdb>_labels.txt` (species names), and colour strips of
the sub-database and the phylum of each leaf on base (`base_subdb.txt`,
`base_phylum.txt`), and phylum, class, order and family on the others
(`<subdb>_phylum.txt` ...). Colours come from the hue family of the sub-database (fungi red to orange,
plant greens, protozoa blue to purple; grey for "no rank"), so every category reads as
the sub-database its genomes are in. A tree whose colours break into many interleaved blocks does
not follow the NCBI taxonomy.

## Reproducibility

The genome samples (hold-out, MetaEuk reference, calibration sets) use fixed seeds over
sorted lists, the greedy selection breaks ties with a fixed seed, and IQ-TREE gets
`-seed`. Two runs on the same inputs give the same genome roles, cut-offs, markers, SCMG
table and alignments.

NCBI changes every day. To rebuild a database, pass the summaries of its build
(`<outdir>/provenance/ncbi/summaries/`) with `--summaries`, its `etetoolkit/taxa.sqlite`
with `--taxdb`, and the same `--taxdump`, recorded in `<outdir>/provenance/taxonomy/taxonomy.tsv`:
the genetic codes always come from the dump.

## Reference

Saary, P. J. (2022). *Recovery and quality estimation of metagenomic assembled genomes of
eukaryotes*. PhD thesis, University of Cambridge (EMBL-EBI), chapter 4.
