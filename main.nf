#!/usr/bin/env nextflow

// Names used in the comments:
//   group      NCBI division: fungi, plant or protozoa
//   catalogue  refseq or genbank
//   clade      sub-database: base, fungi, plant or protozoa (base = every genome)
//   accession  NCBI assembly accession (GCF_... / GCA_...), also the name of each genome file
//   marker     a PANTHER family chosen as placement marker (PTHR...)

include { inputs }           from './subworkflows/inputs'
include { genome_selection } from './subworkflows/genome_selection'
include { proteomes }        from './subworkflows/proteomes'
include { build_db }         from './subworkflows/build_db'
include { traceability }     from './subworkflows/traceability'
include { testing }          from './subworkflows/testing'


workflow build {
    main:
    // 1. Inputs: taxonomy, NCBI assembly summaries, PANTHER
    inputs()

    // 2. Genome selection: candidates, downloads, hold-out, MetaEuk set
    genome_selection(inputs.out.summaries, inputs.out.taxdb)

    // 3. Proteomes: MetaEuk reference and protein prediction
    proteomes(genome_selection.out.group_genomes,
              genome_selection.out.metaeuk_proteomes)

    // 4. Build db: one run per clade, the three groups plus base with every genome
    group_proteomes = proteomes.out.proteomes.groupTuple()           // (group, [proteomes])
    base_proteomes = proteomes.out.proteomes
        .map { _group, proteome -> tuple('base', proteome) }
        .groupTuple()                                                // ('base', [proteomes])
    build_db(group_proteomes.mix(base_proteomes),
             genome_selection.out.ncbi_table,
             proteomes.out.metaeuk_files,
             inputs.out.taxdb,
             inputs.out.panther_db,
             inputs.out.panther)

    // 5. Traceability: build report, provenance, release directory
    traceability(inputs.out.summaries,
                 genome_selection.out.zip_stats,
                 genome_selection.out.ncbi_table,
                 genome_selection.out.downloaded,
                 genome_selection.out.sized,
                 genome_selection.out.classification,
                 proteomes.out.proteomes,
                 build_db.out.searched,
                 build_db.out.report_files,
                 build_db.out.db,
                 inputs.out.panther_db,
                 inputs.out.taxdb,
                 inputs.out.taxonomy_info,
                 genome_selection.out.genetic_codes)
}


workflow {
    if (!params.outdir) {
        error "--outdir is required"
    }
    if (params.mode == 'build') {
        if (!params.taxdump) {
            error "--taxdump is required, e.g. https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_YYYY-MM-01.zip"
        }
        // Nextflow 26 passes --flag false from the command line as the string 'false'
        if (!params.modelfinder.toString().toBoolean() && !params.iqmodel) {
            error "--modelfinder false needs --iqmodel, e.g. LG+R7"
        }
        build()
    } else if (params.mode == 'testing') {
        testing()
    } else {
        error "--mode must be build or testing, not '${params.mode}'"
    }
}
