// Step 2 of docs/workflow.md (2a-2h): candidate genomes, downloads, hold-out and MetaEuk set.

include { genetic_codes }                                                  from '../modules/prepare'
include { filter_release_date; zip_summary; download_genome; genome2faa }  from '../modules/ncbi'
include { size_filter; choose_representative }                             from '../modules/general'

workflow genome_selection {
    take:
    summaries       // (catalogue, group, <catalogue>_<group>.txt)
    taxdb           // taxa.sqlite

    main:
    // the hand-made list of genomes that never enter the database; [] = none. On the command
    // line --exclude_genomes false (or null) arrives as a string
    exclude_genomes = params.exclude_genomes && !(params.exclude_genomes.toString() in ['false', 'null', 'true'])
        ? file(params.exclude_genomes, checkIfExists: true) : []

    // 2a. only with --cutoff_date. in: (catalogue, group, summary), cut-off date
    dated = params.cutoff_date ? filter_release_date(summaries, params.cutoff_date) : summaries

    // 2b. one genome per species. in: (group, [the RefSeq and GenBank summaries of the group], taxdb),
    // GenBank fill groups, min_relative_size and pilot_species ('' = off), exclude_genomes ([] = none)
    zip_summary(dated
        .groupTuple(by: 1)                                           // ([catalogues], group, [summaries])
        .map { _catalogues, group, group_summaries -> tuple(group, group_summaries) }
        .combine(taxdb),
        params.genbank_fill_groups,
        params.min_relative_size ?: '',
        params.pilot_species ?: '',
        exclude_genomes)
    // the candidates of the three groups in one file: group,taxid,accession
    ncbi_table = zip_summary.out.table.map { _group, candidates_csv -> candidates_csv }
        .collectFile(name: 'ncbi_table.txt', sort: { candidates_csv -> candidates_csv.name }).first()

    // 2c. in: taxdump, ncbi_table; out: genetic_codes.tsv (accession, taxid, nuclear code)
    genetic_codes(file(params.taxdump), ncbi_table)

    // 2d. in: (group, taxid, accession); out: (group, <accession>.fna)
    candidates = ncbi_table.splitCsv().map { row -> tuple(row[0], row[1], row[2]) }.unique()
    download_genome(candidates)

    // 2e. in: (group, genome), minbp; out: the same, only for genomes >= minbp
    size_filter(download_genome.out, params.minbp)

    // 2f. in: ncbi_table, [every genome that passed the size filter], taxdb, minbp,
    // hold-out size/min species/protected ranks, calibration and MetaEuk sizes
    choose_representative(ncbi_table, size_filter.out.map { _group, genome -> genome }.collect(), taxdb,
        params.minbp, params.holdout_size, params.holdout_min_species, params.holdout_protect_ranks,
        params.calibration_size, params.metaeuk_size)

    // 2g. the roles of every genome: [accession: [roles]] (a MetaEuk genome has two rows).
    // The genomes keep the path size_filter gave them, no copies: a change in one group
    // (e.g. --exclude_genomes) does not move the files of the others, and -resume keeps
    // their proteomes and sub-databases
    roles = choose_representative.out.info.map { classification_csv ->
        classification_csv.splitCsv(sep: '\t')
            .groupBy { row -> row[1] }
            .collectEntries { accession, rows -> [(accession): rows.collect { row -> row[2] }] }
    }
    // every genome outside the hold-out: (group, genome)
    group_genomes = size_filter.out.combine(roles)
        .filter { _group, genome, role -> role[genome.baseName] && !('testing' in role[genome.baseName]) }
        .map { group, genome, _role -> tuple(group, genome) }
        .ifEmpty { error "No genome left for the databases after the hold-out: every candidate went to it (--holdout_size ${params.holdout_size})" }

    // 2h. in: one genome of the MetaEuk set; out: its annotated RefSeq proteome
    genome2faa(size_filter.out.combine(roles)
        .filter { _group, genome, role -> 'metaeuk' in (role[genome.baseName] ?: []) }
        .map { _group, genome, _role -> genome })

    emit:
    ncbi_table = ncbi_table                                          // ncbi_table.txt
    zip_stats = zip_summary.out.stats                                // <group>_stats.tsv
    genetic_codes = genetic_codes.out                                // genetic_codes.tsv
    downloaded = download_genome.out                                 // (group, genome)
    sized = size_filter.out                                          // (group, genome)
    classification = choose_representative.out.info                  // species_classification.csv
    group_genomes = group_genomes                                    // (group, genome) outside the hold-out
    metaeuk_proteomes = genome2faa.out                               // <accession>.faa
    exclude_genomes = exclude_genomes                                // exclude_genomes.tsv, or [] without a list
}
