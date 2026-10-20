// Step 2 of docs/workflow.md (2a-2h): candidate genomes, downloads, hold-out and MetaEuk set.

include { genetic_codes }                                                  from '../modules/prepare'
include { filter_release_date; zip_summary; download_genome; genome2faa }  from '../modules/ncbi'
include { size_filter; choose_representative; split_clades }               from '../modules/general'

workflow genome_selection {
    take:
    summaries       // (catalogue, group, <catalogue>_<group>.txt)
    taxdb           // taxa.sqlite

    main:
    // 2a. only with --cutoff_date. in: (catalogue, group, summary), cut-off date
    dated = params.cutoff_date ? filter_release_date(summaries, params.cutoff_date) : summaries

    // 2b. one genome per species. in: (group, [the RefSeq and GenBank summaries of the group], taxdb),
    // GenBank fill groups, min_relative_size and pilot_species ('' = off)
    zip_summary(dated
        .groupTuple(by: 1)                                           // ([catalogues], group, [summaries])
        .map { _catalogues, group, group_summaries -> tuple(group, group_summaries) }
        .combine(taxdb),
        params.genbank_fill_groups,
        params.min_relative_size ?: '',
        params.pilot_species ?: '')
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
    // hold-out size/min species/protected ranks, calibration and MetaEuk sizes, MetaEuk per phylum
    choose_representative(ncbi_table, size_filter.out.map { _group, genome -> genome }.collect(), taxdb,
        params.minbp, params.holdout_size, params.holdout_min_species, params.holdout_protect_ranks,
        params.calibration_size, params.metaeuk_size, params.metaeuk_per_phylum)

    // 2g. the legacy emit names: "training" = fna/remaining, "bitscore" = fna/training;
    // together they are every genome that is not in the hold-out
    database_genomes = choose_representative.out.training.flatten()
        .mix(choose_representative.out.bitscore.flatten())
        .ifEmpty { error "No genome left for the databases after the hold-out: every candidate went to it (--holdout_size ${params.holdout_size})" }
        .collect()
    split_clades(database_genomes, choose_representative.out.info)

    // 2h. in: one genome of the MetaEuk set; out: its annotated RefSeq proteome
    genome2faa(choose_representative.out.metaeuk.flatten())

    emit:
    ncbi_table = ncbi_table                                          // ncbi_table.txt
    zip_stats = zip_summary.out.stats                                // <group>_stats.tsv
    genetic_codes = genetic_codes.out                                // genetic_codes.tsv
    downloaded = download_genome.out                                 // (group, genome)
    sized = size_filter.out                                          // (group, genome)
    classification = choose_representative.out.info                  // species_classification.csv
    group_genomes = split_clades.out.protozoa                        // (group, [genomes])
        .mix(split_clades.out.plant)
        .mix(split_clades.out.fungi)
    metaeuk_proteomes = genome2faa.out                               // <accession>.faa
}
