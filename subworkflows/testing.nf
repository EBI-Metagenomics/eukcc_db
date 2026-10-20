// testing: MAGs simulated from the hold-out genomes of a build, EukCC on its database

include { download_genome }                                                from '../modules/ncbi'
include { simulate_mag; eukcc_single; sum_quality; compare_quality }       from '../modules/testing'

workflow testing {
    main:
    if (!params.db) {
        error "--mode testing needs --db <eukcc2_db_name>"
    }
    db_dir = file(params.db, type: 'dir', checkIfExists: true)
    classification = file(params.classification ?: "${params.db}/provenance/species_classification.csv", checkIfExists: true)

    // the hold-out genomes of the build, from the same cache the build downloaded them to.
    // species_classification.csv columns: group, accession, role, taxid, lineage
    holdout = channel.of(classification)
        .splitCsv(sep: '\t')
        .filter { row -> row[2] == 'testing' }
        .map { row -> tuple(row[0], row[3], row[1]) }               // (group, taxid, accession)
    holdout_genomes = download_genome(holdout).map { _group, genome -> genome }

    // contaminant: the genome itself ("self") and/or the next hold-out genome
    // in accession order ("other", thesis tables A1/A2)
    modes = params.contamination_source == 'both' ? ['self', 'other'] : [params.contamination_source]
    // (MAG name, genome, contaminant, completeness, contamination): every combination of
    // genome, mode, completeness, contamination and replicate, then a fixed draw of eval_n_mags
    simulations = holdout_genomes
        .toSortedList { first, second -> first.name <=> second.name }
        .flatMap { sorted_genomes ->
            def grid = []
            sorted_genomes.eachWithIndex { genome, index ->
                modes.each { mode ->
                    def contaminant = mode == 'self' ? genome : sorted_genomes[(index + 1) % sorted_genomes.size()]
                    params.eval_completeness.each { completeness ->
                        params.eval_contamination.each { contamination ->
                            (1..params.eval_n_repeats).each { replicate ->
                                grid << tuple("${genome.baseName}_${completeness}_${contamination}_${replicate}_${mode}",
                                              genome, contaminant, completeness, contamination)
                            }
                        }
                    }
                }
            }
            // a fixed draw of eval_n_mags combinations: seeded shuffle of the sorted grid
            Collections.shuffle(grid, new Random(params.eval_seed as long))
            grid.take(params.eval_n_mags as int)
        }
    // out: (MAG name, MAG fasta, true completeness/contamination)
    simulate_mag(simulations)
    // in: (MAG name, MAG fasta, truth, database); out: (MAG name, truth, EukCC estimate)
    eukcc_single(simulate_mag.out.combine(channel.of(db_dir)))
    // in: (MAG name, truth, estimate); out: one row comparing them
    sum_quality(eukcc_single.out)
    // in: [every row], species_classification.csv
    compare_quality(sum_quality.out.collect(), classification)
}
