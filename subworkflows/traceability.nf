// Step 5 of docs/workflow.md (5a-5b): build report, provenance and the release directory.

include { build_report; package_db } from '../modules/db'

workflow traceability {
    take:
    summaries        // (catalogue, group, <catalogue>_<group>.txt)
    zip_stats        // <group>_stats.tsv
    ncbi_table       // ncbi_table.txt
    downloaded       // (group, genome)
    sized            // (group, genome)
    classification   // species_classification.csv
    proteomes        // (group, <accession>.faa)
    searched         // "clade<TAB>accession" of every genome searched with --cut_ga
    report_files     // per-clade tables from build_db
    dbs              // db_<clade>/
    panther_db       // (all.hmm + 4 indices)
    taxdb            // taxa.sqlite
    taxonomy_info    // taxonomy.tsv
    genetic_codes    // genetic_codes.tsv

    main:
    // 5a. where every genome went, the counts of each sub-database, lineage coverage.
    // Each list of accessions is written to one file
    build_report(
        summaries.map { _catalogue, _group, summary -> summary }.collect(),
        zip_stats.collect(),
        ncbi_table,
        downloaded.map { _group, genome -> genome.baseName }.collectFile(name: 'downloaded.txt', newLine: true, sort: true),
        sized.map { _group, genome -> genome.baseName }.collectFile(name: 'sized.txt', newLine: true, sort: true),
        classification,
        proteomes.map { _group, proteome -> proteome.baseName }.collectFile(name: 'predicted.txt', newLine: true, sort: true),
        searched.collectFile(name: 'searched.txt', newLine: true, sort: true),
        report_files.collect(),
        panther_db.map { all_hmm, _h3f, _h3i, _h3m, _h3p -> all_hmm },
        taxdb,
        genetic_codes,
        params.minbp
    )

    // "tool: version" of every tool, read from the containers (topic channel)
    versions = channel.topic('versions')
        .map { _process, tool, version -> "${tool}: ${version}".toString() }
        .unique()
        .collectFile(name: 'versions.yml', newLine: true, sort: true)

    // tables copied to eukcc2_db_<name>/provenance/ (without the supermatrices and SCMG tables)
    provenance = report_files
        .filter { report_file -> !report_file.name.endsWith('.fasta') && !report_file.name.endsWith('.gz') }
        .mix(build_report.out.exclusions, build_report.out.summary, build_report.out.coverage)
        .mix(ncbi_table, genetic_codes, taxonomy_info)
        .mix(zip_stats)
        .collect()
    // release name: --name, otherwise YYYY_MM of --cutoff_date, otherwise of the day of the run
    def name = params.name ?: (params.cutoff_date ?: params.run_stamp).toString().take(7).replace('-', '_')
    // 5b. in: [db_<clade>/], taxdb, [provenance tables], species_classification.csv, versions.yml, release name,
    // run parameters, [assembly_summary.txt]
    package_db(dbs.collect(), taxdb, provenance, classification, versions, name, params,
               summaries.map { _catalogue, _group, summary -> summary }.collect())
}
