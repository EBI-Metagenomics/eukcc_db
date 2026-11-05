// Sub-database directory, release package and build report.

process db {
    tag "${clade}"
    label 'process_single'
    container 'biocontainers/python:3.9--1'

    input:
    tuple val(clade),
          path(scmgs),           // SCMG table: genome,profile rows
          path(treealn),         // supermatrix
          path(tree),            // rooted tree
          path(hmmdb),           // SCMG HMMs + indices
          path(placement_hmm),   // placement HMMs + indices
          path(bb_aln),          // [blocks]
          path(bb_profiles),     // [block HMMs]
          path(taxinfo)          // taxid,accession of every genome searched
    path(metaeukdb)              // [MetaEuk files]
    val(model)                   // substitution model of the tree
    val(panther_version)         // PANTHER version (--panther_version), written to db_info.json

    output:
    path("db_${clade}")          // db_<clade>/ with its db_info.json

    script:
    """
    create_db.py --debug \\
        --clade ${clade} \\
        --backbone_tree ${tree} \\
        --backbone_alignment ${treealn} \\
        --scmgs  ${scmgs} \\
        --hmm_db db_PANTHER_GA.hmm* \\
        --placement_hmm backbone_profiles.hmm* \\
        --iqtree ${model} \\
        --panther  ${panther_version} \\
        --places_aln ${bb_aln} \\
        --places_hmm ${bb_profiles} \\
        --taxinfo ${taxinfo} \\
        --metaeuk  microeuks_metaeuk.db*
    # empty text library (~2 GB): hmmfetch reads just the pressed .h3* files
    truncate -s 0 db_${clade}/scmgs/${clade}_PANTHERv${panther_version}.hmm
    """
}

// Where every candidate genome went, the counts of each sub-database, and which
// groups (fungi, plant, protozoa) and phyla each sub-database covers
process build_report {
    label 'process_single'
    container 'biocontainers/python:3.9--1'

    input:
    path(summaries, stageAs: 'summaries/*')                      // the assembly_summary.txt files
    path(zip_stats, stageAs: 'zip_stats/*')                      // <group>_stats.tsv of zip_summary
    path(candidates, stageAs: 'ncbi_table.txt')                  // group,taxid,accession of every candidate
    path(downloaded, stageAs: 'downloaded.txt')                  // accessions downloaded
    path(sized, stageAs: 'sized.txt')                            // accessions that passed the size filter
    path(classification, stageAs: 'species_classification.csv')  // role of every genome
    path(predicted, stageAs: 'predicted.txt')                    // accessions with a predicted proteome
    path(searched, stageAs: 'searched.txt')                      // clade<TAB>accession searched with --cut_ga
    path(per_subdb, stageAs: 'subdb/*')                          // per-clade tables of build_db
    path(panther_hmm)                                            // all.hmm (counts the PANTHER families)
    path(taxdb, stageAs: 'taxa.sqlite')                          // taxa.sqlite
    path(codes, stageAs: 'genetic_codes.tsv')                    // accession, taxid, nuclear genetic code
    val(minbp)                                                   // minimum genome size in bp (--minbp)
    path(exclude_genomes, stageAs: 'exclude_genomes.tsv')        // accessions dropped by zip_summary; [] = none (--exclude_genomes)

    output:
    path("exclusions.tsv"), emit: exclusions                     // where every candidate was removed and why
    path("build_summary.tsv"), emit: summary                     // count after each step, per clade
    path("coverage.tsv"), emit: coverage                         // groups and phyla each clade covers

    script:
    def exclude = exclude_genomes ? '--exclude-genomes exclude_genomes.tsv' : ''
    """
    build_report.py \\
        --summaries summaries/* \\
        --zip-stats zip_stats/* \\
        --candidates ncbi_table.txt \\
        --downloaded downloaded.txt \\
        --sized sized.txt \\
        --classification species_classification.csv \\
        --predicted predicted.txt \\
        --searched searched.txt \\
        --subdb-dir subdb \\
        --panther-families \$(grep -c '^NAME' ${panther_hmm}) \\
        --minbp ${minbp} \\
        --taxdb taxa.sqlite \\
        --genetic-codes genetic_codes.tsv ${exclude}
    """
}

// The release directory: the four sub-databases, the taxonomy EukCC must use, the
// provenance of the build and a README.md (version YYYY_MM, sources, parameters).
// Fails the run if the structural checks do not pass.
process package_db {
    label 'process_single'
    container 'biocontainers/python:3.9--1'

    input:
    path(dbs)                                                    // [db_<clade>/]
    path(taxdb, stageAs: 'taxa.sqlite')                          // taxa.sqlite
    path(provenance, stageAs: 'prov/*')                          // provenance tables
    path(classification, stageAs: 'species_classification.csv')  // role of every genome
    path(versions, stageAs: 'versions.yml')                      // tool versions
    val(name)                                                    // release name (--name)
    val(run_params)                                              // every parameter of the run, written to params.json
    path(summaries, stageAs: 'summaries/*')                      // the assembly_summary.txt files (md5 in README.md)

    output:
    path("eukcc2_db_${name}")                                    // the release directory

    script:
    // single-quoted for the shell: escape the quotes inside the JSON
    def json = groovy.json.JsonOutput.toJson(run_params).replace("'", "'\\''")
    """
    D=eukcc2_db_${name}
    mkdir -p \$D/etetoolkit \$D/provenance
    cp -rL ${dbs} \$D/
    cp -L taxa.sqlite \$D/etetoolkit/taxa.sqlite
    cp -L prov/* species_classification.csv versions.yml \$D/provenance/

    echo '${json}' | python -m json.tool > \$D/provenance/params.json
    printf 'pipeline: %s\\ncommit: %s\\nnextflow: %s\\nrun: %s\\n' \\
        '${workflow.manifest.name} ${workflow.manifest.version}' \\
        '${workflow.commitId ?: 'local'}' \\
        '${nextflow.version}' \\
        '${workflow.runName} ${workflow.start}' \\
        > \$D/provenance/run.txt

    check_db.py \$D --classification species_classification.csv
    release_readme.py --db \$D --summaries summaries/* --output \$D/README.md
    """

    stub:
    """
    D=eukcc2_db_${name}
    mkdir -p \$D/etetoolkit \$D/provenance
    cp -rL ${dbs} \$D/
    cp -L prov/* species_classification.csv versions.yml \$D/provenance/
    touch \$D/README.md
    """
}
