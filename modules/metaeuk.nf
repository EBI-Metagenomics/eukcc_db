// Protein prediction with MetaEuk against the annotated RefSeq proteomes.

process metaeuk_db {
    label 'process_low'
    container 'biocontainers/metaeuk:4.a0f584d--h95f258a_2'

    input:
    path(fasta)                                                                            // the proteomes of the MetaEuk set, in one file

    output:
    path("microeuks_metaeuk.db*"), emit: db                                                // [microeuks_metaeuk.db*]
    tuple val("${task.process}"), val('metaeuk'), env('METAEUK_VERSION'), topic: versions  // tool version, set in the script

    script:
    """
    METAEUK_VERSION=\$(metaeuk version)
    metaeuk createdb \\
        ${fasta} microeuks_metaeuk.db
    """

    stub:
    """
    METAEUK_VERSION=\$(metaeuk version)
    touch microeuks_metaeuk.db microeuks_metaeuk.db.dbtype microeuks_metaeuk.db.index
    """
}

// Proteins of one genome, translated with its nuclear genetic code
process metaeuk_pred {
    tag "${clade}/${fasta.baseName} (code ${code})"
    label 'process_medium'
    container 'biocontainers/metaeuk:4.a0f584d--h95f258a_2'

    input:
    tuple val(clade), path(fasta), val(code), path(db)  // (group, <accession>.fna, genetic code, [microeuks_metaeuk.db*])

    output:
    tuple val(clade), path("${fasta.baseName}.faa")     // (group, <accession>.faa)

    script:
    def acc = fasta.baseName
    """
    metaeuk easy-predict \\
        --protein 1 \\
        --translation-table ${code} \\
        --threads ${task.cpus} \\
        ${fasta} \\
        microeuks_metaeuk.db \\
        ${acc}_metaeuk \\
        tmp_folder
    # metaeuk names are really long: <accession>_1, <accession>_2, ...
    awk -v p=${acc} '/^>/ {print ">" p "_" (++n); next} {print}' ${acc}_metaeuk.fas > ${acc}.faa
    rm ${acc}_metaeuk.fas
    rm -rf tmp_folder
    """

    stub:
    """
    printf '>${fasta.baseName}_1\\nMKVL\\n' > ${fasta.baseName}.faa
    """
}
