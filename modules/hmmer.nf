// HMMER: the pressed PANTHER and calibrated libraries, searches, best hit per genome and marker,
// placement and completeness HMMs.

// The PANTHER family HMMs (books/*/hmmer.hmm, no subfamilies), concatenated and pressed.
// No storeDir: it reuses whatever is in the store without checking it, so a broken copy
// (a run killed while the GB-sized files were being moved, or the empty files of a stub
// run) breaks every hmmsearch of every later run. Rebuilt per run, cached by -resume
process concat_panther {
    label 'process_single'
    container 'biocontainers/hmmer:3.3--he1b5a44_0'

    input:
    path(panther)                                                                                              // PANTHER directory (books/)

    output:
    tuple path("all.hmm"), path("all.hmm.h3f"), path("all.hmm.h3i"), path("all.hmm.h3m"), path("all.hmm.h3p")  // (all.hmm + 4 hmmpress indices): every family HMM

    script:
    """
    export LC_ALL=C
    find -L ${panther}/books -mindepth 2 -maxdepth 2 -name hmmer.hmm | sort | xargs cat > all.hmm
    hmmpress all.hmm
    """

    stub:
    """
    touch all.hmm all.hmm.h3f all.hmm.h3i all.hmm.h3m all.hmm.h3p
    """
}

// Calibration genomes against every PANTHER family
process hmmsearch_training {
    tag "${clade}/${faa.baseName}"
    label 'process_medium'
    container 'biocontainers/hmmer:3.3--he1b5a44_0'

    input:
    tuple val(clade), path(faa), path(hmm), path(hmmindex1), path(hmmindex2), path(hmmindex3), path(hmmindex4)  // (clade, proteome of one calibration genome, all.hmm + 4 hmmpress indices)
    val(evalue)                                                                                                 // hmmsearch -E (--evalue)

    output:
    tuple val(clade), path("training_${faa.baseName}.tbl"), emit: tbl                                           // (clade, hmmsearch table)
    tuple val("${task.process}"), val('hmmer'), env('HMMER_VERSION'), topic: versions                           // tool version, set in the script

    script:
    """
    HMMER_VERSION=\$(hmmsearch -h | sed -n 2p | cut -d" " -f3)
    hmmsearch --cpu ${task.cpus} \\
        -o /dev/null \\
        --tblout training_${faa.baseName}.tbl \\
        -E ${evalue} \\
        ${hmm} ${faa}
    """

    stub:
    """
    HMMER_VERSION=\$(hmmsearch -h | sed -n 2p | cut -d" " -f3)
    touch training_${faa.baseName}.tbl
    """
}

process hmmpress_GA {
    tag "${clade}"
    label 'process_single'
    container 'biocontainers/hmmer:3.3--he1b5a44_0'

    input:
    tuple val(clade), path(hmm)                                                                                                       // (clade, all_GA.hmm)

    output:
    tuple val(clade), path(hmm, includeInputs: true), path("${hmm}.h3f"), path("${hmm}.h3i"), path("${hmm}.h3m"), path("${hmm}.h3p")  // (clade, all_GA.hmm + 4 hmmpress indices)

    script:
    """
    hmmpress ${hmm}
    """

    stub:
    """
    touch ${hmm}.h3f ${hmm}.h3i ${hmm}.h3m ${hmm}.h3p
    """
}

// Every genome of a sub-database against its calibrated library
process hmmsearch_GA {
    tag "${clade}/${faa.baseName}"
    label 'process_medium'
    container 'biocontainers/hmmer:3.3--he1b5a44_0'

    input:
    tuple val(clade), path(faa), path(hmm), path(hmmindex1), path(hmmindex2), path(hmmindex3), path(hmmindex4)  // (clade, proteome of one genome, all_GA.hmm + 4 hmmpress indices)

    output:
    tuple val(clade), val("${faa.baseName}"), path(faa), path("hmmerGA_${faa.baseName}.tbl")                    // (clade, accession, proteome, hmmsearch --cut_ga table)

    script:
    """
    hmmsearch --cpu ${task.cpus} \\
        -o /dev/null \\
        --cut_ga \\
        --tblout hmmerGA_${faa.baseName}.tbl \\
        ${hmm} ${faa}
    """

    stub:
    """
    printf '${faa.baseName}_1 - PTHR00001.orig.30.pir - 1e-30 100.0 0.0\\n${faa.baseName}_1 - PTHR00002.orig.30.pir - 1e-30 90.0 0.0\\n' > hmmerGA_${faa.baseName}.tbl
    """
}

// Highest-scoring protein of each genome for one marker
process group_top_hit {
    tag "${clade}/${hmmname}"
    label 'process_single'
    container 'biocontainers/pyfaidx:0.5.9.1--pyh864c0ab_1'

    input:
    tuple val(clade), val(hmmname), val(names), path(faa), path(tbl)  // (clade, marker, [accessions], [proteomes], [hmmsearch --cut_ga tables])

    output:
    tuple val(clade), val(hmmname), path("*_${hmmname}.faa")          // (clade, marker, [best protein of each genome])

    script:
    """
    mkdir fastas
    mv *.faa  fastas/
    mkdir tbls
    mv *.tbl tbls/

    for fasta in fastas/*; do
        NAME=\$(basename \$fasta | sed 's/.faa//')
        echo \$NAME
        top_hit.py  \\
            --faa  fastas/\${NAME}.faa \\
            --hits tbls/hmmerGA_\${NAME}.tbl \\
            --profile ${hmmname} \\
            --name \${NAME} \\
            \${NAME}_${hmmname}.faa
    done
    """

    stub:
    """
    for f in *.faa; do n=\$(basename \$f .faa); printf ">\$n\\nMKVL\\n" > \${n}_${hmmname}.faa; done
    """
}

process concat_faa {
    tag "${clade}/${hmmname}"
    label 'process_single'
    container 'biocontainers/python:3.9--1'

    input:
    tuple val(clade), val(hmmname), path(faas)              // (clade, marker, [best protein of each genome])

    output:
    tuple val(clade), val(hmmname), path("${hmmname}.faa")  // (clade, marker, <marker>.faa)

    script:
    """
    export LC_ALL=C
    cat *.faa > ${hmmname}.faa
    """
}

// Placement HMM of one trimmed marker alignment
process hmmbuild {
    tag "${clade}/${name}"
    label 'process_single'
    container 'biocontainers/hmmer:3.3--he1b5a44_0'

    input:
    tuple val(clade), val(name), path(fasta)  // (clade, marker, block)

    output:
    tuple val(clade), path("bb_${name}.hmm")  // (clade, HMM of the block)

    script:
    """
    hmmbuild bb_${name}.hmm ${fasta}
    """

    stub:
    """
    touch bb_${name}.hmm
    """
}

// The calibrated HMMs of the placement markers (<clade>_search.hmm)
process hmmfetch {
    tag "${clade}"
    label 'process_single'
    container 'biocontainers/hmmer:3.3--he1b5a44_0'

    input:
    tuple val(clade), path(hmm), path(hmmindex1), path(hmmindex2), path(hmmindex3), path(hmmindex4), path(profiles)  // (clade, all_GA.hmm + 4 hmmpress indices, profiles_<clade>.csv)

    output:
    tuple val(clade), path("backbone_profiles*")                                                                     // (clade, placement HMMs + indices)

    script:
    """
    awk '{print \$0".orig.30.pir"}' ${profiles} |cut -f 2 -d,> p.txt
    hmmfetch -f ${hmm} p.txt > backbone_profiles.hmm
    # remove  orig.30.pir
    sed -i '/^NAME/{s/.orig.30.pir//}' backbone_profiles.hmm
    hmmpress backbone_profiles.hmm
    """

    stub:
    """
    touch backbone_profiles.hmm backbone_profiles.hmm.h3f backbone_profiles.hmm.h3i backbone_profiles.hmm.h3m backbone_profiles.hmm.h3p
    """
}

// The calibrated HMMs of every profile in the SCMG table: what EukCC searches a MAG with
process hmmfetch_db {
    tag "${clade}"
    label 'process_single'
    container 'biocontainers/hmmer:3.3--he1b5a44_0'

    input:
    tuple val(clade), path(hmm), path(hmmindex1), path(hmmindex2), path(hmmindex3), path(hmmindex4), path(profiles)  // (clade, all_GA.hmm + 4 hmmpress indices, scmgs_<clade>.gz)

    output:
    tuple val(clade), path("db_PANTHER*")                                                                            // (clade, HMMs of every SCMG candidate + indices)

    script:
    """
    zcat ${profiles} | cut -f 2 -d, | \\
    sort | uniq | awk '{print \$0".orig.30.pir"}' >  p.txt
    hmmfetch -f ${hmm} p.txt > tmp.hmm

    # remove  orig.30.pir
    sed '/^NAME/{s/.orig.30.pir//}' tmp.hmm > db_PANTHER_GA.hmm
    rm tmp.hmm
    rm p.txt

    hmmpress db_PANTHER_GA.hmm
    """

    stub:
    """
    touch db_PANTHER_GA.hmm db_PANTHER_GA.hmm.h3f db_PANTHER_GA.hmm.h3i db_PANTHER_GA.hmm.h3m db_PANTHER_GA.hmm.h3p
    """
}
