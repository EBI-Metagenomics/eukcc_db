// Evaluation: MAGs simulated from the hold-out genomes, EukCC on the new
// database, estimate against truth.

process simulate_mag {
    tag "${prefix}"
    label 'process_single'
    container 'biocontainers/python:3.9--1'

    input:
    tuple val(prefix), path(genome), path(contaminant, stageAs: 'contaminant/*'), val(compl), val(cont)  // (MAG name, genome, contaminant genome, completeness %, contamination %)

    output:
    tuple val(prefix), path("${prefix}_mag.fasta"), path("${prefix}_stats.tsv")                          // (MAG name, MAG fasta, true completeness and contamination)

    script:
    // a seed derived from the MAG name: the same grid gives the same MAGs
    def seed = prefix.hashCode() & 0x7fffffff
    """
    simulate_mag.py  --genome ${genome} \\
                --fraction ${compl}/${cont} \\
                --debug \\
                --prefix ${prefix} \\
                --contaminate ${contaminant} \\
                --seed ${seed}
    """

    stub:
    """
    printf '>reference_1\\nACGT\\n' > ${prefix}_mag.fasta
    printf 'simulation_genome\\tsimulation_contaminant\\tsimulated_completeness\\tsimulated_contamination\\tsimulation_prefix\\n${genome.name}\\t${contaminant.name}\\t${compl}\\t${cont}\\t${prefix}\\n' > ${prefix}_stats.tsv
    """
}

process eukcc_single {
    tag "${prefix}"
    label 'process_medium'
    container 'biocontainers/eukcc:2.2.0--pyhdfd78af_0'

    input:
    tuple val(prefix), path(mag), path(stats), path(db)          // (MAG name, MAG fasta, true values, database directory)

    output:
    tuple val(prefix), path(stats), path("eukcc_${prefix}.tsv")  // (MAG name, true values, EukCC estimate)

    script:
    // a failed estimate is reported as NA, so it counts as a failure in the report
    """
    eukcc single --db ${db} --out out --threads ${task.cpus} ${mag} || true
    if [ -s out/eukcc.tsv ]; then
        cp out/eukcc.tsv eukcc_${prefix}.tsv
    else
        printf 'fasta\\tcompleteness\\tcontamination\\tncbi_lng\\n${mag}\\tNA\\tNA\\tNA\\n' > eukcc_${prefix}.tsv
    fi
    """

    stub:
    """
    printf 'fasta\\tcompleteness\\tcontamination\\tncbi_lng\\n${mag}\\t65.0\\t3.0\\t1-2759\\n' > eukcc_${prefix}.tsv
    """
}

process sum_quality {
    tag "${prefix}"
    label 'process_single'
    container 'biocontainers/python:3.9--1'

    input:
    tuple val(prefix), path(simulation), path(estimation)  // (MAG name, true values, EukCC estimate)

    output:
    path("summed_${prefix}.tsv")                           // one row with both

    script:
    // csv.DictWriter ends lines with \r\n: strip them
    """
    quality_summary.py -e ${estimation} \\
        -s ${simulation}
    tr -d '\\r' < summed_*.tsv > tmp.tsv
    rm summed_*.tsv
    mv tmp.tsv summed_${prefix}.tsv
    """
}

process compare_quality {
    label 'process_single'
    container 'biocontainers/python:3.9--1'

    input:
    path(tables, stageAs: 'summed/*')  // every summed_<name>.tsv
    path(classification)               // role of every genome

    output:
    path("quality_*.tsv")              // quality_per_mag.tsv, quality_summary.tsv

    script:
    """
    compare_quality.py --classification ${classification} summed/*.tsv
    """
}
