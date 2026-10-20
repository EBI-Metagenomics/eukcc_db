process trimal_named {
    tag "${clade}/${name}"
    label 'process_single'
    container 'biocontainers/trimal:1.4.1--h4ac6f70_9'

    input:
    tuple val(clade), val(name), path(fasta)                                             // (clade, marker, <marker>.aln)

    output:
    tuple val(clade), val(name), path("${name}_trimmed.aln"), emit: aln                  // (clade, marker, <marker>_trimmed.aln)
    tuple val("${task.process}"), val('trimal'), env('TRIMAL_VERSION'), topic: versions  // tool version, set in the script

    script:
    """
    TRIMAL_VERSION=\$(trimal --version | grep -o "v[^ ]*")
    trimal -in ${fasta} \\
        -out ${name}_trimmed.aln \\
        -gappyout
    """

    stub:
    """
    TRIMAL_VERSION=\$(trimal --version | grep -o "v[^ ]*")
    cp ${fasta} ${name}_trimmed.aln
    """
}
