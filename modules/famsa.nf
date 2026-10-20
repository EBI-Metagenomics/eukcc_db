process famsa {
    tag "${clade}/${name}"
    label 'process_low'
    container 'biocontainers/famsa:1.6.1--hc9558a2_0'

    input:
    tuple val(clade), val(name), path(fasta)                                           // (clade, marker, <marker>.faa)

    output:
    tuple val(clade), val(name), path("${name}.aln"), emit: aln                        // (clade, marker, <marker>.aln)
    tuple val("${task.process}"), val('famsa'), env('FAMSA_VERSION'), topic: versions  // tool version, set in the script

    script:
    """
    FAMSA_VERSION=\$(famsa 2>&1 | grep -m1 -o "ver\\. [^ ]*" | cut -d" " -f2)
    famsa -t ${task.cpus} ${fasta} ${name}.aln
    """

    stub:
    """
    FAMSA_VERSION=\$(famsa 2>&1 | grep -m1 -o "ver\\. [^ ]*" | cut -d" " -f2)
    cp ${fasta} ${name}.aln
    """
}
