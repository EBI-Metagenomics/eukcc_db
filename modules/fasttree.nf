// Used only with --iqtree false
process fasttree {
    tag "${clade}"
    label 'process_low'
    container 'biocontainers/fasttree:2.1.10--h779adbc_6'

    input:
    tuple val(clade), path(fasta)                                                            // (clade, supermatrix)

    output:
    tuple val(clade), path("*.tree"), emit: tree                                             // (clade, tree)
    tuple val("${task.process}"), val('fasttree'), env('FASTTREE_VERSION'), topic: versions  // tool version, set in the script

    script:
    """
    FASTTREE_VERSION=\$( (fasttree 2>&1 || true) | sed -n 1p | cut -d" " -f3 )
    fasttree -log tree.log ${fasta} > fasttree.tree
    """

    stub:
    """
    FASTTREE_VERSION=\$( (fasttree 2>&1 || true) | sed -n 1p | cut -d" " -f3 )
    echo "(a,b);" > fasttree.tree
    """
}
