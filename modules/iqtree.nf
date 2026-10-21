// Substitution model chosen once, on one sub-database (base by default), and used for
// every tree. -m MF selects the model without building the tree that MFP would add.
// -safe: the base supermatrix (~1,700 leaves) underflows the fast likelihood kernel
// ("Numerical underflow (lh-derivative)"); the safe kernel is slower but does not.
process modelfinder {
    tag "${clade}"
    label 'process_high'
    label 'process_long'
    container 'biocontainers/iqtree:2.3.6--hdbdd923_0'

    input:
    tuple val(clade), path(fasta)             // (clade, supermatrix) of --modelfinder_clade
    val(seed)                                 // IQ-TREE -seed (--seed)

    output:
    path("model.txt"), emit: model            // the best-fit model
    path("modelfinder.iqtree"), emit: report  // ModelFinder report

    script:
    """
    iqtree -s ${fasta} \\
        -m MF -T ${task.cpus} \\
        -safe \\
        --seqtype AA \\
        -seed ${seed} \\
        --prefix modelfinder
    sed -n 's/^Best-fit model: \\([^ ]*\\) chosen.*/\\1/p' modelfinder.log > model.txt
    test -s model.txt
    """

    stub:
    """
    echo LG+R7 > model.txt
    touch modelfinder.iqtree
    """
}

process iqtree {
    tag "${clade}"
    label 'process_high'
    label 'process_long'
    container 'biocontainers/iqtree:2.3.6--hdbdd923_0'

    input:
    tuple val(clade), path(fasta), val(model)                                            // (clade, supermatrix, model)
    val(seed)                                                                            // IQ-TREE -seed (--seed)

    output:
    tuple val(clade), path("*.treefile"), emit: tree                                     // (clade, tree)
    tuple val("${task.process}"), val('iqtree'), env('IQTREE_VERSION'), topic: versions  // tool version, set in the script

    script:
    """
    IQTREE_VERSION=\$(iqtree --version | sed -n 1p | grep -o "version [^ ]*" | cut -d" " -f2)
    iqtree -s  ${fasta} \\
        --seqtype AA \\
        -m ${model} \\
        -safe \\
        -T ${task.cpus} \\
        -seed ${seed} \\
        --prefix backbone
    """

    stub:
    """
    IQTREE_VERSION=\$(iqtree --version | sed -n 1p | grep -o "version [^ ]*" | cut -d" " -f2)
    echo "(a,b);" > backbone.treefile
    """
}
