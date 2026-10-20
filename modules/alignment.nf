// Gap filter and concatenation of the trimmed marker alignments.

// Leaves with too many gaps in the concatenation are dropped from every block
process filter_gappy_leaves {
    tag "${clade}"
    label 'process_single'
    container 'biocontainers/biopython:1.79'

    input:
    tuple val(clade), path(alns)                                    // (clade, [<marker>_trimmed.aln of every marker])
    val(max_leaf_gap_fraction)                                      // leaves with more gaps are dropped; '' = report only (--max_leaf_gap_fraction)

    output:
    tuple val(clade), path("filtered/*_trimmed.aln"), emit: alns    // (clade, [the same alignments without the gappy leaves]): the blocks
    tuple val(clade), path("leaf_gaps_${clade}.tsv"), emit: report  // (clade, fraction of gaps of each leaf)

    script:
    def threshold = max_leaf_gap_fraction ? "--max-gap-fraction ${max_leaf_gap_fraction}" : ''
    """
    filter_gappy_leaves.py --clade ${clade} ${threshold} --outdir filtered ${alns}
    """
}

// The supermatrix: marker blocks side by side, missing genomes filled with gaps
process horizontal {
    tag "${clade}"
    label 'process_single'
    container 'biocontainers/python:3.9--1'

    input:
    tuple val(clade), path(fasta)             // (clade, [blocks])

    output:
    tuple val(clade), path("${clade}.fasta")  // (clade, supermatrix)

    script:
    """
    horizontal.py ${clade}.fasta *.aln
    """
}
