// Candidate genomes from NCBI: assembly summaries, one genome per species, downloads.

process get_group_summary {
    tag "${db}/${group}"
    label 'process_single'
    container 'biocontainers/python:3.9--1'

    input:
    tuple val(group), val(db)                              // (group, catalogue): e.g. (fungi, refseq)
    val(ncbi_url)                                          // NCBI FTP root (--ncbi_url)

    output:
    tuple val(db), val(group), path("${db}_${group}.txt")  // (catalogue, group, assembly_summary.txt)

    script:
    """
    wget -q -O ${db}_${group}.txt ${ncbi_url}/genomes/${db}/${group}/assembly_summary.txt
    """
}

// Optional cut-off on seq_rel_date (column 15 of assembly_summary.txt)
process filter_release_date {
    tag "${db}/${group}"
    label 'process_single'
    container 'biocontainers/python:3.9--1'

    input:
    tuple val(db), val(group), path(summary)                     // (catalogue, group, assembly_summary.txt)
    val(cutoff_date)                                             // YYYY-MM-DD (--cutoff_date)

    output:
    tuple val(db), val(group), path("dated_${db}_${group}.txt")  // the same, only rows released up to --cutoff_date

    script:
    def cutoff = cutoff_date.toString().replace('-', '/')
    """
    awk -F'\\t' -v d='${cutoff}' '/^#/ || \$15 <= d' ${summary} > dated_${db}_${group}.txt
    """
}

// One genome per species of a group (RefSeq first), land plants and animals excluded.
// In genbank_fill_groups, a GenBank genome enters only for a phylum without RefSeq.
// Genomes far smaller than their relatives (partial assemblies) are dropped.
// The accessions of --exclude_genomes are dropped before one genome per species.
// Writes the count left after each filter.
process zip_summary {
    tag "${group}"
    label 'process_single'
    container 'biocontainers/ete3:3.1.2'

    input:
    tuple val(group), path(summaries), path(taxdb, stageAs: 'taxa.sqlite')           // (group, [RefSeq and GenBank summaries of the group], taxa.sqlite)
    val(genbank_fill_groups)                                                         // groups where GenBank fills only the phyla without RefSeq (--genbank_fill_groups)
    val(min_relative_size)                                                           // minimum size as a fraction of the genus median; '' = off (--min_relative_size)
    val(pilot_species)                                                               // species kept per division; '' = all (--pilot_species)
    path(exclude_genomes, stageAs: 'exclude/*')                                      // accession<TAB>reason never to keep; [] = none (--exclude_genomes)

    output:
    tuple val(group), path("${group}_combined.csv"), emit: table                     // (group, <group>_combined.csv: group,taxid,accession of every candidate)
    path("${group}_stats.tsv"), emit: stats                                          // count left after each filter
    tuple val("${task.process}"), val('ete3'), env('ETE3_VERSION'), topic: versions  // tool version, set in the script

    script:
    // ete3's NCBITaxa() reads $HOME/.etetoolkit/taxa.sqlite: point HOME at the run's copy
    def fill = group in genbank_fill_groups ? '--genbank_fill' : ''
    def relsize = min_relative_size ? "--min_relative_size ${min_relative_size}" : ''
    def pilot = pilot_species ? "--pilot ${pilot_species}" : ''
    def exclude = exclude_genomes ? "--exclude ${exclude_genomes}" : ''
    """
    ETE3_VERSION=\$(python -c "import ete3; print(ete3.__version__)")
    export HOME=\$PWD LC_ALL=C
    mkdir -p .etetoolkit && ln -sf \$PWD/taxa.sqlite .etetoolkit/taxa.sqlite
    zip_ncbi.py --label ${group} --stats ${group}_stats.tsv ${fill} ${relsize} ${pilot} ${exclude} ${group}_combined.csv *.txt
    """

    stub:
    """
    ETE3_VERSION=\$(python -c "import ete3; print(ete3.__version__)")
    # the accessions of --exclude_genomes (exact match) are out before the first 4 are taken
    { cat exclude/* 2>/dev/null || true; } > excluded.tsv
    awk -F'\\t' -v g=${group} 'FILENAME == "excluded.tsv" {x[\$1]; next}
        !/^#/ && !(\$1 in x) && (\$5 == "reference genome" || \$5 == "representative genome") {print g","\$7","\$1}' excluded.tsv *.txt | head -4 > ${group}_combined.csv
    n=\$(cat ${group}_combined.csv | wc -l)
    printf '${group}\\trows\\t%s\\n${group}\\treference_full\\t%s\\n${group}\\tspecies\\t%s\\n${group}\\texcluded_lineage\\t0\\n${group}\\tcandidates\\t%s\\n' \$n \$n \$n \$n > ${group}_stats.tsv
    """
}

process download_genome {
    tag "${accession}"
    label 'process_single'
    container 'biocontainers/python:3.9--1'
    storeDir "${workflow.stubRun ? "${workflow.workDir}/stub_cache" : params.cachedir}/genomes"
    maxForks params.ncbi_maxforks

    input:
    tuple val(group), val(taxid), val(accession)  // (group, taxid, accession) of one candidate

    output:
    tuple val(group), path("${accession}.fna")    // (group, <accession>.fna)

    script:
    """
    mkdir genome
    download_from_ncbi.py ${accession} genome/
    mv genome/${accession}* ${accession}_dirty.fna.gz
    gunzip *.gz

    # rename fasta headers
    rename_multifasta_prefix.py -f ${accession}_dirty.fna -p ${accession} > ${accession}.fna
    rm ${accession}_dirty.fna
    """

    stub:
    """
    printf '>${accession}_1\\nACGT\\n' > ${accession}.fna
    """
}

// Annotated RefSeq proteome of a genome, for the MetaEuk reference database
process genome2faa {
    tag "${fna.baseName}"
    label 'process_single'
    container 'biocontainers/python:3.9--1'
    maxForks params.ncbi_maxforks

    input:
    path(fna)                    // <accession>.fna of one MetaEuk set genome (only its name is used)

    output:
    path("${fna.baseName}.faa")  // annotated RefSeq proteome

    script:
    def accession = fna.baseName
    """
    mkdir genome
    download_from_ncbi.py --prot ${accession} genome/
    mv genome/${accession}* ${accession}_dirty.faa.gz
    gunzip *.gz

    # rename fasta headers
    rename_multifasta_prefix.py -f ${accession}_dirty.faa -p ${accession} > ${accession}.faa
    rm ${accession}_dirty.faa
    """

    stub:
    """
    printf '>${fna.baseName}_1\\nMKV\\n' > ${fna.baseName}.faa
    """
}
