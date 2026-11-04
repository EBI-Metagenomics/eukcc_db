process size_filter {
    tag "${genome.baseName}"
    label 'process_single'
    container 'biocontainers/python:3.9--1'

    input:
    tuple val(group), path(genome)                                            // (group, <accession>.fna)
    val(minbp)                                                                // minimum genome size in bp (--minbp)

    output:
    tuple val(group), path("${genome}", includeInputs: true), optional: true  // the same, only when the genome has >= minbp

    script:
    """
    BP=\$(grep -v ">" ${genome} | wc | awk '{print \$3-\$1}')
    echo "genome is \${BP} long"
    if [[ \$BP -lt "${minbp}" ]]; then
        echo "Deleting genome"
        rm ${genome}
    else
        echo "Keeping genome"
    fi
    """

    stub:
    """
    true
    """
}

// Global roles: hold-out ("testing", kept out of every database) and the RefSeq genomes
// whose annotated proteomes form the MetaEuk reference ("metaeuk")
process choose_representative {
    label 'process_single'
    container 'biocontainers/ete3:3.1.2'

    input:
    path(info)                                                   // ncbi_table: group,taxid,accession of every candidate
    path(genomes)                                                // every genome that passed size_filter
    path(taxdb, stageAs: 'taxa.sqlite')                          // taxa.sqlite
    val(minbp)                                                   // minimum genome size in bp (--minbp)
    val(holdout_size)                                            // RefSeq species in the hold-out (--holdout_size)
    val(holdout_min_species)                                     // species of every protected taxon kept out of the hold-out (--holdout_min_species)
    val(holdout_protect_ranks)                                   // ranks the hold-out never empties (--holdout_protect_ranks)
    val(calibration_size)                                        // species in the calibration set (--calibration_size)
    val(metaeuk_size)                                            // RefSeq proteomes in the MetaEuk set (--metaeuk_size)

    output:
    path("species_classification.csv"), emit: info               // group, accession, role, taxid, lineage
    path("fna/testing/*.fna"), emit: testing, optional: true     // hold-out
    path("fna/training/*.fna"), emit: bitscore, optional: true   // legacy name: the "training" role
    path("fna/remaining/*.fna"), emit: training, optional: true  // legacy name: the "remaining" role
    path("fna/metaeuk/*.fna"), emit: metaeuk, optional: true     // MetaEuk set

    script:
    """
    export HOME=\$PWD LC_ALL=C
    mkdir -p .etetoolkit && ln -sf \$PWD/taxa.sqlite .etetoolkit/taxa.sqlite
    mkdir genomes
    mv ${genomes} genomes/

    choose_representative.py --info ${info} \\
         --debug \\
         --testing \\
         --balance_level phylum \\
         --minbp ${minbp} \\
         --n_testing ${holdout_size} \\
         --testing_min_bucket ${holdout_min_species} \\
         --testing_protect ${[holdout_protect_ranks].flatten().join(',')} \\
         --n_training ${calibration_size} \\
         --n_metaeuk ${metaeuk_size} \\
         abc

    mkdir -p fna/remaining fna/testing fna/training fna/metaeuk
    awk '{print "cp -av genomes/"\$2".fna fna/"\$3"/" }' species_classification.csv | sh
    """

    stub:
    // first genome per group is the hold-out, second also goes to MetaEuk
    """
    mkdir -p genomes fna/remaining fna/testing fna/training fna/metaeuk
    mv ${genomes} genomes/
    awk -F, '{n[\$1]++; r = n[\$1] == 1 ? "testing" : (n[\$1] == 2 ? "training" : "remaining");
              print \$1"\\t"\$3"\\t"r"\\t"\$2"\\t1;2759"; if (n[\$1] == 2) print \$1"\\t"\$3"\\tmetaeuk\\t"\$2"\\t1;2759"}' ${info} \\
        > species_classification.csv
    awk '{print "cp -av genomes/"\$2".fna fna/"\$3"/" }' species_classification.csv | sh
    """
}

// Every genome except the hold-out, grouped by sub-database
process split_clades {
    label 'process_single'
    container 'biocontainers/python:3.9--1'

    input:
    path(genomes)                                                              // every genome outside the hold-out
    path(info)                                                                 // role of every genome

    output:
    tuple val("plant"), path("plant/*"), emit: plant, optional: true           // ("plant", [genomes])
    tuple val("fungi"), path("fungi/*"), emit: fungi, optional: true           // ("fungi", [genomes])
    tuple val("protozoa"), path("protozoa/*"), emit: protozoa, optional: true  // ("protozoa", [genomes])

    script:
    """
    mkdir genomes plant protozoa fungi
    mv *.fna genomes/

    grep -v testing ${info} | cut -f 1,2 | sort | uniq |\\
        awk '{print "cp -av genomes/"\$2".fna "\$1"/" }' | sh
    """
}

// Calibration set of one sub-database: ~calibration_size species balanced by phylum ("training")
process choose_representative_no_testing {
    tag "${clade}"
    label 'process_single'
    container 'biocontainers/ete3:3.1.2'

    input:
    tuple val(clade), path(genomes), path(info), path(taxdb, stageAs: 'taxa.sqlite')  // (clade, [proteomes of the clade], ncbi_table, taxa.sqlite)
    val(calibration_size)                                                             // species in the calibration set (--calibration_size)

    output:
    tuple val(clade), path("species_classification_${clade}.csv"), emit: info         // (clade, role of every genome)
    tuple val(clade), path("fna/training/*"), emit: bitscore                          // (clade, [proteomes of the calibration set])
    tuple val(clade), path("fna/remaining/*"), emit: training, optional: true         // (clade, [the other proteomes])

    script:
    """
    export HOME=\$PWD LC_ALL=C
    mkdir -p .etetoolkit && ln -sf \$PWD/taxa.sqlite .etetoolkit/taxa.sqlite
    mkdir genomes
    mv *.faa genomes

    choose_representative.py --info ${info} \\
         --debug \\
         --balance_level phylum \\
         --minbp 0 \\
         --n_training ${calibration_size} \\
         abc

    mkdir -p fna/remaining fna/testing fna/training fna/metaeuk
    awk '{print "cp -av genomes/"\$2".faa fna/"\$3"/" }' species_classification.csv | sh
    mv species_classification.csv species_classification_${clade}.csv
    """

    stub:
    // up to 2 genomes calibrate, the rest are "remaining"
    """
    mkdir -p genomes fna/training fna/remaining
    mv *.faa genomes
    i=0; for f in genomes/*.faa; do
        i=\$((i+1)); a=\$(basename \$f .faa)
        if [ \$i -le 2 ]; then r=training; else r=remaining; fi
        printf "${clade}\\t\$a\\t\$r\\t1\\t1\\n" >> species_classification_${clade}.csv
        cp -av \$f fna/\$r/
    done
    """
}

// Per-profile bit score cut-off that maximises single copies across taxonomic groups
process optimize_bitscore {
    tag "${clade}"
    label 'process_single'
    container 'biocontainers/ete3:3.1.2'

    input:
    tuple val(clade), path(tbls), path(taxinfo), path(taxdb, stageAs: 'taxa.sqlite')  // (clade, [hmmsearch tables of the calibration set], ncbi_table, taxa.sqlite)

    output:
    tuple val(clade), path("bitscore_cutoffs_${clade}.csv")                           // (clade, one bit score cut-off per profile)

    script:
    """
    export HOME=\$PWD LC_ALL=C
    mkdir -p .etetoolkit && ln -sf \$PWD/taxa.sqlite .etetoolkit/taxa.sqlite
    mkdir hmmer
    mv *.tbl hmmer/
    optimize_bitscore_groups.py --taxinfo ${taxinfo} \\
         --debug \\
        hmmer/*.tbl
    mv bitscore_cutoffs.csv bitscore_cutoffs_${clade}.csv
    """

    stub:
    """
    printf 'profile\\tcutoff\\tscore\\tmax_singletons\\tspecies\\nPTHR00001\\t20\\t1\\t1\\t1\\n' > bitscore_cutoffs_${clade}.csv
    """
}

// The calibrated cut-offs written as GA lines. Profiles without a cut-off are left out.
process insert_bitscore {
    tag "${clade}"
    label 'process_single'
    container 'biocontainers/python:3.9--1'

    input:
    tuple val(clade), path(bitscore), path(panther)  // (clade, bitscore_cutoffs_<clade>.csv, PANTHER directory)

    output:
    tuple val(clade), path("all_GA.hmm")             // (clade, the calibrated library)

    script:
    """
    modify_hmm.py --info ${bitscore} \\
        --panther ${panther}/books \\
        --debug all_GA.hmm
    """

    stub:
    """
    touch all_GA.hmm
    """
}

// Greedy selection of the placement markers
process choose_tree_profiles {
    tag "${clade}"
    label 'process_single'
    container 'biocontainers/python:3.9--1'

    input:
    tuple val(clade), path(hmmtbls)                                  // (clade, [hmmsearch --cut_ga table of every genome])

    output:
    tuple val(clade), path("profiles_${clade}.csv"), emit: profiles  // (clade, one clade,marker row per placement marker)
    path("find_tree_MG_${clade}.log"), emit: log                     // every iteration of the greedy selection

    script:
    // LC_ALL=C pins the glob order, i.e. the order find_tree_MG.py reads the tables,
    // which together with its seed(42) allows reproducible selection of the same markers
    """
    export LC_ALL=C
    mkdir hmmer
    mv *.tbl hmmer/
    find_tree_MG.py --debug  \\
        --clade ${clade} \\
        hmmer/*.tbl > find_tree_MG_${clade}.log 2>&1
    # the script writes a set: sort so the same choice gives the same file
    sort profiles.csv > profiles_${clade}.csv
    """

    stub:
    """
    printf '${clade},PTHR00001\\n${clade},PTHR00002\\n' > profiles_${clade}.csv
    touch find_tree_MG_${clade}.log
    """
}

// genome,profile for every profile with exactly one hit in a genome: the table EukCC
// draws its completeness markers from
process get_SCMG {
    tag "${clade}"
    label 'process_single'
    container 'biocontainers/python:3.9--1'

    input:
    tuple val(clade), path(hmmtbls)              // (clade, [hmmsearch --cut_ga table of every genome])

    output:
    tuple val(clade), path("scmgs_${clade}.gz")  // (clade, SCMG table: genome,profile rows)

    script:
    """
    export LC_ALL=C
    mkdir hmmer
    mv *.tbl hmmer/
    extract_SCMG.py --debug  \\
        hmmer/*.tbl
    # re-compress without the timestamp gzip writes in its header
    gzip -dc scmgs.gz | gzip -n > scmgs_${clade}.gz
    """

    stub:
    """
    echo "x,PTHR00001" | gzip > scmgs_${clade}.gz
    """
}

process reduce_taxinfo {
    tag "${clade}"
    label 'process_single'
    container 'biocontainers/python:3.9--1'

    input:
    tuple val(clade), path(hmmtbls), path(taxinfo)  // (clade, [hmmsearch --cut_ga tables], ncbi_table)

    output:
    tuple val(clade), path("taxinfo_${clade}.csv")  // (clade, taxid,accession of every genome searched)

    script:
    """
    ls *.tbl | cut -f 2-3 -d_ | cut -f 1-2 -d. > acc
    grep -f acc ${taxinfo} | cut -f 2-3 -d, > taxinfo_${clade}.csv
    """
}

process midpoint_root {
    tag "${clade}"
    label 'process_single'
    container 'biocontainers/ete3:3.1.2'

    input:
    tuple val(clade), path(tree)              // (clade, tree)

    output:
    tuple val(clade), path("rooted_${tree}")  // (clade, rooted tree)

    script:
    """
    root_tree.py ${tree} rooted_${tree}
    """

    stub:
    """
    cp ${tree} rooted_${tree}
    """
}

// iTOL annotation files: sub-database and phylum of each leaf on base, phylum to
// family on the others
process itol_annotations {
    tag "${clade}"
    label 'process_single'
    container 'biocontainers/ete3:3.1.2'

    input:
    tuple val(clade), path(tree), path(taxinfo), path(ncbi_table), path(taxdb, stageAs: 'taxa.sqlite')  // (clade, rooted tree, taxinfo_<clade>.csv, ncbi_table, taxa.sqlite)

    output:
    tuple val(clade), path("${clade}*")                                                                 // (clade, [tree and iTOL annotation files])

    script:
    def ranks = clade == 'base' ? 'subdb phylum' : 'phylum class order family'
    """
    cp ${tree} ${clade}_rooted.tre
    itol_annotations.py --tree ${tree} --taxinfo ${taxinfo} --groups ${ncbi_table} \\
        --taxdb taxa.sqlite --prefix ${clade} --levels ${ranks}
    """

    stub:
    """
    cp ${tree} ${clade}_rooted.tre
    touch ${clade}_labels.txt
    """
}
