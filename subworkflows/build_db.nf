// Step 4 of docs/workflow.md (4a-4r): one sub-database from the proteomes of its genomes.
// Runs once per clade: base, fungi, plant, protozoa.

include { choose_representative_no_testing as genome_group;
          optimize_bitscore; insert_bitscore; choose_tree_profiles;
          get_SCMG; reduce_taxinfo; midpoint_root; itol_annotations }      from '../modules/general'
include { hmmsearch_training; hmmpress_GA; hmmsearch_GA; group_top_hit;
          concat_faa; hmmbuild; hmmfetch; hmmfetch_db }                    from '../modules/hmmer'
include { famsa }                                                          from '../modules/famsa'
include { trimal_named }                                                   from '../modules/trimal'
include { filter_gappy_leaves; horizontal }                                from '../modules/alignment'
include { modelfinder; iqtree }                                            from '../modules/iqtree'
include { fasttree }                                                       from '../modules/fasttree'
include { db }                                                             from '../modules/db'

workflow build_db {
    take:
    clade_proteomes // (clade, [proteome .faa of every genome of the clade])
    ncbi_table      // ncbi_table.txt: group,taxid,accession of every candidate
    metaeuk_files   // [microeuks_metaeuk.db*]
    taxdb           // taxa.sqlite
    panther_db      // (all.hmm, all.hmm.h3f, .h3i, .h3m, .h3p): every PANTHER family
    panther         // PANTHER directory (books/)

    main:
    // calibration sample (~calibration_size species, balanced by phylum) and the rest.
    genome_group(clade_proteomes.combine(ncbi_table).combine(taxdb), params.calibration_size)
    calibration_proteomes = genome_group.out.bitscore.transpose()   // (clade, proteome)
    other_proteomes = genome_group.out.training.transpose()         // (clade, proteome)

    // per-profile bit score cut-offs
    // in: (clade, calibration proteome, all.hmm + 4 indices)
    hmmsearch_training(calibration_proteomes.combine(panther_db), params.evalue)
    // in: (clade, [hmmsearch tables of the calibration set], ncbi_table, taxdb)
    optimize_bitscore(hmmsearch_training.out.tbl.groupTuple().combine(ncbi_table).combine(taxdb))
    // in: (clade, bitscore_cutoffs_<clade>.csv, PANTHER directory)
    insert_bitscore(optimize_bitscore.out.combine(panther))
    // in: (clade, all_GA.hmm); out: (clade, all_GA.hmm + 4 indices) = the calibrated library
    calibrated = hmmpress_GA(insert_bitscore.out)

    // every genome of the sub-database against its calibrated library
    // in: (clade, proteome, all_GA.hmm + 4 indices)
    hmmsearch_GA(other_proteomes.mix(calibration_proteomes).combine(calibrated, by: 0))
    // out of hmmsearch_GA: (clade, accession, proteome, hmmerGA_<accession>.tbl)
    ga_tables = hmmsearch_GA.out
        .map { clade, _accession, _proteome, ga_table -> tuple(clade, ga_table) }
        .groupTuple()                                                // (clade, [GA tables])

    // placement markers, SCMG table, taxinfo
    // in: (clade, [GA tables])
    choose_tree_profiles(ga_tables)
    get_SCMG(ga_tables)
    // in: (clade, [GA tables], ncbi_table)
    reduce_taxinfo(ga_tables.combine(ncbi_table))

    // profiles_<clade>.csv has one "clade,marker" row per placement marker
    markers = choose_tree_profiles.out.profiles.flatMap { _clade, profiles_csv ->
        profiles_csv.splitCsv().collect { row -> tuple(row[0], row[1]) }
    }
    // every genome for every marker: (clade, marker, accession, proteome, GA table)
    marker_hits = markers.combine(hmmsearch_GA.out, by: 0)

    // best protein per genome and marker -> FAMSA -> trimAl
    // in: (clade, marker, [accessions], [proteomes], [GA tables])
    group_top_hit(marker_hits.groupTuple(by: [0, 1]))
    // in: (clade, marker, [best protein of each genome]); out: (clade, marker, <marker>.faa)
    concat_faa(group_top_hit.out)
    // in: (clade, marker, <marker>.faa); out: (clade, marker, <marker>.aln)
    famsa(concat_faa.out)
    // in: (clade, marker, <marker>.aln); out: (clade, marker, <marker>_trimmed.aln)
    trimal_named(famsa.out.aln)

    // drop leaves with too many gaps, then one HMM per block and the concatenation
    // in: (clade, [trimmed alignment of every marker])
    filter_gappy_leaves(trimal_named.out.aln
        .map { clade, _marker, trimmed_aln -> tuple(clade, trimmed_aln) }
        .groupTuple(),
        params.max_leaf_gap_fraction ?: '')
    // the filtered alignments are the blocks: (clade, marker, <marker>_trimmed.aln)
    blocks = filter_gappy_leaves.out.alns.transpose()
        .map { clade, block -> tuple(clade, block.name - '_trimmed.aln', block) }
    hmmbuild(blocks)
    // in: (clade, [blocks]); out: (clade, <clade>.fasta) = the supermatrix
    horizontal(filter_gappy_leaves.out.alns)

    // placement HMMs and the HMMs of every SCMG candidate
    // in: (clade, all_GA.hmm + 4 indices, profiles_<clade>.csv)
    hmmfetch(calibrated.join(choose_tree_profiles.out.profiles))
    // in: (clade, all_GA.hmm + 4 indices, scmgs_<clade>.gz)
    hmmfetch_db(calibrated.join(get_SCMG.out))

    // tree: the model comes from ModelFinder on one sub-database unless it is off
    if (params.iqtree) {
        // in: (clade, supermatrix) of params.modelfinder_clade only; out: model.txt
        model = params.modelfinder
            ? modelfinder(horizontal.out.filter { clade, _supermatrix -> clade == params.modelfinder_clade }, params.seed)
                .model.map { model_file -> model_file.text.trim() }.first()
            : channel.value(params.iqmodel)
        // in: (clade, supermatrix, model); out: (clade, tree)
        tree = iqtree(horizontal.out.combine(model), params.seed).tree
    } else {
        model = channel.value(params.iqmodel)
        // in: (clade, supermatrix); out: (clade, tree)
        tree = fasttree(horizontal.out).tree
    }
    // in: (clade, tree); out: (clade, rooted tree)
    midpoint_root(tree)
    // in: (clade, rooted tree, taxinfo_<clade>.csv, ncbi_table, taxdb)
    itol_annotations(midpoint_root.out.join(reduce_taxinfo.out).combine(ncbi_table).combine(taxdb))

    // one row per clade, in the order the db process reads it:
    // (clade, SCMG table, supermatrix, rooted tree, SCMG HMMs, placement HMMs, [blocks],
    //  [block HMMs], taxinfo, [MetaEuk files], model)
    db_inputs = get_SCMG.out
        .join(horizontal.out)
        .join(midpoint_root.out)
        .join(hmmfetch_db.out)
        .join(hmmfetch.out)
        .join(filter_gappy_leaves.out.alns)
        .join(hmmbuild.out.groupTuple())
        .join(reduce_taxinfo.out)
        .combine(metaeuk_files.map { metaeuk_db_files -> [metaeuk_db_files] })
        .combine(model)
    db(db_inputs, params.panther_version)

    emit:
    db = db.out
    report_files = genome_group.out.info.map { _clade, classification_csv -> classification_csv }
        .mix(optimize_bitscore.out.map { _clade, cutoffs_csv -> cutoffs_csv })
        .mix(choose_tree_profiles.out.profiles.map { _clade, profiles_csv -> profiles_csv })
        .mix(choose_tree_profiles.out.log)
        .mix(get_SCMG.out.map { _clade, scmg_table -> scmg_table })
        .mix(filter_gappy_leaves.out.report.map { _clade, leaf_gaps_tsv -> leaf_gaps_tsv })
        .mix(horizontal.out.map { _clade, supermatrix -> supermatrix })
    searched = hmmsearch_GA.out.map { clade, accession, _proteome, _ga_table -> "${clade}\t${accession}" }
}
