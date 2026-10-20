// Step 3 of docs/workflow.md (3a-3b): MetaEuk reference and protein prediction.

include { metaeuk_db; metaeuk_pred } from '../modules/metaeuk'

workflow proteomes {
    take:
    group_genomes      // (group, [genomes])
    metaeuk_proteomes  // annotated RefSeq proteomes of the MetaEuk set

    main:
    // 3a. in: every MetaEuk set proteome in one file; out: [microeuks_metaeuk.db*]
    metaeuk_files = metaeuk_db(metaeuk_proteomes.collectFile(name: 'metaeuk_proteins.faa', sort: { proteome -> proteome.name })).db.collect()

    // 3b. in: (group, genome, [MetaEuk files]); out: (group, <accession>.faa)
    metaeuk_pred(group_genomes.transpose().combine(metaeuk_files.map { metaeuk_db_files -> [metaeuk_db_files] }))

    emit:
    proteomes = metaeuk_pred.out     // (group, <accession>.faa)
    metaeuk_files = metaeuk_files    // [microeuks_metaeuk.db*]
}
