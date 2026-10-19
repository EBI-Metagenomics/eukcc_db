// Step 3 of docs/workflow.md (3a-3b): MetaEuk reference and protein prediction.

include { metaeuk_db; metaeuk_pred } from '../modules/metaeuk'

workflow proteomes {
    take:
    group_genomes      // (group, [genomes])
    metaeuk_proteomes  // annotated RefSeq proteomes of the MetaEuk set
    genetic_codes      // genetic_codes.tsv: accession, taxid, nuclear code

    main:
    // 3a. in: every MetaEuk set proteome in one file; out: [microeuks_metaeuk.db*]
    metaeuk_files = metaeuk_db(metaeuk_proteomes.collectFile(name: 'metaeuk_proteins.faa', sort: { proteome -> proteome.name })).db.collect()

    // accession -> genetic code
    code_of_accession = genetic_codes.map { codes_tsv ->
        codes_tsv.readLines().collectEntries { line ->
            def columns = line.tokenize('\t')
            [(columns[0]): columns[2]]
        }
    }
    // (group, genome, genetic code); a code MetaEuk does not support falls back to 1
    genomes_to_predict = group_genomes
        .transpose()
        .combine(code_of_accession)
        .map { group, genome, codes ->
            def code = codes[genome.baseName] ?: '1'
            def metaeuk_codes = ((1..6) + (9..16) + (21..31))*.toString()
            tuple(group, genome, code in metaeuk_codes ? code : '1')
        }
    // 3b. in: (group, genome, genetic code, [MetaEuk files]); out: (group, <accession>.faa)
    metaeuk_pred(genomes_to_predict.combine(metaeuk_files.map { metaeuk_db_files -> [metaeuk_db_files] }))

    emit:
    proteomes = metaeuk_pred.out     // (group, <accession>.faa)
    metaeuk_files = metaeuk_files    // [microeuks_metaeuk.db*]
}
