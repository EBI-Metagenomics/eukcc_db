// Step 1 of docs/workflow.md (1a-1d): the inputs, given by the user or fetched.

include { prepare_taxdb; fetch_panther }  from '../modules/prepare'
include { get_group_summary }             from '../modules/ncbi'
include { concat_panther }                from '../modules/hmmer'

workflow inputs {
    main:
    if (params.cutoff_date && !params.taxdump.toString().contains('taxdump_archive')) {
        def month = params.cutoff_date.toString().take(7)
        log.warn "--cutoff_date ${params.cutoff_date}: the ${params.taxdb ? 'genetic codes' : 'taxonomy and the genetic codes'} " +
                 "still come from --taxdump ${params.taxdump}. The NCBI dump of that month: " +
                 "--taxdump https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_${month}-01.zip"
    }

    // 1a. taxonomy: a ready taxa.sqlite (--taxdb), or one built from the NCBI dump (--taxdump)
    if (params.taxdb) {
        taxdb = channel.value(file(params.taxdb, checkIfExists: true))
        taxonomy_info = channel.empty()
    } else {
        // in: the taxdump file (tar.gz or zip), its source (--taxdump, written to taxonomy.tsv)
        prepare_taxdb(file(params.taxdump), params.taxdump)
        taxdb = prepare_taxdb.out.taxdb
        taxonomy_info = prepare_taxdb.out.info
    }

    // 1b. assembly_summary.txt per division: from --summaries or downloaded from --ncbi_url.
    divisions = channel.fromList(params.divisions) // [group, catalogue] pairs, e.g. ['fungi', 'refseq']
    summaries = params.summaries
        ? divisions.map { group, catalogue -> tuple(catalogue, group, file("${params.summaries}/${catalogue}_${group}.txt", checkIfExists: true)) }
        : get_group_summary(divisions, params.ncbi_url)

    // 1c. PANTHER: a local directory holding books/ (--panther), or the public tarball
    // (--panther_url) of --panther_version
    panther = params.panther ? channel.value(file(params.panther, type: 'dir', checkIfExists: true))
                             : fetch_panther(file(params.panther_url), params.panther_version)

    // 1d. in: PANTHER directory; out: (all.hmm + 4 indices)
    panther_db = concat_panther(panther)

    emit:
    taxdb = taxdb                    // taxa.sqlite
    taxonomy_info = taxonomy_info    // taxonomy.tsv: which dump the taxa.sqlite comes from
    summaries = summaries            // (catalogue, group, <catalogue>_<group>.txt)
    panther = panther                // PANTHER directory (books/)
    panther_db = panther_db          // (all.hmm, .h3f, .h3i, .h3m, .h3p)
}
