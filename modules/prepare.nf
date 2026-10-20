// Reference data: NCBI taxonomy and the PANTHER library.

process prepare_taxdb {
    tag "${dump.name}"
    label 'process_low'
    container 'biocontainers/ete3:3.1.2'

    input:
    path(dump)                        // the NCBI taxonomy dump (taxdump.tar.gz or a taxdmp_YYYY-MM-01.zip)
    val(source)                       // where the dump came from (--taxdump), written to taxonomy.tsv

    output:
    path("taxa.sqlite"), emit: taxdb  // taxa.sqlite ete3
    path("taxonomy.tsv"), emit: info  // source, md5 and date of the dump

    script:
    """
    taxdump.py targz ${dump} dump.tar.gz
    python -c "from ete3 import NCBITaxa; NCBITaxa(dbfile='taxa.sqlite', taxdump_file='dump.tar.gz')" > build.log 2>&1
    rm dump.tar.gz
    {
        printf 'source\\t%s\\n' '${source}'
        printf 'md5\\t%s\\n' "\$(md5sum ${dump} | cut -d' ' -f1)"
        printf 'ncbi_date\\t%s\\n' "\$(taxdump.py date ${dump})"
        printf 'built\\t%s\\n' "\$(date -u +%Y-%m-%d)"
        printf 'ete3\\t%s\\n' "\$(python -c 'import ete3; print(ete3.__version__)')"
    } > taxonomy.tsv
    """

    stub:
    """
    touch taxa.sqlite
    printf 'source\\t${source}\\n' > taxonomy.tsv
    """
}

// Nuclear genetic code of every candidate genome (nodes.dmp, column 7): ciliates,
// some diplomonads and kinetoplastids do not use the standard code. Reported in coverage.tsv only
process genetic_codes {
    label 'process_single'
    container 'biocontainers/python:3.9--1'

    input:
    path(dump)                 // the NCBI taxonomy dump
    path(ncbi_table)           // group,taxid,accession of every candidate

    output:
    path("genetic_codes.tsv")  // accession, taxid, nuclear genetic code

    script:
    """
    taxdump.py extract ${dump} nodes.dmp
    awk -F'\\t[|]\\t' 'FNR == NR { split(\$0, a, ","); acc[a[2]] = acc[a[2]] " " a[3]; next }
                       (\$1 in acc) { n = split(acc[\$1], x, " "); for (i = 1; i <= n; i++) print x[i] "\\t" \$1 "\\t" \$7 }' \\
        ${ncbi_table} nodes.dmp | sort > genetic_codes.tsv
    rm nodes.dmp
    """

    stub:
    """
    awk -F, '{print \$3 "\\t" \$2 "\\t1"}' ${ncbi_table} > genetic_codes.tsv
    """
}

// PANTHER library from the public tarball
process fetch_panther {
    tag "${panther_version}"
    label 'process_single'
    container 'biocontainers/python:3.9--1'
    storeDir "${workflow.stubRun ? "${workflow.workDir}/stub_cache" : params.cachedir}/panther/${panther_version}"

    input:
    path(tarball)                      // the PANTHER library tarball (--panther_url)
    val(panther_version)               // PANTHER version (--panther_version)

    output:
    path("PANTHER${panther_version}")  // PANTHER<version>/ with books/ (no subfamilies)

    script:
    """
    mkdir x
    tar xzf ${tarball} -C x --exclude='*/SF*'
    mv "\$(dirname "\$(find x -type d -name books | head -1)")" PANTHER${panther_version}
    rm -rf x
    """

    stub:
    """
    mkdir -p PANTHER${panther_version}/books
    """
}
