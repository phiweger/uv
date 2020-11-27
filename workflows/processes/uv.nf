process reading_frames {
    input:
    tuple(val(name), path(genome))

    output:
    tuple(val(name), path(genome), path("${name}.reading_frames.faa"))

    """
    if [[ ${genome} =~ \\.gz\$ ]]
        then
            gunzip -c ${genome} > unpacked
        else
            cat ${genome} > unpacked
    fi

    prodigal -i unpacked -a ${name}.reading_frames.faa > /dev/null
    """
}


process search_protein_db {
    memory "${params.maxram} GB"  
    // Limits number of instances run in parallel

    input:
    tuple(val(name), path(genome), path(proteins))
    path(db)
    val(prefix)

    output:
    tuple(val(name), path(genome), path(proteins), path('aln.m8'))

    """
    mmseqs easy-search --min-seq-id 0.95 -c 0.95 ${proteins} ${prefix} aln.m8 tmp
    """
}


process intervals {
    publishDir "${params.results}/${name}", mode: 'copy', overwrite: true, pattern: '*.bed'
    input:
    tuple(val(name), path(genome), path(proteins), path(hits), path(names))

    output:
    tuple(val(name), path("results/*.fasta"), emit: sequences)
    tuple(val(name), path("putative.bed"), emit: mask)

    """
    ${params.path_bin_uv}/find_phage_breakpoints.py --genome ${genome} --frames ${proteins} --hits ${hits} --threshold 1 --names ${names} --outdir results

    cat results/*.bed > tmp
    bedtools sort -i tmp > putative.bed
    """
}


process careful_frames {
    /*
    Phanotate will translate frames with eg GTG as the first codon -- this will
    be translated into methionin when in first position, as are other codons
    other then ATG -- https://www.biostars.org/p/364080/

    The translation to aa does not bother correcting this, bc/ it should not
    make a huge difference to subsequent HMM domain search.
    */
    errorStrategy 'ignore'  // if eg N in sequence

    input:
    tuple(val(name), path(genome))

    output:
    tuple(val(name), path('frames.faa'))

    """
    phanotate.py --outfmt fasta ${genome} > frames.fna
    map_into_global_coordinates.py -i frames.fna -o frames.recoord.fna

    ${params.path_bin_uv}/translate.py -i frames.recoord.fna -o frames.faa
    """
    // pip3 install phanotate
    // conda install -y -c bioconda trnascan-se
    // phanotate.py -o test.fasta -f fasta maybe_phage.fasta
    // 792 ORFs vs 564 from Prodigal
}


process qc {
    /*
    > Note: CheckV will not detect proviruses if host regions have already been removed (e.g. using VIBRANT or VirSorter) -- https://bitbucket.org/berkeleylab/checkv/src/master/

    In these cases, genome completeness is our measure of quality.
    */
    publishDir "${params.results}/${name}"
    memory '8 GB'  // limits number of instances run in parallel
    // container 'nanozoo/checkv:0.7.0--097a445'

    input:
    tuple(val(name), path(genomes), path(db))

    output:
    tuple(val(name), path('qc.tar.gz'))

    // optional true

    """
    cat ${genomes} > all
    checkv end_to_end -t 8 -d ${db} all qc 2>&1 > /dev/null
    
    tar -zcvf qc.tar.gz qc 
    """
}


process filter_qc {
    publishDir "${params.results}/${name}", mode: 'copy', overwrite: true
    // echo true

    input:
    tuple(val(name), path(qc_results), path(contig_map))
    
    output:
    tuple(val(name), path('sequences/*.fasta'), emit: sequences) \
        optional true
    tuple(val(name), path('phages.bed'), emit: mask) \
        optional true
    tuple(val(name), path('contigs.txt'), emit: positive) \
        optional true

    """
    tar -zxvf ${qc_results}
    
    # if no .fna, don't throw error
    # stackoverflow, 12829081
    cat qc/*.fna > all 2> /dev/null
    
    if [ -s all ]; then 
        ${params.path_bin_uv}/filter_checkv.py \
            --id ${name} \
            --qc-results qc \
            --min-viral-genes ${params.min_viral_genes} \
            --sequences all \
            --contig-map ${contig_map} \
            --minlen ${params.min_phage_len} \
            --outdir sequences --force

        mv sequences/phages.bed phages.bed
        cut -f1 phages.bed | sort | uniq > contigs.txt
    fi
    """
}


process search_hmms {
    /*
    On thresholds:
    https://hmmer-web-docs.readthedocs.io/en/latest/searches.html#thresholds
    */
    input:
    tuple(val(name), path(proteins), path(hmms))

    output:
    tuple(val(name), path('proteins.txt'))

    """
    hmmsearch -E ${params.min_e_hmm} --noali --cpu 8 --tblout proteins.tblout ${hmms} ${proteins} > /dev/null
    grep -v '^#' proteins.tblout | sed 's/ * / /g' | cut -f1,3,8 -d ' ' --output-delimiter '\t' > proteins.txt
    """
    // Parse hmsearch result like a boss:
    // https://madsalbertsen.github.io/multi-metagenome/docs/step5.html
}


process interpret_hmms {
    input:
    tuple(val(name), path(hits), path(groups))

    output:
    tuple(val(name), path("*.bed"))

    shell:
    '''
    !{params.path_bin_uv}/lookup_hmms.py --hits !{hits} --groups !{groups} -o tmp
    # we need a unique name
    checksum=$(md5sum tmp | awk '{ printf $1 }')
    mv tmp ${checksum}.bed
    '''
}


process collect_hmms {
    publishDir "${params.results}/${name}", mode: 'copy', overwrite: true

    input:
    tuple(val(name), path(beds), path(names))

    output:
    tuple(val(name), path('annotation.bed'))

    """
    cat *.bed > all
    bedtools sort -i all > sorted
    ${params.path_bin_uv}/deduplicate_and_rename.py -i sorted -o annotation.bed --names ${names}
    """
}


process minlen {
    input:
    tuple(val(name), path(genome))

    output:
    tuple(val(name), path('reference.fasta'))

    script:
    """
    ${params.path_bin_uv}/minlen.py -i ${genome} -o reference.fasta --minlen ${params.min_contig_len}
    """
}


process rename_contigs {
    /*
    Will unzip, too.
    */
    input:
    tuple(val(name), path(genome))

    output:
    tuple(val(name), path("${name}.renamed_contigs.fasta"), emit: genome)
    tuple(val(name), path("${name}.contig_names.txt"), emit: names)

    """
    ${params.path_bin_uv}/rename.py -i ${genome} \
        --sequences ${name}.renamed_contigs.fasta \
        --names ${name}.contig_names.txt
    """
}

