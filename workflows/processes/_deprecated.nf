process dereplicate {
    /*
    git clone https://github.com/rrwick/Assembly-dereplicator
    cp Assembly-dereplicator/dereplicator.py /usr/local/bin

    TODO: copy to workflow/bin?
    */
    input:
    path(fragments)

    output:
    path('dereplicated/*.fasta')

    """
    mkdir fragments
    # To avoid argument list too long issue use for loop:
    # stackoverflow, 11289551
    for i in *.fasta; do mv \$i fragments; done
    TH=\$(awk 'BEGIN {print 1 - ${params.dereplicate}}')
    dereplicator.py --threshold \$TH fragments dereplicated
    """
}


process search_phage_sketches {
    // echo true

    input:
    path(candidate_phages)
    path(ix)

    output:
    path('phage_search.csv', emit: found)
    path('phage_search.sig', emit: sketches)

    """
    cat ${candidate_phages} > all
    search_phage_sketches.py \
        --candidates all \
        -k ${params.k_finegrained} \
        --scaled ${params.rate_finegrained} \
        --min-containment ${params.min_containment} \
        --index ${ix} \
        --prefix phage_search \
        --minlen ${params.min_phage_len}
    """
}


process checksum {
    input:
    path(genome)

    output:
    tuple(env(checksum), path(genome))

    shell:
    '''
    checksum=$(md5sum !{genome} | awk '{ printf $1 }')
    '''
}


process sketch {
    input:
    tuple(val(name), path(genome))

    output:
    tuple(val(name), path("${genome}.sig"))

    """
    sourmash compute -k 21 --scaled 100 ${genome}
    """
    // split in a sepatate process to paralellize or split channel btw/
    // fasta and bed
    // also, create gzipped SBT index
}


process search_sketch {
    /*
    Search a sketch database using queries. This process is very inefficient
    bc/ it will load the database every time! This is a known problem.

    https://github.com/dib-lab/sourmash/issues/475

    For now we bow our head in acceptance and hope for the future.
    */
    // echo true
    memory '8 GB'  // limits number of instances run in parallel
    
    input:
    tuple(val(name), path(sketches), path(db))

    output:
    tuple(val(name), path('found.txt')) optional true

    shell:
    '''
    # sourmash signature cat -o all ${sketches}
    
    for i in !{sketches}; do
        sourmash search --containment \
            -o ${i}.csv --threshold !{params.min_containment} \
            $i !{db}

        if [[ $(wc -l < ${i}.csv) -ge 2 ]]; then
            head -n2 ${i}.csv | tail -n1 | cut -f1 | cut -d, -f2 > ${i}.txt
        fi
    done

    # Did we find anything?
    count=$(ls -1 *.txt 2> /dev/null | wc -l)
    if [ $count != 0 ]; then
        cat *.txt > found.txt
    fi 
    '''
    // Format of the sourmash search result:
    // 0.549...,uvig_394964   SRR3132147_71 length_58589_VirSorter_...
    // --best-only does not work w/ --containment? Bug?
}


process collect_phage_genome {
    publishDir "${params.results}/similar_phages"

    input:
    tuple(val(name), val(found), path(genomes), path(index))

    output:
    tuple(val(name), path("${found}.fasta"))

    """
    samtools faidx ${genomes} ${found} > ${found}.fasta
    """
}


process exact_align {
    /*
    - TODO: https://github.com/sanger-pathogens/pymummer
    - TODO: mask the exact alignment at > 99%
    */
    input:
    tuple(val(name), path(target), path(query))

    output:
    path('out.coords')

    """
    if [[ ${target} =~ \\.gz\$ ]]
        then
            gunzip -c ${target} > unpacked
        else
            cat ${target} > unpacked
    fi

    nucmer unpacked ${query} > out.delta
    show-coords out.delta > out.coords
    """
}


process index {
    input:
    path(genomes)

    output:
    tuple(path(genomes), path("${genomes}.fai"))

    """
    samtools faidx ${genomes}
    """
}


process sketch_finegrained {
    input:
    tuple(val(name), path(genome))

    output:
    tuple(val(name), path("${name}.sig"))

    """
    sourmash compute -k ${params.k_finegrained} --scaled 100 --name ${name} -o ${name}.sig ${genome}
    """
}

