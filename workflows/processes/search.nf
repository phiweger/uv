process phage_containment {
    // echo true
    publishDir "${params.results}"

    input:
    path(query_sketches)
    path(phage_sketches)
    path(tree)

    output:
    path('phage_search.containment.csv')

    """
    ${workflow.projectDir}/bin/phage_containment.py \
        --query-sketches ${query_sketches} \
        --phage-sketches ${phage_sketches} \
        --tree ${tree} \
        -k ${params.k_finegrained} \
        --scaled ${params.rate_finegrained} \
        --deduplicate ${params.deduplicate} \
        --prefix phage_search
    """
}


process search_hash {
    publishDir "${params.results}/${name}", mode: 'copy', overwrite: true
    container 'nanozoo/sourmash:3.5.0--bd13d14'
    cpus 8

    input:
        tuple(val(name), path(prophages), path(ix))
    
    output:
        tuple(val(name), path('db_hits_prophage.csv')) optional true

    """
    for i in ${prophages}; do
        sourmash compute -p ${task.cpus} -k ${params.ksize} --scaled ${params.scaled} -o \${i}.sig \$i
        sourmash search --containment -o \${i}.hits.csv --threshold ${params.threshold} \${i}.sig ${ix}
    done

    # remove csv header similarity,name,filename,md5
    cat *.hits.csv | sed '/.*md5.*/d' | cut -f1-2 -d, > db_hits_prophage.csv
    [[ -s found.csv ]] || rm db_hits_prophage.csv
    """
}


process search_kmer {
    publishDir "${params.results}/${name}", mode: 'copy', overwrite: true
    container 'nanozoo/metagraph_dna:0.1--184de8f'
    cpus 8

    input:
        tuple(val(name), path(crispr), path(graph), path(anno))

    output:
        tuple(val(name), path('db_hits_crispr.csv')) optional true

    shell:
    '''
    metagraph query --num-top-labels 1 --suppress-unlabeled \
        -p !{task.cpus} \
        --discovery-fraction 0.9 \
        -i !{graph} \
        -a !{anno} \
        !{crispr} \
    > tmp

    cut -f3 tmp | sed "s/:/\\n/g" | sort | uniq > db_hits_crispr.csv
    # uniq -c
    '''
}
