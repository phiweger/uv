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