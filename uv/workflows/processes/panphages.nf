
process assemble {
    // echo true
    container 'nanozoo/flye:2.8--95b6dca'
    // publishDir params.outdir, mode: 'copy', overwrite: true
    cpus 8

    input:
    path(fragments)

    output:
    path('assembly/assembly_graph.gfa', emit: graph)
    path('assembly/assembly.fasta', emit: contigs)

    """
    find . -name "*.fasta" -exec cat '{}' ';' > all.fasta
    
    # There are duplicate contig names, so we rename them (running count)
    awk '/^>/{print ">" ++i; next}{print}' < all.fasta > all.renamed.fasta
    
    flye --min-overlap ${params.min_assembly_overlap} --subassemblies all.renamed.fasta --threads 8 --out-dir assembly
    # --iterations 0
    """
}


process phage_subgraphs {
    input:
    path(assembly)

    output:
    path('subgraphs/*.fasta')

    """
    mkdir subgraphs
    decompose_assembly.py --graph ${assembly} --outdir subgraphs --minlen ${params.min_phage_len}
    """
}



process sketch_subgraphs {
    /*
    Create a single signature file containing n sketches, one for each subgraph
    */
    input:
    path(sketches)

    output:
    path('panphage.sig')

    """
    sourmash compute -k ${params.k_finegrained} --scaled ${params.rate_finegrained} -o panphage.sig ${sketches}
    """
}


process mec_tree {
    /*
    Tree created from distance of scaled MinHash containment distances for
    mobile elements, namely plasmids and phages.
    */
    // echo true

    input:
    path(containment)

    output:
    path('mec_tree.nwk')

    """
    containment_tree.py -i ${containment} -o mec_tree.nwk
    """
}