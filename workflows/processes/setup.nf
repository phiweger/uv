process create_hash_ix {
    container 'nanozoo/sourmash:3.5.0--bd13d14'
    cpus 8

    input:
        path(genomes)

    output:
        path('genomes.sbt.zip')

    """
    sourmash compute -p ${task.cpus} --singleton -k ${params.ksize} --scaled ${params.scaled} -o sigs ${genomes}
    sourmash index genomes.sbt.zip sigs
    """
}


process create_kmer_ix {
    // container 'nanozoo/metagraph_dna:0.1--184de8f'
    // docker run --cpus=8 --memory=12g -v $PWD:/data ratschlab/metagraph:latest build -k 21 -v -p 8 --mem-cap-gb 10 -o /data/graph /data/GPD_sequences.fa.gz
    container 'ratschlab/metagraph:latest'
    cpus 8
    memory '12 GB'

    input:
        path(genomes)

    output:
        tuple(path('graph.dbg'), path('graph.column.annodbg'))

    """
    metagraph build -k 21 -v -p ${task.cpus} --mem-cap-gb ${params.maxram} -o graph ${genomes}
    metagraph annotate -v -p ${task.cpus} -i graph.dbg --anno-header -o graph ${genomes}
    """
    // k=21
    // > CRISPR repeats typically range in size from 28 to 37 base pairs (bps), though there can be as few as 23 bp and as many as 55 bp. [...] The size of spacers in different CRISPR arrays is typically 32 to 38 bp (range 21 to 72 bp) -- https://en.wikipedia.org/wiki/CRISPR
}

