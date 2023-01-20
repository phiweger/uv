nextflow.enable.dsl = 2


include { uv; annotate } from './workflows/uv.nf'
// include { crispr } from './workflows/crispr.nf'
// include { setup } from './workflows/setup.nf'
// include { search } from './workflows/search.nf'


workflow {
    /*
    uv -- Finding prophages using brute force

    Name from the fact that UV light induces prophages in bacterial isolates,
    ie it makes them switch from a temperate to a lysogen lifestyle whereby
    they "reveal themselves".

    Starting with a bacterial genome assembly, this workflow identifies
    regions that are likely (pro)phages.

    The workflow returns:

    - a (conservative) mask around regions of increased phage protein homology
    - a list of identified phages
    - a list of protein domains present across phage regions
    */

    genomes = channel.fromPath(params.genomes)
                     .splitCsv(header: true)
                     .map{ row -> tuple(row.name, row.path) }

    // phages = channel.fromPath(params.phages)
    // setup(phages)

    //sigs = channel.fromPath(params.signatures)

    // Find prophages
    uv(genomes)
    // prophages = uv.out.sequences.transpose().combine(setup.out.hash_ix)
    //prophages = uv.out.sequences.transpose().combine(sigs)

    // Annotate them carefully
    if (params.annotate) {
        annotate(uv.out.sequences.transpose(), uv.out.names)
    }

    //crispr(genomes).view()
    // phage_fragments = crispr.out.combine(setup.out.kmer_ix)

    // search(prophages, crispr.out)
    // search(prophages)
}








