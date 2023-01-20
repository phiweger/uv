include { 
    create_hash_ix
    create_kmer_ix
    } from './processes/setup.nf'


workflow setup {
    take:
        genomes

    main:
        create_hash_ix(genomes)
        // create_kmer_ix(genomes)

    emit:
        hash_ix = create_hash_ix.out
        kmer_ix = create_kmer_ix.out
}