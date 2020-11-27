include {
    phigaro
    phastaf
    } from './processes/benchmark.nf'


workflow benchmark {
    take:
    genomes

    main:
    phigaro(genomes.take(2))
    phastaf(genomes.take(2))

}