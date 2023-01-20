include { 
    search_kmer
    search_hash
    } from './processes/search.nf'


workflow search {
    take:
        prophages
        fragments

    main:
        graph = channel.fromPath(
            '/Volumes/shed/gut_phage_database/graph.dbg',
            checkIfExists: true).view()
        anno = channel.fromPath(
            '/Volumes/shed/gut_phage_database/graph.column.annodbg',
            checkIfExists: true).view()
        
        search_hash(prophages)
        search_kmer(fragments.combine(graph).combine(anno).view())
}
