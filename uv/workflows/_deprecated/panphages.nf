include { phage_containment } from './processes/search.nf'

include { 
    assemble
    phage_subgraphs
    phage_containment
    mec_tree
    } from './processes/panphages.nf'


workflow panphages {
    take:
    sketchlist
    candidates
    tree

    main:
    assemble(candidates)
    phage_subgraphs(assemble.out.graph)
    sketch_subgraphs(phage_subgraphs.out)
    phage_containment(sketchlist, sketch_subgraphs.out, tree)
    mec_tree(phage_containment.out)
    
    emit:
    phage_containment.out
}
