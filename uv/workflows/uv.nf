include { 
    reading_frames
    search_protein_db
    intervals
    careful_frames
    qc
    filter_qc
    search_hmms
    interpret_hmms
    collect_hmms
    rename_contigs
    minlen
    } from './processes/uv.nf'


workflow uv {
    take:
    genome

    main:
    // Housekeeping
    prefix = 'targetDB'
    db_proteins = channel.fromPath(
        "${params.uvdb}/phage_proteins/targetDB*",
        checkIfExists: true)
    
    checkv_db = channel.fromPath(
        "${params.uvdb}/checkv",
        checkIfExists: true)

    // Action
    rename_contigs(minlen(genome))
    reading_frames(rename_contigs.out.genome)

    // Do these reading frames match known phage proteins?
    search_protein_db(reading_frames.out, db_proteins.collect(), prefix)
    // Extract the intervals around sequences of hits
    intervals(search_protein_db.out.join(rename_contigs.out.names))
    
    // Assess the quality of these intervals and filter accordingly
    qc(intervals.out.sequences.combine(checkv_db))
    filter_qc(qc.out.join(rename_contigs.out.names))

    emit:
    sequences        = filter_qc.out.sequences
    mask_broad       = intervals.out.mask
    mask_tight       = filter_qc.out.mask
    positive         = filter_qc.out.positive
    names            = rename_contigs.out.names
}


workflow annotate {
    /*
    Carefully annotate the sequences identified as being of phage origin.
    */
    take:
    fragment
    names

    main:
    hmms = channel.fromPath(
        "${params.uvdb}/pvog.hmm",
        checkIfExists: true)
    hmm_groups = channel.fromPath(
        "${params.uvdb}/pvog.groups.csv",
        checkIfExists: true)

    // Reannotate reading frames using phage-specific ORF caller ...
    careful_frames(fragment)
    // ... and see which typical protein families/ domains can be identified.
    
    // Apply virus model(s) to the candidate intervals and filter them
    search_hmms(careful_frames.out.combine(hmms))
    interpret_hmms(search_hmms.out.combine(hmm_groups))
    collect_hmms(interpret_hmms.out.groupTuple().join(names))

    emit:
    annotation = interpret_hmms.out
}

