// Param "path_bin_<name>" has to be unique across the composite workflow!
// This block needs to come before include {...}
if (params.standalone) {
    params.path_bin_uv  = "${workflow.projectDir}/bin"
} else {
    params.path_bin_uv  = "${workflow.projectDir}/submodules/uv/bin"
}


include { 
    reading_frames
    search_protein_db
    search_protein_db_careful
    search_tails
    segment
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
    segment(search_protein_db.out.join(rename_contigs.out.names))

    // Assess the quality of these intervals and filter accordingly
    qc(segment.out.sequences.combine(checkv_db))
    filter_qc(qc.out.join(rename_contigs.out.names))

    emit:
    sequences  = filter_qc.out.sequences
    mask_broad = segment.out.mask
    mask_tight = filter_qc.out.mask
    positive   = filter_qc.out.positive
    names      = rename_contigs.out.names
}


workflow annotate {
    /*
    Carefully annotate the sequences identified as being of phage origin.
    */
    take:
    fragments
    names

    main:
    hmms = channel.fromPath("${params.uvdb}/pvog.hmm", checkIfExists: true)
    hmm_groups = channel.fromPath("${params.uvdb}/pvog.groups.csv", checkIfExists: true)

    prefix = 'targetDB'
    fp = "${params.uvdb}/phage_proteins/targetDB*"
    db_proteins = channel.fromPath(fp, checkIfExists: true) | collect

    db_tails = channel.fromPath("${params.uvdb}/tails.faa")

    // Reannotate reading frames using phage-specific ORF caller ...
    frames = fragments | careful_frames
    // ... and see which typical protein families/ domains can be identified.
    
    // Apply virus model(s) to the candidate intervals and filter them
    domains = frames | combine(hmms) | search_hmms 
    domains | combine(hmm_groups) | interpret_hmms | groupTuple | join(names) | collect_hmms
    
    all_frames = careful_frames.out.map { it -> it[1] } | collect
    foo = careful_frames.out | groupTuple | view
    
    search_protein_db_careful(all_frames, db_proteins, prefix)
    foo | combine(db_tails) | search_tails

    emit:
    annotation = interpret_hmms.out
}

