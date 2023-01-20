process find_crispr {
    /*
    Check if file empty -- stackoverflow.com/questions/9964823
    */
    // container 'nanozoo/prokka:1.13.4--d6a71cb'
    // publishDir "${params.results}/${name}", mode: 'copy', overwrite: true

    input:
        tuple(val(name), path(genome))

    output:
        tuple(val(name), path('crispr_spacers.fa')) optional true

    """
    minced -spacers ${genome} crispr

    [[ -s crispr_spacers.fa ]] || rm crispr_spacers.fa
    """
}


workflow crispr {
    take:
        genome

    main:
        find_crispr(genome)

    emit:
        find_crispr.out

}
