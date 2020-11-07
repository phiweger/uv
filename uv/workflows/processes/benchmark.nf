process phigaro {
    container 'multifractal/phigaro:0.5.2'

    input:
    tuple(val(name), path(genome))

    output:
    tuple(val(name), path('phage.bed'))

    """
    phigaro --no-cleanup --delete-shorts -e html tsv bed -f ${genome} -o phage
    """
}


process phastaf {
    input:
    tuple(val(name), path(genome))

    output:
    tuple(val(name), path('out/phage.bed'))

    """
    phastaf --outdir out ${genome}
    """
}
