def encode_dayhoff(seq):
    '''
    Turn a protein sequence into its corresponding Dayhoff encoding. Return
    None if the encoding is not unique (for details on this see comment in
    scripts/encoding.py)

    https://en.wikipedia.org/wiki/Margaret_Oakley_Dayhoff

    Usage:

    seq = 'MAFSAEDVLKEYDRRRRMEALLLSLYYPNDRKLLDYKEWSPPRVQVEC'
    encode_dayhoff(seq)
    # 'ebfbbcceedcfcddddecbeeebeffbccddeecfdcfbbbdececa'
    '''
    dayhoff = {
        'C' + 'U': 'a', 
        'GSTAP': 'b',
        'DENQ' + 'Z': 'c',
        'RHK' + 'O': 'd',
        'LVMI' + 'J': 'e',
        'YFW': 'f'}

    encoding = ''

    for letter in seq:
        found = 0
        for k in dayhoff:    
            if letter in k:
                # Will only return one letter per loop bc/ keys unique
                encoding += dayhoff[k]
                found = 1
        if not found:
            return None
    assert len(encoding) == len(seq)
    return encoding
