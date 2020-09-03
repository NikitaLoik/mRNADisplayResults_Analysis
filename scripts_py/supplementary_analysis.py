def au_count(dna_sequence):
    au_count = (
        dna_sequence.count('A')
        + dna_sequence.count('U')
        + dna_sequence.count('T'))

    return au_count