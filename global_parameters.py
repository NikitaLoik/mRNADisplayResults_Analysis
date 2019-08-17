QUALITY_SCORE_STRING = '''
    !"#$%&'()*+,-./0123456789:;<=>?@
    ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`
    abcdefghijklmnopqrstuvwxyz{|}~
    '''

QUALITY_SCORE = 29  # threshold quality_score must be between 0 and 93

TRANSLATION_CODE = {
    'AAA':'K','AAC':'N','AAG':'K','AAU':'N',
    'ACA':'T','ACC':'T','ACG':'T','ACU':'T',
    'AGA':'R','AGC':'S','AGG':'R','AGU':'S',
    'AUA':'I','AUC':'I','AUG':'M','AUU':'I',
    
    'CAA':'Q','CAC':'H','CAG':'Q','CAU':'H',
    'CCA':'P','CCC':'P','CCG':'P','CCU':'P',
    'CGA':'R','CGC':'R','CGG':'R','CGU':'R',
    'CUA':'L','CUC':'L','CUG':'L','CUU':'L',
    
    'GAA':'E','GAC':'D','GAG':'E','GAU':'D',
    'GCA':'A','GCC':'A','GCG':'A','GCU':'A',
    'GGA':'G','GGC':'G','GGG':'G','GGU':'G',
    'GUA':'V','GUC':'V','GUG':'V','GUU':'V',
    
    'UAA':'#','UAC':'Y','UAG':'*','UAU':'Y',
    'UCA':'S','UCC':'S','UCG':'S','UCU':'S',
    'UGA':'&','UGC':'C','UGG':'W','UGU':'C',
    'UUA':'L','UUC':'F','UUG':'L','UUU':'F'
    }


TRANSCRIPTION_CODE = {
    'A':'A','C':'C','G':'G','T':'U','U':'T'
    }

CDNA_MIN_LENGTH = 24
CDNA_MAX_LENGTH = 240

START_SEQUENCE = 'TAATACGACTCACTATAGGGTTAACTTTAAGAAGGAGATATACATATG'    # NNK - T7g10M.F48 
STOP_SEQUENCE = 'TGCGGCAGCGGCAGCGGCAGCTAGGACGGGGGGCGGAAA' #NNK - CGS3an13.R39
# START_SEQUENCE = 'TAATACGACTCACTATAGGGTTGAACTTTAAGTAGGAGATATATCCATG'   #NNU - T7-CH-F49
# STOP_SEQUENCE = 'TGTGGGTCTGGGTCTGGGTCTTAGGACGGGGGGCGGAAA'  #NNU - CGS3-CH-R39
# START_SEQUENCE = 'ATG' # Met codon
# STOP_SEQUENCE = 'TGCGGCAGC'# Akane seams to have trimmed siquences
# START_SEQUENCE = 'TAGGGTTAACTTTAAGAAGGAGATATACATATG'# Oxford, Akane and Tom
# STOP_SEQUENCE = 'TGCGGC'# Oxford, Akane and Tom
# STOP_SEQUENCE = 'TAG' # amber stop codon