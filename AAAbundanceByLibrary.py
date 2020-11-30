def CodonAndAAAbundanceByLybrary(LibrariesList, LibratiesRatioList):
    
ComplementaryDNA = {'A':'T','C':'G','G':'C','T':'A'}

ComplementaryRNA = {'A':'U','C':'G','G':'C','U':'A'}

ExtendedDNA = {'A':'A','C':'C','G':'G','T':'T',
                'R':('A','G'),'Y':('C','T'),'S':('G','C'),'W':('A','T'),'K':('G','T'),'M':('A','C'),
                'B':('C','G','T'),'D':('A','G','T'),'H':('A','C','T'),'V':('A','C','G'),
                'N':('A','C','G','T')}
                
transcription = {'A':'A','C':'C','G':'G','T':'U','U':'U'}

translation = {'AAA':'K','AAC':'N','AAG':'K','AAU':'N',
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
                
                'UAA':'ochre','UAC':'Y','UAG':'amber','UAU':'Y',
                'UCA':'S','UCC':'S','UCG':'S','UCU':'S',
                'UGA':'opal','UGC':'C','UGG':'W','UGU':'C',
                'UUA':'L','UUC':'F','UUG':'L','UUU':'F'}

# (1.A) expand the libraries

# (1.B) transcribe libraries into RNA sequences

# (2.A) convert the libraries into lists of sequences
