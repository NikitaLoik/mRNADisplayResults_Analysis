peptide = raw_input('insert a peptide string ')

nnuLibTrnsl = {'UUU':'F','CUU':'L','AUU':'I','AUG':'M','GUU':'V','UCU':'S','CCU':'P','ACU':'T',
'GCU':'A','UAU':'Y','UAG':'stop','CAU':'H','AAU':'N','GAU':'D','UGU':'C','CGU':'R','AGU':'S','GGU':'G'}

DNA = {'A':'T','C':'G','G':'C','T':'A'}



nnuLibRevTrnsl = dict([[v,k] for k,v in nnuLibTrnsl.items()])

def complement(dna):
    DNA = {'A':'T','C':'G','G':'C','T':'A'}
    cdna = ''
    for i in range (len(dna)):
        a = DNA[dna[i]]
        cdna = cdna + a
    return cdna
    

i = 0
rna = ''

while i < len(peptide):
    a = nnuLibRevTrnsl[peptide[i]]
    rna = rna + a
    i += 1

dna = rna.replace('U','T')

cdna = complement(dna)

rcdna = cdna[::-1]
StartBit = 'AGACCCAGACCCAGACCCACA'
EndBit = 'GGATATATCTCCTACTTAAAG'

primer = StartBit + rcdna + EndBit

print primer