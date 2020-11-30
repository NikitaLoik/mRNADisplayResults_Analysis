DNA = {'A':'T','C':'G','G':'C','T':'A'}

DNAEntropy = {}

DNAEnthalpy = {}

DNAFreeEnergy = {'AA':1.00,'AT':-0.88,'TA':-0.58,'CA':-1.45,'GT':-1.44,'CT':-1.28,'GA':-1.30,'CG':-2.17,'GC':-2.24,'GG':-1.84}

print DNAFreeEnergy['AA']

dna = 'AGACCCAGACCCAGACCCACA'

i = 0
j = 2
AnnealingFreeEnergy = 0

while j <= len(dna):
    AnnealingFreeEnergy = AnnealingFreeEnergy + DNAFreeEnergy[dna[i:j]]
    i += 1
    j += 1
print AnnealingFreeEnergy