#-------------------------------------------------------------------------------
# utilises ONLY ONE stop codons, returns ONLY ONE ORF
def OpenReadingFrame(DNASequence):
    StartCodon = 'ATG'
    StopCodon = 'TAG'
    MinLength = 12
    MaxLength = 120
    SubString = DNASequence
    
    while StartCodon in SubString:
        StartIndex = SubString.find(StartCodon)
        StartString = SubString[SubString.find(StartCodon):]
        while StopCodon in StartString:
            StopIndex = StartString.rfind(StopCodon)
            ORF = StartString[:StopIndex+3]
            if MinLength <= len(ORF) and len(ORF) <= MaxLength and len(ORF)%3 == 0:
                return str(ORF)
            StartString = StartString[:StopIndex+2]
        SubString = SubString[StartIndex+1:]
        
a = 'ATGCGTCATAGGCTTCCTGTTGATGTTTCTAGGTATATTTGTCGTTAGTGCGGCAGCGGCAGCGGCAGCTAG'
print OpenReadingFrame(a)