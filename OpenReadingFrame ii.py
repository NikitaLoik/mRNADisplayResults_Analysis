# utilises different stop codons, finds all possible Open Reading Frames
#-------------------------------------------------------------------------------
def OpenReadingFrame(DNASequence):
    StartCodon = 'ATG'
    StopCodonList = ['TAA','TAG','TGA']
    MinLength = 12
    MaxLength = 120
    ORFList = []
    i = 0
    
    for StopCodon in StopCodonList:
        SubString = DNASequence
        while StartCodon in SubString:
            StartIndex = SubString.find(StartCodon)
            StartString = SubString[SubString.find(StartCodon):]
        
            while StopCodon in StartString:
                StopIndex = StartString.rfind(StopCodon)
                i += 1
                ORF = StartString[:StopIndex+3]
                if MinLength <= len(ORF) and len(ORF) <= MaxLength and len(ORF)%3 == 0:
                    ORFList.append(ORF)
                StartString = StartString[:StopIndex+2]
            SubString = SubString[StartIndex+1:]
            
    FinalisedORFList = '\n'.join(ORFList)
    return FinalisedORFList      
#-------------------------------------------------------------------------------        
a = 'ATGCGTCATAGGCTTCCTGTTGATGTTTCTAGGTATATTTGTCGTTAGTGCGGCAGCGGCAGCGGCAGCTAG'
print OpenReadingFrame(a)