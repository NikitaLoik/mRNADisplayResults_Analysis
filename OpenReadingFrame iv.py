#-------------------------------------------------------------------------------
# utilises ONLY ONE stop codons, returns ONLY ONE ORF
def OpenReadingFrame(DNASequence):
    StartCodon = 'ATG'
    StopCodon = 'TAG'
    MinLength = 12
    MaxLength = 120
    SubString = DNASequence
    ORF = ''
    RF = ''
    if StartCodon in SubString:
        ORF = ORF + StartCodon     
        RFStartIndex = SubString.find(StartCodon) + 3
        RFStopIndex = SubString.find(StartCodon) + 3
        while RF != StopCodon:
            RF = SubString[RFStartIndex:RFStopIndex]
            ORF = ORF + RF
            RFStartIndex += 3
            RFStopIndex += 3
            
    return ORF
        
        
        

        
a = 'ATGAAATAG'
print OpenReadingFrame(a)