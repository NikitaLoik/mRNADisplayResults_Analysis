#-------------------------------------------------------------------------------
# utilise ONLY ONE stop codons, returns ONLY ONE ORF
def OpenReadingFrame(DNASequence):
    StartCodon = 'ATG'
    StopCodon = 'TAG'
    MinLength = 42
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
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# translate DNA sequence
def Translation(DNASequence):
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
                    'UUA':'L','UUC':'F','UUG':'L','UUU':'F'
                    }
    transcription = {'A':'A','C':'C','G':'G','T':'U','U':'T'}
    
    RNASequence = ''
    for character in DNASequence:
        RNASequence = RNASequence + transcription.get(character,'X')
#    print RNASequence
        
    PeptideSequence = ''
    
    while len(RNASequence) != 0:
        PeptideSequence = PeptideSequence + translation.get(RNASequence[0:3],"Don't fuck with me!")
        RNASequence = RNASequence[3:]
    return PeptideSequence
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# read data file (.fastq extension)
#RawDataFile = open("/Users/NikitaLoik/Documents/R6.fastq", 'r')
#Lines = RawDataFile.readlines()
#RawDataFile.close
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# creat a .CSV file with sequence list, DNA sequences and peptide sequences
#DNASequenceFile = open('TestR6.csv', 'w')
#SelectionRound = 6
#SequenceCounter = 0
#for Line in Lines:
#    ORF = OpenReadingFrame(Line)
#    if ORF != None:
#        PeptideSequence = Translation(ORF)
#        SequenceCounter += 1
#        DNASequenceFile.write('selection round # ' + str(SelectionRound) + ',' +
#                            'ORF ' + str(SequenceCounter) + ',' +
#                            ORF + ',' +
#                            PeptideSequence + '\n')
#DNASequenceFile.close 
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# return a list of lists with peptide-sequences and their frequencies, sorted by frequency in descending order
def SortedPeptideSequencesList(fastqFileLocation):
    RawDataFile = open(fastqFileLocation, 'r')
    Lines = RawDataFile.readlines()
    RawDataFile.close
    PeptideSequences = {}
    PeptideSequencesList = []
    
    SelectionRoundNumber = fastqFileLocation[fastqFileLocation.find('.')-1]
    
    # populate the dictionary, so that Peptides are the keys and 
    for Line in Lines:
        ORF = OpenReadingFrame(Line)
        if ORF != None:
            Peptide = Translation(ORF)
            if Peptide not in PeptideSequences:
                PeptideSequences[str(Peptide)] = 1
            else:
                PeptideSequences[str(Peptide)] = PeptideSequences[str(Peptide)] + 1

    # convert the dictionary into the list of lists
    for key, value in PeptideSequences.iteritems():
        PeptideSequencesList.append([str(SelectionRoundNumber), key, value])
    # sort the PeptideSequenceList by peptide sequence occurence in descendent order
    SortedPeptideSequences = sorted(PeptideSequencesList, key = lambda x: x[2], reverse = True)
    return SortedPeptideSequences
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# creat a .CSV file with peptide-sequences list and their frequency, sorted by frequency in descending order
def SelectionRoundSortedSequenceListGenerator(SortedPeptideSequecesFileName, fastqFileLocation):
    
    SelectionRoundNumber = fastqFileLocation[fastqFileLocation.find('.')-1]
    
    SortedSequenceFile = open(SortedPeptideSequecesFileName, 'w')
    SortedPeptideSequences = SortedPeptideSequencesList(fastqFileLocation)
    
    TotalPeptideNumber = 0
    for Data in SortedPeptideSequences:
        TotalPeptideNumber = TotalPeptideNumber + Data[2]
        
    SortedSequenceFile.write('selection round # ' + str(SelectionRoundNumber) + '\n' +
                            'total sequence # ' + str(TotalPeptideNumber) + '\n')
                                            
    UniqueSequenceNumber = 0
    for Data in SortedPeptideSequences:
        UniqueSequenceNumber += 1
        PeptideSequence = Data[1]    
        SequenceFraction = float(Data[2])/float(TotalPeptideNumber)
        
        SortedSequenceFile.write('ORF ' + str(UniqueSequenceNumber) + ',' +
                                str(SelectionRoundNumber) + ',' +
                                PeptideSequence + ',' +
                                str(Data[2]) + ',' +
                                '{:.3%}'.format(SequenceFraction) + '\n')
    SortedSequenceFile.close    
#-------------------------------------------------------------------------------


#_____________________________RUNNING THE FUNCTION_____________________________#
#___SortedPeptideSequecesFileName, fastqFileLocation___
SelectionRoundSortedSequenceListGenerator('Nasir3Mar2016 PB1 R5.csv', '/Users/NikitaLoik/SugaLabDataAnalysis/Nasir3Mar2016/PB1/PB1 R5.fastq')
