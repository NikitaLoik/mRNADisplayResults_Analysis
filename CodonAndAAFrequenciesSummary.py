# -*- coding: utf-8 -*-
#-------------------------------------------------------------------------------
# utilise ONLY ONE stop codons, returns ONLY ONE ORF
def OpenReadingFrame(DNASequence):
    StartCodon = 'ATG'
    StopCodon = 'TAG'
    MinLength = 36
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
                    
                    'UAA':'+','UAC':'Y','UAG':'*','UAU':'Y',
                    'UCA':'S','UCC':'S','UCG':'S','UCU':'S',
                    'UGA':'#','UGC':'C','UGG':'W','UGU':'C',
                    'UUA':'L','UUC':'F','UUG':'L','UUU':'F'}
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
# return tuple of two lists: (1) CodonFrequenciesList and (2) AAFrequenciesList
def CodonAndAAFrequencies(fastqFileLocation):
    RawDataFile = open(fastqFileLocation, 'r')
    Lines = RawDataFile.readlines()
    RawDataFile.close
    CodonOccurences = {'AAA':0,'AAC':0,'AAG':0,'AAT':0,
                        'ACA':0,'ACC':0,'ACG':0,'ACT':0,
                        'AGA':0,'AGC':0,'AGG':0,'AGT':0,
                        'ATA':0,'ATC':0,'ATG':0,'ATT':0,
                        
                        'CAA':0,'CAC':0,'CAG':0,'CAT':0,
                        'CCA':0,'CCC':0,'CCG':0,'CCT':0,
                        'CGA':0,'CGC':0,'CGG':0,'CGT':0,
                        'CTA':0,'CTC':0,'CTG':0,'CTT':0,
                        
                        'GAA':0,'GAC':0,'GAG':0,'GAT':0,
                        'GCA':0,'GCC':0,'GCG':0,'GCT':0,
                        'GGA':0,'GGC':0,'GGG':0,'GGT':0,
                        'GTA':0,'GTC':0,'GTG':0,'GTT':0,
                        
                        'TAA':0,'TAC':0,'TAG':0,'TAT':0,
                        'TCA':0,'TCC':0,'TCG':0,'TCT':0,
                        'TGA':0,'TGC':0,'TGG':0,'TGT':0,
                        'TTA':0,'TTC':0,'TTG':0,'TTT':0}
                        
    AAOccurences = {'A':0,'C':0,'D':0,'E':0,
                    'F':0,'G':0,'H':0,'I':0,
                    'K':0,'L':0,'M':0,'N':0,
                    'P':0,'Q':0,'R':0,'S':0,
                    'T':0,'V':0,'W':0,'Y':0,
                    '*':0, '+':0, '#':0}
    
    CodonFrequenciesList = []
    TotalNumberOfCodons = 0
    CodonFrequencies = {}
    AAFrequenciesList =[]
    TotalNumberOfAA = 0
    AAFrequencies = {}
    
    SelectionRoundNumber = fastqFileLocation[fastqFileLocation.find('.')-1]
    
    # finde ORF in very line, count codons in every ORF, count every AA in every translated ORF
    for Line in Lines:
        ORF = OpenReadingFrame(Line)
        if ORF != None:
            Index = 0
            while Index < (len(ORF) - 3):
                CodonOccurences[ORF[Index:Index + 3]] += 1
                Index += 3
            Peptide = Translation(ORF)
            for AminoAcid in range(len(Peptide)-1):
#                print Peptide[AminoAcid]
                AAOccurences[Peptide[AminoAcid]] += 1
    
    # find Total Number of Codons Occurrences and Total Number of AA Occurences    
    for Codon in CodonOccurences:
        TotalNumberOfCodons += CodonOccurences[Codon]
    for AA in AAOccurences:
        TotalNumberOfAA += AAOccurences[AA]
    
    #transform Codons Occurences into Codons Frequencies and AA Occurences into AA Frequencies
    for Codon in CodonOccurences:
        CodonFrequencies[Codon] = CodonOccurences[Codon]/TotalNumberOfCodons
        
    for AA in AAOccurences:
        AAFrequencies[AA] = AAOccurences[AA]/TotalNumberOfAA
                        
    # convert the dictionary into the list of lists
    for key, value in CodonFrequencies.iteritems():
        CodonFrequenciesList.append([str(SelectionRoundNumber), key, value])
    for key, value in AAFrequencies.iteritems():
        AAFrequenciesList.append([str(SelectionRoundNumber), key, value])
               
    # sort the PeptideSequenceList by peptide sequence occurence in descendent order
    return (CodonFrequenciesList,TotalNumberOfCodons, AAFrequenciesList, TotalNumberOfAA)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def CodonAndAAFrequenciesSummary(DataFolderLocation, CodonAndAAFrequenciesFileName):

    import os
        
    # (1) create ConcatenatedCodonFrequenciesList and ConcatenatedAAFrequenciesList to store Results from all the rounds of selection       
    ConcatenatedCodonFrequenciesList = {}
    ConcotenatedTotalNumberOfCodonsList = {}
    ConcatenatedAAFrequenciesList = {}
    ConcatenatedTotalNumberOfAANumber = {}
    
    # open files in directory one by one
    for file in os.listdir(DataFolderLocation):
        CurrentFile = ''
        if file.endswith('.fastq'): # this conditional is necessary; without it some shit appears in the beginning of the file list
            CurrentFile = os.path.join(DataFolderLocation, file)
    # print CurrentFile
            
    # (1.A) extract round number from the file name        
            SelectionRoundNumber = CurrentFile[CurrentFile.find('.')-1]
    # print SelectionRoundNumber
                    
    # (1.B) populate ConcatenatedResultsList                
            ConcatenatedCodonFrequenciesList[SelectionRoundNumber] = CodonAndAAFrequencies(CurrentFile)[0]
            ConcatenatedAAFrequenciesList[SelectionRoundNumber] = CodonAndAAFrequencies(CurrentFile)[2]
            
    # (1.C) populate ConcatenatedCodonNumberList and ConcatenatedAANumberList
            ConcotenatedTotalNumberOfCodonsList[SelectionRoundNumber] = CodonAndAAFrequencies(CurrentFile)[1]
            ConcatenatedTotalNumberOfAANumber[SelectionRoundNumber] = CodonAndAAFrequencies(CurrentFile)[3]
            
            
    # print ConcatenatedResultsList
    
    # (2.A) create sorted lists ConcatenatedCodonFrequenciesList, ConcotenatedCodonNumberList, ConcatenatedAAFrequenciesList, ConcatenatedAANumberList
    RoundsList = sorted(ConcatenatedCodonFrequenciesList.keys())
    SortedConcatenatedCodonFrequenciesList = sorted(ConcatenatedCodonFrequenciesList.keys())
    SortedConcotenatedTotalNumberOfCodonsList = sorted(ConcotenatedTotalNumberOfCodonsList.keys())
    SortedConcatenatedAAFrequenciesList = sorted(ConcatenatedAAFrequenciesList.keys())
    SortedConcatenatedTotalNumberOfAANumber = sorted(ConcatenatedTotalNumberOfAANumber.keys())

    #-------------------------------------------------------------------------------
    
    #-------------------------------------------------------------------------------
    CodonsList = ['AAA','AAC','AAG','AAT',
                'ACA','ACC','ACG','ACT',
                'AGA','AGC','AGG','AGT',
                'ATA','ATC','ATG','ATT',
            
                'CAA','CAC','CAG','CAT',
                'CCA','CCC','CCG','CCT',
                'CGA','CGC','CGG','CGT',
                'CTA','CTC','CTG','CTT',
                        
                'GAA','GAC','GAG','GAT',
                'GCA','GCC','GCG','GCT',
                'GGA','GGC','GGG','GGT',
                'GTA','GTC','GTG','GTT',
                        
                'TAA','TAC','TAG','TAT',
                'TCA','TCC','TCG','TCT',
                'TGA','TGC','TGG','TGT',
                'TTA','TTC','TTG','TTT']
    AAList = ['A','C','D','E',
            'F','G','H','I',
            'K','L','M','N',
            'P','Q','R','S',
            'T','V','W','Y',
            '§', '*', '#']
    
    
    CSVFileName = CodonAndAAFrequenciesFileName + '.csv'
    CodonAndAAFrequenciesFile = open(CSVFileName, 'w')
    
    CodonAndAAFrequenciesFile.write('codon' + ',')
    for Round in RoundsList:
        CodonAndAAFrequenciesFile.write('round # ' + Round + ' abundance (%)' + ',')
    CodonAndAAFrequenciesFile.write('\n')
    
    for Codon in CodonsList:
        CodonAndAAFrequenciesFile.write(Codon + ',')
        for Round in SortedConcatenatedCodonFrequenciesList:
            CodonAndAAFrequenciesFile.write('{:.3%}'.format(SortedConcatenatedCodonFrequenciesList[Round][Codon]) + ',')
        CodonAndAAFrequenciesFile.write('\n')
        
    CodonAndAAFrequenciesFile.write('total #' + ',')
    for Round in RoundsList:
        CodonAndAAFrequenciesFile.write(str(SortedConcotenatedTotalNumberOfCodonsList[Round]) + ',')
    CodonAndAAFrequenciesFile.write('\n\n\n')
    
    
    CodonAndAAFrequenciesFile.write('AA' + ',')
    for Round in RoundsList:
        CodonAndAAFrequenciesFile.write('round # ' + Round + ' abundance (%)' + ',')
    CodonAndAAFrequenciesFile.write('\n')
    
    for AA in AAList:
        CodonAndAAFrequenciesFile.write(Codon + ',')
        for Round in SortedConcatenatedAAFrequenciesList:
            CodonAndAAFrequenciesFile.write('{:.3%}'.format(SortedConcatenatedAAFrequenciesList[Round][AA]) + ',')
        CodonAndAAFrequenciesFile.write('\n')
        
    CodonAndAAFrequenciesFile.write('total #' + ',')
    for Round in RoundsList:
        CodonAndAAFrequenciesFile.write(str(SortedConcatenatedTotalNumberOfAANumber[Round]) + ',')
    CodonAndAAFrequenciesFile.write('\n')
            
    CodonAndAAFrequenciesFile.close()
    #-------------------------------------------------------------------------------
    
    #-------------------------------------------------------------------------------
#    import matplotlib.pyplot as plt
    
#    for PeptideNumber in range(len(TopNPeptidesList)):
#        plt.plot(RoundsList, ListOfPeptidesFractionsByRound[PeptideNumber])
    
#    plt.xlabel('Selection Round #', fontsize=14)
#    plt.ylabel('Peptide Fraction', fontsize=14)
#    legend = plt.legend(TopNPeptidesList, loc='upper center', bbox_to_anchor=(0.5, -0.15),
#            fancybox=True, shadow=False, ncol=2)
#    PNGFileName = SummaryFileName + '.png'
#    plt.savefig(PNGFileName, bbox_extra_artists=[legend],bbox_inches='tight', dpi = 300)
#    plt.show()
#-------------------------------------------------------------------------------


#_____________________________RUNNING THE FUNCTION_____________________________#

#___DataFolderLocation, CodonAndAAFrequenciesFileName___
CodonAndAAFrequenciesSummary('/Users/NikitaLoik/NST', 'NasirNSTTest')

