# -*- coding: utf-8 -*-
#-------------------------------------------------------------------------------
def DNACodingSequence(DNASequence, QualityScoreSequence, StartSequence, StopSequence):
#utilises ONLY ONE StopSequence, returns ONLY ONE CodingSequence
    
    QualityScoreString = """!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"""
    ThresholdQualityScore = 29 # ThresholdQualityScore must be between 0 and 93
    ThresholdQualityString = QualityScoreString[ThresholdQualityScore:]
    
    MinLength = 24
    MaxLength = 240
            
    StartIndex = DNASequence.find(StartSequence) + len(StartSequence)
    StopIndex = DNASequence.rfind(StopSequence)
    CodingSequence =  DNASequence[StartIndex:StopIndex]
    if MinLength <= len(CodingSequence) and len(CodingSequence) <= MaxLength and len(CodingSequence)%3 == 0:
        for Character in QualityScoreSequence[StartIndex:StopIndex]:
            if Character not in ThresholdQualityString:
                return None
        return str(CodingSequence)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def Translation(CodingSequence):
#translates DNA sequence

    TranslationCode = {'AAA':'K','AAC':'N','AAG':'K','AAU':'N',
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
                    'UUA':'L','UUC':'F','UUG':'L','UUU':'F'}
    # UAA (ochre) — #
    # UAG (amber) — *
    # UGA (opal) — &
                    
    TranscriptionCode = {'A':'A','C':'C','G':'G','T':'U','U':'T'}
      
    RNASequence = ''
    for Nucleotide in CodingSequence:
        RNASequence = RNASequence + TranscriptionCode.get(Nucleotide,'X')
    #converts DNA to RNA
    #print RNASequence
        
    PeptideSequence = ''
    while len(RNASequence) != 0:
        PeptideSequence = PeptideSequence + TranslationCode.get(RNASequence[0:3],'Do not fuck with me!')
        RNASequence = RNASequence[3:]
    return PeptideSequence
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def SingleSelectionRoundSummary(fastqFileLocation):
#returns a list of lists with peptide-sequences and their frequencies, sorted by frequency in descending order
    
    RawDataFile = open(fastqFileLocation, 'r')
    Lines = RawDataFile.readlines()
    RawDataFile.close
    
    #StartSequence = 'ATG' # Met codon
    #StopSequence = 'TAG' # amber stop codon
    
    StartSequence = 'TAATACGACTCACTATAGGGTTAACTTTAAGAAGGAGATATACATATG'    # NNK - T7g10M.F48 
    StopSequence = 'TGCGGCAGCGGCAGCGGCAGCTAGGACGGGGGGCGGAAA' #NNK - CGS3an13.R39 
    # StartSequence = 'TAATACGACTCACTATAGGGTTGAACTTTAAGTAGGAGATATATCCATG'   #NNU - T7-CH-F49
    # StopSequence = 'TGTGGGTCTGGGTCTGGGTCTTAGGACGGGGGGCGGAAA'  #NNU - CGS3-CH-R39
    
    SingleSelectionRoundSummary = {}
    #creates empty SingleSelectionRoundSummary dictionary to store the results from a single round of selection
    #SingleSelectionRoundSummary = {PeptideY:    {CodingSequenceYZ:    OccurenceYZ}}
        
    #populates SingleSelectionRoundSummary dictionary with the results from a single round of selection
    for i,Line in enumerate(Lines):
        if StartSequence in Line and StopSequence in Line:
            CodingSequence = DNACodingSequence(Line, Lines[i + 2], StartSequence, StopSequence)
            if CodingSequence != None:
                PeptideSequence = Translation(CodingSequence)
                if PeptideSequence not in SingleSelectionRoundSummary:
                    SingleSelectionRoundSummary[str(PeptideSequence)] = {str(CodingSequence) : 1}
                else:
                    if CodingSequence not in SingleSelectionRoundSummary[str(PeptideSequence)]:
                        SingleSelectionRoundSummary[str(PeptideSequence)][str(CodingSequence)] = 1
                    else:
                        SingleSelectionRoundSummary[str(PeptideSequence)][str(CodingSequence)] += 1

    return SingleSelectionRoundSummary
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def HammingDistance(SequenceI, SequenceII):
    
    if len(SequenceI) < len(SequenceII):
        SequenceI = SequenceI + (len(SequenceII) - len(SequenceI)) * '%'
    elif len(SequenceI) > len(SequenceII):
        SequenceII = SequenceII + (len(SequenceI) - len(SequenceII)) * '%'
    
    HammingDistance = 0
    for i in range(len(SequenceI)):
        if SequenceI[i] == SequenceII[i]:
            HammingDistance = HammingDistance
        else:
            HammingDistance = HammingDistance + 1
            
    return HammingDistance
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def CompleteSelectionSummary(fastqDataFolderLocation):
# returns a SelectionSummary dictionary with the following structure
# SelectionSummary = {SelectionRoundX:    {PeptideSequenceXY:    {CodingSequenceXYZ:    OccurenceXYZ}}}

    CompleteSelectionSummary = {}
    # creates empty SelectionSummary dictionary to store the results from all the rounds of selection

    import os           
    for file in os.listdir(fastqDataFolderLocation):
        
        FileLocation = os.path.join(fastqDataFolderLocation, file)
          
        if file.endswith('.fastq'): # this conditional is necessary; without it some shit appears in the beginning of the file list
            RoundNumberFirstDigit = file[file.find('.')-2]
            RoundNumberSecondDigit = file[file.find('.')-1]
            if RoundNumberFirstDigit == '0':
                RoundNumber = int(RoundNumberSecondDigit)
                #print RoundNumber
            elif RoundNumberFirstDigit != '0':
                RoundNumber = int(file[file.find('.')-2 : file.find('.')])
                #print RoundNumber
        #(1.A) extracts the round number from the file name (file name should have two digit number before full stop — '00.') 
                
            SelectionRoundSummary = SingleSelectionRoundSummary(FileLocation)
            #(1.B) extracts single round results 
                    
            CompleteSelectionSummary[RoundNumber] = SelectionRoundSummary
            #(1.C) populate ConcatenatedResultsList
            #print ConcatenatedResultsList
            
    return CompleteSelectionSummary
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------    
def SelectionSummaryReport(fastqDataFolderLocation, BaseRoundIndex, TopPeptidesNumber, SelectionSummaryReportFileName):
    
    SelectionSummary = CompleteSelectionSummary(fastqDataFolderLocation)
    
    SortedRoundsList = sorted(SelectionSummary.keys()) #sorts the dictionary by keys, returns list of sorted keys
    #(1.A) creates sorted RoundsList
    
    DistinctPeptideSequencesNumbersByRound = {}
    for Round in SortedRoundsList:
        PeptideSequencesNumbersInRound = {}
        for PeptideSequence in SelectionSummary[Round]:
            PeptideSequencesNumbersInRound[PeptideSequence] = sum(SelectionSummary[Round][PeptideSequence].values())
        DistinctPeptideSequencesNumbersByRound[Round] = PeptideSequencesNumbersInRound
    #(1.B) create list of IndividualPeptideSequencesOccurenceByRound
    
    TotalPeptideSequencesNumbersByRound = {}
    for Round in SortedRoundsList:
        TotalPeptideSequencesNumbersByRound[Round] = sum(DistinctPeptideSequencesNumbersByRound[Round].values())
    #(1.C) create list of IndividualPeptideSequencesOccurenceByRound
    
            
    DistinctPeptideSequencesNumbersInBaseRound = DistinctPeptideSequencesNumbersByRound[BaseRoundIndex]
    SortedPeptidesSequencesList = sorted(DistinctPeptideSequencesNumbersInBaseRound, key = DistinctPeptideSequencesNumbersInBaseRound.get, reverse = True)
    print SortedPeptidesSequencesList
    SortedTopPeptidesSequencesList = SortedPeptidesSequencesList[0 : (TopPeptidesNumber)]
    # creates SortedTopPeptideSequencesList
    
    CloneNumberList = {}
    for PeptideSequence in SortedTopPeptidesSequencesList:
        PeptideSequenceOccurenceNumber = DistinctPeptideSequencesNumbersByRound[6].get(PeptideSequence, 0)
        PeptideSequenceDNANumber = len(SelectionSummary[6][PeptideSequence])
        CloneNumberList[PeptideSequence] = (PeptideSequenceDNANumber,PeptideSequenceOccurenceNumber)
    print CloneNumberList
        
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
    CSVSelectionSummaryReportFileName = SelectionSummaryReportFileName + '.csv'
    SelectionSummaryReportFile = open(CSVSelectionSummaryReportFileName, 'w')
    
    SelectionSummaryReportFile.write('peptide sequence' + ',')
    for Round in SortedRoundsList:
        SelectionSummaryReportFile.write('round # ' + str(Round) + ' occurence (#)' + ',')
    SelectionSummaryReportFile.write('\n')
    
    for PeptideSequence in SortedTopPeptidesSequencesList:
        SelectionSummaryReportFile.write(PeptideSequence + ',')
        for Round in SortedRoundsList:
            SelectionSummaryReportFile.write(str(DistinctPeptideSequencesNumbersByRound[Round].get(PeptideSequence, 0)) + ',')
        SelectionSummaryReportFile.write('\n')
        
    SelectionSummaryReportFile.write('total #' + ',')
    for Round in SortedRoundsList:
        SelectionSummaryReportFile.write(str(TotalPeptideSequencesNumbersByRound[Round]) + ',')
    SelectionSummaryReportFile.write('\n\n\n')
    
    SelectionSummaryReportFile.write('peptide sequence' + ',')
    for Round in SortedRoundsList:
        SelectionSummaryReportFile.write('round # ' + str(Round) + ' fraction (%)' + ',')
    SelectionSummaryReportFile.write('\n')
    
    for PeptideSequence in SortedTopPeptidesSequencesList:
        SelectionSummaryReportFile.write(PeptideSequence + ',')
        for Round in SortedRoundsList:
            PeptideSequenceFraction = float((DistinctPeptideSequencesNumbersByRound[Round].get(PeptideSequence, 0)))/float(TotalPeptideSequencesNumbersByRound[Round])
            SelectionSummaryReportFile.write('{:.3%}'.format(PeptideSequenceFraction) + ',')
        SelectionSummaryReportFile.write('\n')
            
    SelectionSummaryReportFile.close()
#-------------------------------------------------------------------------------
    import matplotlib.pyplot as plt
    
    plt.style.use('fivethirtyeight') # just to create 'ggplot' style
    
    for PeptideSequence in SortedTopPeptidesSequencesList:
        PeptideSequencesFractionsByRound = []
        for Round in SortedRoundsList:
            PeptideSequencesFractionsByRound += [float((DistinctPeptideSequencesNumbersByRound[Round].get(PeptideSequence, 0)))/float(TotalPeptideSequencesNumbersByRound[Round])]
        
        x = SortedRoundsList
        y = PeptideSequencesFractionsByRound
        plt.plot(x, y, 'o-', lw = 2, markersize = 3)
    
    plt.xlabel('Selection Round #', fontsize=14)
    plt.ylabel('Peptide Fraction', fontsize=14)
    
    legend = plt.legend(SortedTopPeptidesSequencesList, loc='upper center', bbox_to_anchor=(0.5, -0.15),
            fancybox=True, shadow=False, ncol=2)
    
    PNGFileName = SelectionSummaryReportFileName + '.png'
    plt.savefig(PNGFileName, bbox_extra_artists=[legend],bbox_inches='tight', dpi = 300)
    plt.show()
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def RelatedPeptideSequences(fastqDataFolderLocation, BaseRoundIndex):
    
    SelectionSummary = CompleteSelectionSummary(fastqDataFolderLocation)

    DistinctPeptideSequencesNumbersByRound = {}
    for Round in SelectionSummary:
        PeptideSequencesNumbersInRound = {}
        for PeptideSequence in SelectionSummary[Round]:
            PeptideSequencesNumbersInRound[PeptideSequence] = sum(SelectionSummary[Round][PeptideSequence].values())
        DistinctPeptideSequencesNumbersByRound[Round] = PeptideSequencesNumbersInRound
    
    DistinctPeptideSequencesNumbersInBaseRound = DistinctPeptideSequencesNumbersByRound[BaseRoundIndex]
    SortedPeptidesSequencesList = sorted(DistinctPeptideSequencesNumbersInBaseRound, key = DistinctPeptideSequencesNumbersInBaseRound.get, reverse = True)
        
    TopPeptideSequence = SortedPeptidesSequencesList[0]
    
    HammingDistanceThreshold = 2
    
    RelatedPeptideSequencesList = []
    
    for Round in SelectionSummary:
        PeptideSequencesDictionary = SelectionSummary[Round]
        for PeptideSequence in PeptideSequencesDictionary:
            if HammingDistance(TopPeptideSequence, PeptideSequence) <= HammingDistanceThreshold and PeptideSequence not in RelatedPeptideSequences:
                RelatedPeptideSequencesList += [PeptideSequence]
    return RelatedPeptideSequencesList
#-------------------------------------------------------------------------------

#_____________________________RUNNING THE FUNCTION_____________________________#

#___DataFolderLocation, BaseSelectionRoundNumber, TopNPeptidesNumber, SummaryFileName___
SelectionSummaryReport('/Users/NikitaLoik/SequenceTest', 6, 25, '6June2016Test NNK QS29')