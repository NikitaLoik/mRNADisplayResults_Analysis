# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import math
import datetime
import networkx as nx
################################################################################
################################################################################
def TodaysDate():
        
    today = datetime.date.today()
    TodaysDate = today.strftime('%d%b%Y')
    
    return TodaysDate

################################################################################
################################################################################
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
################################################################################
################################################################################
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
        RNASequence += TranscriptionCode.get(Nucleotide,'X')
    #converts DNA to RNA
    #print RNASequence
        
    Peptide = ''
    while len(RNASequence) != 0:
        Peptide += TranslationCode.get(RNASequence[0:3],'Do not fuck with me!')
        RNASequence = RNASequence[3:]
    return Peptide
################################################################################
################################################################################
def SingleSelectionRoundSummary(fastqFileLocation):
#returns a list of lists with peptide-sequences and their frequencies, sorted by frequency in descending order
    
    RawDataFile = open(fastqFileLocation, 'r')
    Lines = RawDataFile.readlines()
    RawDataFile.close
    
    #StartSequence = 'ATG' # Met codon
    #StopSequence = 'TAG' # amber stop codon
    
    StartSequence = 'TAATACGACTCACTATAGGGTTAACTTTAAGAAGGAGATATACATATG'    # NNK - T7g10M.F48 
    StopSequence = 'TGCGGCAGCGGCAGCGGCAGCTAGGACGGGGGGCGGAAA' #NNK - CGS3an13.R39 
    #StartSequence = 'TAATACGACTCACTATAGGGTTGAACTTTAAGTAGGAGATATATCCATG'   #NNU - T7-CH-F49
    #StopSequence = 'TGTGGGTCTGGGTCTGGGTCTTAGGACGGGGGGCGGAAA'  #NNU - CGS3-CH-R39
    
    SingleSelectionRoundSummary = {}
    #creates empty SingleSelectionRoundSummary dictionary to store the results from a single round of selection
    #SingleSelectionRoundSummary = {PeptideY:    {CodingSequence_YZ:    Occurence_YZ}}
        
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
################################################################################
################################################################################
def HammingDistance(Sequence1, Sequence2):
    
    if len(Sequence1) < len(Sequence2):
        Sequence1 = Sequence1 + (len(Sequence2) - len(Sequence1)) * '%'
    elif len(Sequence1) > len(Sequence2):
        Sequence2 = Sequence2 + (len(Sequence1) - len(Sequence2)) * '%'
    
    HammingDistance = 0
    for i in range(len(Sequence1)):
        if Sequence1[i] == Sequence2[i]:
            HammingDistance = HammingDistance
        else:
            HammingDistance = HammingDistance + 1
            
    return HammingDistance
################################################################################
################################################################################    
def HammingDistanceBasedFormating(Sequence1, Sequence2):
    
    if len(Sequence1) < len(Sequence2):
        Sequence1 = Sequence1 + (len(Sequence2) - len(Sequence1)) * '-'
    elif len(Sequence1) > len(Sequence2):
        Sequence2 = Sequence2 + (len(Sequence1) - len(Sequence2)) * '-'
    
    HammingDistance = 0
    FormatedSequence2 = ''
    for i in range(len(Sequence1)):
        if Sequence1[i] == Sequence2[i]:
            FormatedSequence2 += Sequence2[i].lower()
            HammingDistance = HammingDistance
        else:
            FormatedSequence2 += Sequence2[i]
            HammingDistance = HammingDistance + 1            
    return FormatedSequence2
################################################################################
################################################################################
def CompleteSelectionSummary(fastqDataFolderLocation):
# returns a SelectionSummary dictionary with the following structure
# SelectionSummary = {SelectionRound_X:    {PeptideXY:    {CodingDNA_XYZ:    Occurence_XYZ}}}

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
################################################################################
################################################################################   
def PeptidesOccurences_BY_Round(fastqDataFolderLocation):
    SelectionSummary = CompleteSelectionSummary(fastqDataFolderLocation)
    
    PeptidesOccurences_BY_Round = {}
    for Round in SelectionSummary:
        PeptidesOccurences_IN_Round = {}
        for Peptide in SelectionSummary[Round]:
            PeptidesOccurences_IN_Round[Peptide] = sum(SelectionSummary[Round][Peptide].values())
        PeptidesOccurences_BY_Round[Round] = PeptidesOccurences_IN_Round
        
    return PeptidesOccurences_BY_Round
################################################################################
def DNAsOccurences_BY_Round(fastqDataFolderLocation):
    SelectionSummary = CompleteSelectionSummary(fastqDataFolderLocation)
    
    DNAsOccurences_BY_Round = {}
    for Round in SelectionSummary:
        DNAsOccurences_IN_Round = {}
        for Peptide in SelectionSummary[Round]:
            for DNA in SelectionSummary[Round][Peptide]:
                DNAsOccurences_IN_Round[DNA] = SelectionSummary[Round][Peptide][DNA]
        DNAsOccurences_BY_Round[Round] = DNAsOccurences_IN_Round

    return DNAsOccurences_BY_Round
################################################################################ 
def TotalReads_BY_Round(fastqDataFolderLocation):
    SelectionSummary = CompleteSelectionSummary(fastqDataFolderLocation)
    Peptides_BY_Round = PeptidesOccurences_BY_Round(fastqDataFolderLocation)
    
    TotalReads_BY_Round = {}
    for Round in SelectionSummary:
        TotalReads_BY_Round[Round] = sum(Peptides_BY_Round[Round].values())
        
    return TotalReads_BY_Round
################################################################################
################################################################################ 
def BaseRoundSortedPeptidesList(fastqDataFolderLocation, BaseRoundIndex):
    Peptides_BY_Round = PeptidesOccurences_BY_Round(fastqDataFolderLocation)  
            
    PeptidesOccurencesInBaseRound = Peptides_BY_Round[BaseRoundIndex]
    BaseRoundSortedPeptidesList = sorted(PeptidesOccurencesInBaseRound, key = PeptidesOccurencesInBaseRound.get, reverse = True)
    
    return BaseRoundSortedPeptidesList
################################################################################
################################################################################
def BaseRoundSortedDNAsList(fastqDataFolderLocation, BaseRoundIndex):
    DNAs_BY_Round = DNAsOccurences_BY_Round(fastqDataFolderLocation)  
            
    DNAsOccurences_IN_BaseRound = DNAs_BY_Round[BaseRoundIndex]
    BaseRoundSortedDNAsList = sorted(DNAsOccurences_IN_BaseRound, key = DNAsOccurences_IN_BaseRound.get, reverse = True)
    
    return BaseRoundSortedDNAsList
################################################################################
################################################################################
def DNAClonesOccurences_BY_Round_BY_Peptide(fastqDataFolderLocation):
    SelectionSummary = CompleteSelectionSummary(fastqDataFolderLocation)
    
    DNAClonesOccurences_BY_Round_BY_Peptide = {}
    for Round in SelectionSummary:
        DNAClonesOccurences_BY_Peptide = {}
        for Peptide in SelectionSummary[Round]:
            DNAClonesOccurences_BY_Peptide[Peptide] = len(SelectionSummary[Round][Peptide])
        DNAClonesOccurences_BY_Round_BY_Peptide[Round] = DNAClonesOccurences_BY_Peptide
        
    return DNAClonesOccurences_BY_Round_BY_Peptide
################################################################################
################################################################################
def PeptidesAppearances_BY_Round(BaseRoundSortedPeptidesList, PeptidesOccurences_BY_Round):
    
    PeptidesAppearances_BY_Round = {}
    
    for Peptide in BaseRoundSortedPeptidesList:
        PeptidesAppearances_BY_Round[Peptide] = []
        for Round in PeptidesOccurences_BY_Round:
            if Peptide in PeptidesOccurences_BY_Round[Round]:
                PeptidesAppearances_BY_Round[Peptide] += [Round]
    return PeptidesAppearances_BY_Round
################################################################################
################################################################################      
def DNAsAppearances_BY_Round(BaseRoundSortedDNAsList, DNAsOccurences_BY_Round):
    
    DNAsAppearances_BY_Round = {}
    
    for DNA in BaseRoundSortedDNAsList:
        DNAsAppearances_BY_Round[DNA] = []
        for Round in DNAsOccurences_BY_Round:
            if DNA in DNAsOccurences_BY_Round[Round]:
                DNAsAppearances_BY_Round[DNA] += [Round]
    return DNAsAppearances_BY_Round
################################################################################
################################################################################  
def SelectionSummaryReport(fastqDataFolderLocation, BaseRoundIndex, NumberOfTopPeptides, SelectionSummaryReportFileName):
    
    today = TodaysDate() 
    
    SelectionSummaryFileNameCSV = str(today) + 'SelectionSummary' + SelectionSummaryReportFileName + '.csv'
    SelectionSummaryReportFile = open(SelectionSummaryFileNameCSV, 'w')
    
    SelectionSummary = CompleteSelectionSummary(fastqDataFolderLocation)
    SortedRoundsList = sorted(SelectionSummary.keys())
    
    Peptides_BY_Round = PeptidesOccurences_BY_Round(fastqDataFolderLocation)
    TotalPeptides_BY_Round = TotalReads_BY_Round(fastqDataFolderLocation)
    
    BaseRoundSortedPeptides = BaseRoundSortedPeptidesList(fastqDataFolderLocation, BaseRoundIndex)
    BaseRoundTopSortedPeptides = BaseRoundSortedPeptides[0 : (NumberOfTopPeptides)] 
    
    SelectionSummaryReportFile.write('peptide sequence' + ',')
    for Round in SortedRoundsList:
        SelectionSummaryReportFile.write('round # ' + str(Round) + ' occurence (#)' + ',')
    SelectionSummaryReportFile.write('\n')
    
    for Peptide in BaseRoundTopSortedPeptides:
        SelectionSummaryReportFile.write(Peptide + ',')
        for Round in SortedRoundsList:
            SelectionSummaryReportFile.write(str(Peptides_BY_Round[Round].get(Peptide, 0)) + ',')
        SelectionSummaryReportFile.write('\n')
        
    SelectionSummaryReportFile.write('total #' + ',')
    for Round in SortedRoundsList:
        SelectionSummaryReportFile.write(str(TotalPeptides_BY_Round[Round]) + ',')
    SelectionSummaryReportFile.write('\n\n\n')
    
    SelectionSummaryReportFile.write('peptide sequence' + ',')
    for Round in SortedRoundsList:
        SelectionSummaryReportFile.write('round # ' + str(Round) + ' fraction (%)' + ',')
    SelectionSummaryReportFile.write('\n')
    
    for Peptide in BaseRoundTopSortedPeptides:
        SelectionSummaryReportFile.write(Peptide + ',')
        for Round in SortedRoundsList:
            PeptideFraction = float((Peptides_BY_Round[Round].get(Peptide, 0)))/float(TotalPeptides_BY_Round[Round])
            SelectionSummaryReportFile.write('{:.3%}'.format(PeptideFraction) + ',')
        SelectionSummaryReportFile.write('\n')
            
    SelectionSummaryReportFile.close()
    
#-------------------------------------------------------------------------------
   
    plt.style.use('fivethirtyeight') # just to create 'ggplot' style
    
    Xs = []
    Ys = []
    for Peptide in BaseRoundTopSortedPeptides:
        PeptidesFractions_BY_Round = []
        for Round in SortedRoundsList:
            PeptidesFractions_BY_Round += [float((Peptides_BY_Round[Round].get(Peptide, 0)))/float(TotalPeptides_BY_Round[Round])]
        
        x = SortedRoundsList
        y = PeptidesFractions_BY_Round
        Xs += x
        Ys += y
        
        plt.plot(x, y,
                    'o-',
                    lw = 2.0,
                    ms = 4.0,
                    mew = 0.1,
                    mec = '#191919')
    
    XMin = min(Xs) - 0.05*(max(Xs) - min(Xs))
    XMax = max(Xs) + 0.05*(max(Xs) - min(Xs))
    YMin = min(Ys) - 0.05*(max(Ys) - min(Ys))
    YMax = max(Ys) + 0.05*(max(Ys) - min(Ys))
    
    plt.axis([XMin, XMax, YMin, YMax])
    
    plt.xlabel('Selection Round #', fontsize=14)
    plt.ylabel('Peptide Fraction', fontsize=14)
    
    legend = plt.legend(BaseRoundTopSortedPeptides, loc='upper center', bbox_to_anchor=(0.5, -0.15),
                        fancybox=True, shadow=False, ncol=2)
    
    SelectionSummaryFileNamePNG = str(today) + 'SelectionSummary' + SelectionSummaryReportFileName + '.png'
    
    plt.savefig(SelectionSummaryFileNamePNG, bbox_extra_artists=[legend],bbox_inches='tight', dpi = 300)
    plt.show()
    plt.close()
################################################################################
################################################################################    
def DNAMutantsAnalysis(fastqDataFolderLocation, BaseRoundIndex, NumberOfTopPeptides, SelectionSummaryReportFileName):
    
    today = TodaysDate() 
    
    DNAMutantsAnalysisFileNameCSV =  str(today) + 'DNAsMutantsAnalysis' + SelectionSummaryReportFileName + '.csv'
    DNAMutantsAnalysisFile = open(DNAMutantsAnalysisFileNameCSV, 'w')
    
    SelectionSummary = CompleteSelectionSummary(fastqDataFolderLocation)
    SortedRoundsList = sorted(SelectionSummary.keys())
    
    Peptides_BY_Round = PeptidesOccurences_BY_Round(fastqDataFolderLocation)
    TotalPeptides_BY_Round = TotalReads_BY_Round(fastqDataFolderLocation)
    
    BaseRoundSortedPeptides = BaseRoundSortedPeptidesList(fastqDataFolderLocation, BaseRoundIndex)
    BaseRoundTopSortedPeptides = BaseRoundSortedPeptides[0 : (NumberOfTopPeptides)]
    
    DNAClones_BY_Round_BY_Peptide = DNAClonesOccurences_BY_Round_BY_Peptide(fastqDataFolderLocation)
    
    DNAMutantsAnalysisFile.write('peptide sequence' + ',')
    for Round in SortedRoundsList:
        DNAMutantsAnalysisFile.write('round # ' + str(Round) + ' DNA clones (#)' + ',')
    DNAMutantsAnalysisFile.write('\n')
    
    for Peptide in BaseRoundTopSortedPeptides:
        DNAMutantsAnalysisFile.write(Peptide + ',')
        for Round in SortedRoundsList:
            DNAMutantsAnalysisFile.write(str(DNAClones_BY_Round_BY_Peptide[Round].get(Peptide, 0)) + ',')
        DNAMutantsAnalysisFile.write('\n')
        
#    DNAClonesSummaryReportFile.write('total #' + ',')
#    for Round in SortedRoundsList:
#        DNAClonesSummaryReportFile.write(str(TotalPeptides_BY_Round[Round]) + ',')
#    DNAClonesSummaryReportFile.write('\n\n\n')
    
#    DNAClonesSummaryReportFile.write('peptide sequence' + ',')
#    for Round in SortedRoundsList:
#        DNAClonesSummaryReportFile.write('round # ' + str(Round) + ' fraction (%)' + ',')
#    DNAClonesSummaryReportFile.write('\n')
    
#    for Peptide in BaseRoundTopSortedPeptides:
#        DNAClonesSummaryReportFile.write(Peptide + ',')
#        for Round in SortedRoundsList:
#            PeptideFraction = float((Peptides_BY_Round[Round].get(Peptide, 0)))/float(TotalPeptides_BY_Round[Round])
#            DNAClonesSummaryReportFile.write('{:.3%}'.format(PeptideFraction) + ',')
#        DNAClonesSummaryReportFile.write('\n')
            
    DNAMutantsAnalysisFile.close()
    
#-------------------------------------------------------------------------------        
    plt.style.use('fivethirtyeight') # just to create 'ggplot' style
    
    Xs = []
    Ys = []    
    for Peptide in BaseRoundTopSortedPeptides:
        SortedDNAClones_BY_Peptide = []
        for Round in SortedRoundsList:
            SortedDNAClones_BY_Peptide += [DNAClones_BY_Round_BY_Peptide[Round].get(Peptide, 0)]
        
        x = SortedRoundsList
        y = SortedDNAClones_BY_Peptide
        Xs += x
        Ys += y
        plt.plot(x, y,
                'o-',
                lw = 2.0,
                ms = 4.0,
                mew = 0.1,
                mec = '#191919')
                
    XMin = min(Xs) - 0.05*(max(Xs) - min(Xs))
    XMax = max(Xs) + 0.05*(max(Xs) - min(Xs))
    YMin = min(Ys) - 0.05*(max(Ys) - min(Ys))
    YMax = max(Ys) + 0.05*(max(Ys) - min(Ys))
    
    plt.axis([XMin, XMax, YMin, YMax])
    
    plt.xlabel('Selection Round #', fontsize=14)
    plt.ylabel('DNA Clones #', fontsize=14)
    
    legend = plt.legend(BaseRoundTopSortedPeptides, loc='upper center', bbox_to_anchor=(0.5, -0.15),
                        fancybox=True, shadow=False, ncol=2)
    
    DNAMutantsAnalysisFileNamePNG = str(today) + 'DNAsMutantsAnalysis' + SelectionSummaryReportFileName + '.png'
    
    plt.savefig(DNAMutantsAnalysisFileNamePNG, bbox_extra_artists=[legend],bbox_inches='tight', dpi = 300)
    plt.show()
    plt.close()
#-------------------------------------------------------------------------------    
    plt.style.use('fivethirtyeight')
    
#    PeptideDNAClonesNumber_IN_BaseRound = []
#    PeptideOccurence_IN_BaseRound = []

    RoundIndex = BaseRoundIndex

    Xs = []
    Ys = []        
    for Peptide in DNAClones_BY_Round_BY_Peptide[RoundIndex]:
        PeptideDNAClonesNumber_IN_BaseRound = math.log(DNAClones_BY_Round_BY_Peptide[RoundIndex].get(Peptide, 0), 2)
        PeptideOccurence_IN_BaseRound = math.log(Peptides_BY_Round[RoundIndex].get(Peptide, 0), 2)
    
        x = PeptideDNAClonesNumber_IN_BaseRound
        y = PeptideOccurence_IN_BaseRound
        Xs += [x]
        Ys += [y]
        
        plt.plot(x, y,
                'o',
                ms = 5.0,
                mew = 0.1,
                mec = '#191919')

    XMin = min(Xs) - 0.05*(max(Xs) - min(Xs))
    XMax = max(Xs) + 0.05*(max(Xs) - min(Xs))
    YMin = min(Ys) - 0.05*(max(Ys) - min(Ys))
    YMax = max(Ys) + 0.05*(max(Ys) - min(Ys))
    
    plt.axis([XMin, XMax, YMin, YMax])
        
    XLabel = 'log$2$ (DNA Clones #)' #$_$ makes subscript possible
    plt.xlabel(XLabel, fontsize = 14)
    YLabel = 'log$2$ (Peptide Occurence)'
    plt.ylabel(YLabel, fontsize = 14)
    
#    legend = plt.legend(BaseRoundTopSortedPeptides, loc='upper center', bbox_to_anchor=(0.5, -0.15),
#                        fancybox=True, shadow=False, ncol=2)
    
    DNAClonesAnalysisFileNamePNG = str(today) + 'DNAsMutantsAnalysisRegression' + 'R' + str(Round) + SelectionSummaryReportFileName + '.png'
    
    plt.savefig(DNAClonesAnalysisFileNamePNG, bbox_extra_artists=[legend], bbox_inches='tight', dpi = 300)
    plt.show()
    plt.close()
    
################################################################################
################################################################################      
def PeptidesMutantsSummaryReport(fastqDataFolderLocation, BaseRoundIndex, PeptidesClonesSummaryReportFileName):
    
    today = TodaysDate()
    
    Peptides_BY_Round = PeptidesOccurences_BY_Round(fastqDataFolderLocation)
    TotalPeptides_BY_Round = TotalReads_BY_Round(fastqDataFolderLocation)

    BaseRoundSortedPeptides = BaseRoundSortedPeptidesList(fastqDataFolderLocation, BaseRoundIndex)
    PeptidesAppearances = PeptidesAppearances_BY_Round(BaseRoundSortedPeptides, Peptides_BY_Round)
    
    BaseRoundPeptidesForest = nx.Graph()
    for Peptide in BaseRoundSortedPeptides:
        BaseRoundPeptidesForest.add_node(Peptide,
                                        Occurence = Peptides_BY_Round[BaseRoundIndex][Peptide],
                                        FirstAppearance = min(PeptidesAppearances[Peptide]))
    # adds nodes with attributes of Occurence in the BaseRound and FirstAppearance to the BaseRoundPeptidesForest
    UsedNodes = []
    for Peptide1 in BaseRoundSortedPeptides:
        UsedNodes += [Peptide1]
        for Peptide2 in BaseRoundSortedPeptides:
            if Peptide2 not in UsedNodes and HammingDistance(Peptide1, Peptide2) == 1:
                BaseRoundPeptidesForest.add_edge(Peptide1,Peptide2,
                                                MutationsNumber = 1)
    # adds edges to the BaseRoundPeptidesForest so that it can be made into a directed graph
    
    BaseRoundPeptidesTrees = list(nx.connected_component_subgraphs(BaseRoundPeptidesForest, copy = True))
    
    PeptidesClonesSummaryReportFileNameCSV =  str(today) + 'PeptidesMutantsSummary' + PeptidesClonesSummaryReportFileName + '.csv'
    PeptidesClonesSummaryReportFile = open(PeptidesClonesSummaryReportFileNameCSV, 'w')
    
    Positions = {}
    X_0_Coordinate = 1
    Y_0_Coordinate = 0
    Y_X0_Coordinate = 0
    
    for Tree in BaseRoundPeptidesTrees:
                        
        RootPeptide = max(nx.get_node_attributes(Tree, 'Occurence'),
                        key = lambda(Peptide): Tree[Peptide])
        
        TreePeptides = {}
        TreePeptides[RootPeptide] = [0, '', 0, Tree.node[RootPeptide]['Occurence'], Tree.node[RootPeptide]['FirstAppearance']]
        TreePeptidesList = Tree.nodes()
        TreePeptidesList.remove(RootPeptide)
        
        for Peptide in TreePeptidesList:
            PeptidePredecessor = nx.shortest_path(Tree, source = Peptide, target = RootPeptide, weight = None)[1]
            PredecessorOccurence = Tree.node[PeptidePredecessor]['Occurence']
            PeptideOccurence = Tree.node[Peptide]['Occurence']
            
            TreePeptides[Peptide] = [PeptidePredecessor, PredecessorOccurence, PeptideOccurence]
        
        Peptides_BY_DistanceToTheRoot = {}
        for Peptide in Tree.nodes():
            DistanceToTheRoot = nx.shortest_path_length(Tree, source = Peptide, target = RootPeptide, weight = None)
            if DistanceToTheRoot not in Peptides_BY_DistanceToTheRoot:
                Peptides_BY_DistanceToTheRoot[DistanceToTheRoot] = [Peptide]
            else:
                Peptides_BY_DistanceToTheRoot[DistanceToTheRoot] += [Peptide]
                
        MaxPeptideMutantsNumber = max(map(lambda(k): len(Peptides_BY_DistanceToTheRoot[k]), Peptides_BY_DistanceToTheRoot))
        
        SortedPeptides_BY_DistanceToTheRoot = {}
        
        for DistanceToTheRoot in Peptides_BY_DistanceToTheRoot:
            
            EquidistantPeptides = Peptides_BY_DistanceToTheRoot[DistanceToTheRoot]
            
            EquidistantPeptides = sorted(EquidistantPeptides, key = lambda (Peptide): (TreePeptides[Peptide][2]), reverse = True)
#            EquidistantPeptides = sorted(EquidistantPeptides, key = lambda (Peptide): (TreePeptides[Peptide][1]), reverse = True)
            EquidistantPeptides = sorted(EquidistantPeptides, key = lambda (Peptide): (TreePeptides[Peptide][0]), reverse = False)
            
            AdditionalElements = MaxPeptideMutantsNumber - len(EquidistantPeptides)
            SortedPeptides_BY_DistanceToTheRoot[DistanceToTheRoot] = EquidistantPeptides + AdditionalElements * ['']
            
            if len(Tree.nodes()) > 1:
                for Peptide in EquidistantPeptides:
                    XCoordinate = X_0_Coordinate + DistanceToTheRoot
                    YCoordinate = Y_0_Coordinate - EquidistantPeptides.index(Peptide)
                    Positions[Peptide] = (XCoordinate, YCoordinate)
            elif len(Tree.nodes()) == 1:
                for Peptide in EquidistantPeptides:
                    XCoordinate = 0
                    YCoordinate = Y_X0_Coordinate
                    Positions[Peptide] = (XCoordinate, YCoordinate)
        
                            
        if len(Tree.nodes()) > 1:
            X_0_Coordinate += max(Peptides_BY_DistanceToTheRoot.keys()) + 1
        
        if len(Tree.nodes()) == 1:
            Y_X0_Coordinate -= 1
            
            
#        print SortedPeptides_BY_DistanceToTheRoot
        

#-------------------------------------------------------------------------------

        if len(Tree.nodes()) > 1:
            for DistanceToTheRoot in SortedPeptides_BY_DistanceToTheRoot:
                PeptidesClonesSummaryReportFile.write(str(DistanceToTheRoot) + ' mutations' + ',' + 'frequency' + ',' + 'clones #' + ',' )
            PeptidesClonesSummaryReportFile.write('\n')
        
        for i in range(MaxPeptideMutantsNumber):
            for MutationsNumber in SortedPeptides_BY_DistanceToTheRoot:                        
                Peptide = SortedPeptides_BY_DistanceToTheRoot[MutationsNumber][i]
                
                if Peptide != '':
                    FormatedPeptide = HammingDistanceBasedFormating(RootPeptide, Peptide)
                    ClonesNumber = str(len(Tree.neighbors(Peptide)))
                    PeptideFraction = ('{:.2%}'.format(float((Peptides_BY_Round[BaseRoundIndex].get(Peptide, 0)))/float(TotalPeptides_BY_Round[BaseRoundIndex])))
                else:
                    FormatedPeptide = ''
                    ClonesNumber = ''
                    PeptideFraction = ''
                                                        
                PeptidesClonesSummaryReportFile.write(FormatedPeptide + ',' +
                            PeptideFraction + ',' +
                            ClonesNumber + ',')
            PeptidesClonesSummaryReportFile.write('\n')
        PeptidesClonesSummaryReportFile.write('\n')
                    
    PeptidesClonesSummaryReportFile.close()
        
    BaseRoundPeptidesGraph = nx.Graph()    
    BaseRoundPeptidesGraph.add_nodes_from(BaseRoundSortedPeptides)
     
    Sizes = []
    for Peptide in BaseRoundPeptidesGraph.nodes():
        Sizes.append(5 * math.log(Peptides_BY_Round[BaseRoundIndex][Peptide], 10) + 2)
    Colours = []
    for Peptide in BaseRoundPeptidesGraph.nodes():
        Colours.append(min(PeptidesAppearances[Peptide]))
    
#    print(BaseRoundPeptidesGraph.nodes())
    i = 0
    for Peptide1 in BaseRoundSortedPeptides:
        for Peptide2 in BaseRoundSortedPeptides:
            if HammingDistance(Peptide1, Peptide2) == 1:
                BaseRoundPeptidesGraph.add_edge(Peptide1,Peptide2)
        i += 1

    XMin = min(map(lambda(Peptide): Positions[Peptide][0], Positions)) - 1
    XMax = max(map(lambda(Peptide): Positions[Peptide][0], Positions)) + 1
    YMin = min(map(lambda(Peptide): Positions[Peptide][1], Positions)) - 1
    YMax = max(map(lambda(Peptide): Positions[Peptide][1], Positions)) + 1
    
    nx.draw_networkx(BaseRoundPeptidesGraph,
                    pos = Positions,
                    node_size = Sizes,
                    node_color = Colours,
                    cmap = 'Paired',
                    linewidths = 0.2,
                    width = 0.2,
                    with_labels = False,
                    font_size = 6)

    plt.axis('off')
    plt.axis([XMin, XMax, YMin, YMax])
    
    PeptidesMutantsSummaryReportFileNamePNG = str(today) + 'PeptidesMutantsSummary' + PeptidesClonesSummaryReportFileName + '.png'
    plt.savefig(PeptidesMutantsSummaryReportFileNamePNG, dpi = 500)
    
    plt.show()
    plt.close()
    
################################################################################
################################################################################      
def DNAsMutantsSummaryReport(fastqDataFolderLocation, BaseRoundIndex, DNAsClonesSummaryReportFileName):
    
    today = TodaysDate()
    
    DNAs_BY_Round = DNAsOccurences_BY_Round(fastqDataFolderLocation)
    TotalDNAs_BY_Round = TotalReads_BY_Round(fastqDataFolderLocation)

    BaseRoundSortedDNAs = BaseRoundSortedDNAsList(fastqDataFolderLocation, BaseRoundIndex)
    DNAsAppearances = DNAsAppearances_BY_Round(BaseRoundSortedDNAs, DNAs_BY_Round)
    
    BaseRoundSortedPeptides = BaseRoundSortedPeptidesList(fastqDataFolderLocation, BaseRoundIndex)
    
    BaseRoundDNAsForest = nx.Graph()
    for DNA in BaseRoundSortedDNAs:
        BaseRoundDNAsForest.add_node(DNA,
                                        Occurence = DNAs_BY_Round[BaseRoundIndex][DNA],
                                        FirstAppearance = min(DNAsAppearances[DNA]))
    # adds nodes with attributes of Occurence in the BaseRound and FirstAppearance to the BaseRoundDNAsForest
    UsedNodes = []
    for DNA1 in BaseRoundSortedDNAs:
        UsedNodes += [DNA1]
        for DNA2 in BaseRoundSortedDNAs:
            if DNA2 not in UsedNodes and HammingDistance(DNA1, DNA2) == 1:
                BaseRoundDNAsForest.add_edge(DNA1,DNA2,
                                                MutationsNumber = 1)
    # adds edges to the BaseRoundDNAsForest so that it can be made into a directed graph
    
    BaseRoundDNAsTrees = list(nx.connected_component_subgraphs(BaseRoundDNAsForest, copy = True))
    
    DNAsClonesSummaryReportFileNameCSV =  str(today) + 'DNAsMutantsSummary' + DNAsClonesSummaryReportFileName + '.csv'
    DNAsClonesSummaryReportFile = open(DNAsClonesSummaryReportFileNameCSV, 'w')
    
    Positions = {}
    X_0_Coordinate = 1
    Y_0_Coordinate = 0
    Y_X0_Coordinate = 0
    
    for Tree in BaseRoundDNAsTrees:
                        
        RootDNA = max(nx.get_node_attributes(Tree, 'Occurence'),
                        key = lambda(DNA): Tree[DNA])
        
        TreeDNAs = {}
        TreeDNAs[RootDNA] = [0, '', 0, Tree.node[RootDNA]['Occurence'], Tree.node[RootDNA]['FirstAppearance']]
        TreeDNAsList = Tree.nodes()
        TreeDNAsList.remove(RootDNA)
        
        for DNA in TreeDNAsList:
            DNAPredecessor = nx.shortest_path(Tree, source = DNA, target = RootDNA, weight = None)[1]
            PredecessorOccurence = Tree.node[DNAPredecessor]['Occurence']
            DNAOccurence = Tree.node[DNA]['Occurence']
            
            TreeDNAs[DNA] = [DNAPredecessor, PredecessorOccurence, DNAOccurence]
        
        DNAs_BY_DistanceToTheRoot = {}
        for DNA in Tree.nodes():
            DistanceToTheRoot = nx.shortest_path_length(Tree, source = DNA, target = RootDNA, weight = None)
            if DistanceToTheRoot not in DNAs_BY_DistanceToTheRoot:
                DNAs_BY_DistanceToTheRoot[DistanceToTheRoot] = [DNA]
            else:
                DNAs_BY_DistanceToTheRoot[DistanceToTheRoot] += [DNA]
                
        MaxDNAMutantsNumber = max(map(lambda(k): len(DNAs_BY_DistanceToTheRoot[k]), DNAs_BY_DistanceToTheRoot))
        
        SortedDNAs_BY_DistanceToTheRoot = {}
        
        for DistanceToTheRoot in DNAs_BY_DistanceToTheRoot:
            
            EquidistantDNAs = DNAs_BY_DistanceToTheRoot[DistanceToTheRoot]
            
            EquidistantDNAs = sorted(EquidistantDNAs, key = lambda (DNA): (TreeDNAs[DNA][2]), reverse = True)
#            EquidistantDNAs = sorted(EquidistantDNAs, key = lambda (DNA): (TreeDNAs[DNA][1]), reverse = True)
            EquidistantDNAs = sorted(EquidistantDNAs, key = lambda (DNA): (TreeDNAs[DNA][0]), reverse = False)
            
            AdditionalElements = MaxDNAMutantsNumber - len(EquidistantDNAs)
            SortedDNAs_BY_DistanceToTheRoot[DistanceToTheRoot] = EquidistantDNAs + AdditionalElements * ['']
            
            if len(Tree.nodes()) > 1:
                for DNA in EquidistantDNAs:
                    XCoordinate = X_0_Coordinate + DistanceToTheRoot
                    YCoordinate = Y_0_Coordinate - EquidistantDNAs.index(DNA)
                    Positions[DNA] = (XCoordinate, YCoordinate)
            elif len(Tree.nodes()) == 1:
                for DNA in EquidistantDNAs:
                    XCoordinate = 0
                    YCoordinate = Y_X0_Coordinate
                    Positions[DNA] = (XCoordinate, YCoordinate)
        
                            
        if len(Tree.nodes()) > 1:
            X_0_Coordinate += max(DNAs_BY_DistanceToTheRoot.keys()) + 1
        
        if len(Tree.nodes()) == 1:
            Y_X0_Coordinate -= 1
            
            
#        print SortedDNAs_BY_DistanceToTheRoot
        

#-------------------------------------------------------------------------------

        if len(Tree.nodes()) > 1:
            for DistanceToTheRoot in SortedDNAs_BY_DistanceToTheRoot:
                DNAsClonesSummaryReportFile.write(str(DistanceToTheRoot) + ' mutations' + ',' + 'frequency' + ',' + 'clones #' + ',' )
            DNAsClonesSummaryReportFile.write('\n')
        
        for i in range(MaxDNAMutantsNumber):
            for MutationsNumber in SortedDNAs_BY_DistanceToTheRoot:                        
                DNA = SortedDNAs_BY_DistanceToTheRoot[MutationsNumber][i]
                
                if DNA != '':
                    FormatedDNA = HammingDistanceBasedFormating(RootDNA, DNA)
                    ClonesNumber = str(len(Tree.neighbors(DNA)))
                    DNAFraction = ('{:.2%}'.format(float((DNAs_BY_Round[BaseRoundIndex].get(DNA, 0)))/float(TotalDNAs_BY_Round[BaseRoundIndex])))
                else:
                    FormatedDNA = ''
                    ClonesNumber = ''
                    DNAFraction = ''
                                                        
                DNAsClonesSummaryReportFile.write(FormatedDNA + ',' +
                            DNAFraction + ',' +
                            ClonesNumber + ',')
            DNAsClonesSummaryReportFile.write('\n')
        DNAsClonesSummaryReportFile.write('\n')
                    
    DNAsClonesSummaryReportFile.close()
        
    BaseRoundDNAsGraph = nx.Graph()    
    BaseRoundDNAsGraph.add_nodes_from(BaseRoundSortedDNAs)
     
    Sizes = []
    for DNA in BaseRoundDNAsGraph.nodes():
        Sizes.append(5 * math.log(DNAs_BY_Round[BaseRoundIndex][DNA], 10) + 2)
#    Colours = []
#    for DNA in BaseRoundDNAsGraph.nodes():
#        Colours.append(min(DNAsAppearances[DNA]))
    
    Colours = []
    for DNA in BaseRoundDNAsGraph.nodes():
        Peptide = Translation(DNA)
        Colours.append(BaseRoundSortedPeptides.index(Peptide))
    
#    print(BaseRoundDNAsGraph.nodes())
    i = 0
    for DNA1 in BaseRoundSortedDNAs:
        for DNA2 in BaseRoundSortedDNAs:
            if HammingDistance(DNA1, DNA2) == 1:
                BaseRoundDNAsGraph.add_edge(DNA1,DNA2)
        i += 1

    XMin = min(map(lambda(DNA): Positions[DNA][0], Positions)) - 1
    XMax = max(map(lambda(DNA): Positions[DNA][0], Positions)) + 1
    YMin = min(map(lambda(DNA): Positions[DNA][1], Positions)) - 5
    YMax = max(map(lambda(DNA): Positions[DNA][1], Positions)) + 5
    
    nx.draw_networkx(BaseRoundDNAsGraph,
                    pos = Positions,
                    node_size = Sizes,
                    node_color = Colours,
                    cmap = 'Paired',
                    linewidths = 0.1,
                    width = 0.0,
                    with_labels = False,
                    font_size = 6)

    plt.axis('off')
    plt.axis([XMin, XMax, YMin, YMax])
    
    DNAsMutantsSummaryReportFileNamePNG = str(today) + 'DNAsMutantsSummary' + DNAsClonesSummaryReportFileName + '.png'
    plt.savefig(DNAsMutantsSummaryReportFileNamePNG, dpi = 500)
    
    plt.show()
    plt.close()


#_____________________________RUNNING THE FUNCTIONS_____________________________#

#___DataFolderLocation, BaseSelectionRoundNumber, TopNPeptidesNumber, SummaryFileName___
SelectionSummaryReport('/Users/NikitaLoik/SugaLabDataAnalysis/Chris1Jul2016/160211 LSD1 D generations', 7, 24, 'LSD1 NNK D Chris')

#___DataFolderLocation, BaseSelectionRoundNumber, TopNPeptidesNumber, SummaryFileName___
DNAMutantsAnalysis('/Users/NikitaLoik/SugaLabDataAnalysis/Chris1Jul2016/160211 LSD1 D generations', 7, 24, 'LSD1 NNK D Chris')

#___DataFolderLocation, BaseSelectionRoundNumber, TopNPeptidesNumber, SummaryFileName___
PeptidesMutantsSummaryReport('/Users/NikitaLoik/SugaLabDataAnalysis/Chris1Jul2016/160211 LSD1 D generations', 7, 'LSD1 NNK D Chris')

#___DataFolderLocation, BaseSelectionRoundNumber, SummaryFileName___
DNAsMutantsSummaryReport('/Users/NikitaLoik/SugaLabDataAnalysis/Chris1Jul2016/160211 LSD1 D generations', 7, 'LSD1 NNK D Chris')
