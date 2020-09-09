QUALITY_SCORE_STRING = '''
    !"#$%&'()*+,-./0123456789:;<=>?@
    ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`
    abcdefghijklmnopqrstuvwxyz{|}~
    '''

QUALITY_SCORE = 29  # threshold quality_score must be between 0 and 93

TRANSLATION_CODE = {
    'AAA':'K','AAC':'N','AAG':'K','AAU':'N',
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
    'UUA':'L','UUC':'F','UUG':'L','UUU':'F'
    }


TRANSCRIPTION_CODE = {
    'A':'A','C':'C','G':'G','T':'U','U':'T'
    }

CDNA_MIN_LENGTH = 24
CDNA_MAX_LENGTH = 240

START_SEQUENCE = 'TAATACGACTCACTATAGGGTTAACTTTAAGAAGGAGATATACATATG'    # NNK - T7g10M.F48 
STOP_SEQUENCE = 'TGCGGCAGCGGCAGCGGCAGCTAGGACGGGGGGCGGAAA' #NNK - CGS3an13.R39
# START_SEQUENCE = 'TAATACGACTCACTATAGGGTTGAACTTTAAGTAGGAGATATATCCATG'   #NNU - T7-CH-F49
# STOP_SEQUENCE = 'TGTGGGTCTGGGTCTGGGTCTTAGGACGGGGGGCGGAAA'  #NNU - CGS3-CH-R39
# START_SEQUENCE = 'ATG' # Met codon
# STOP_SEQUENCE = 'TGCGGCAGC'# Akane seams to have trimmed siquences
# START_SEQUENCE = 'TAGGGTTAACTTTAAGAAGGAGATATACATATG'# Oxford, Akane and Tom
# STOP_SEQUENCE = 'TGCGGC'# Oxford, Akane and Tom
# STOP_SEQUENCE = 'TAG' # amber stop codon


TOP_SORTED_PEPTIDES = [
    'VWDPRTFYLSRI', 'WDANTIFIKRV', 'WNPRTIFIKRA', 'VWDPRTFYLSRT',
    'IWDTGTFYLSRT', 'WWNTRSFYLSRI', 'FWDPRTFYLSRI', 'VWDPSTFYLSRI',
    'KWDTRTFYLSRY', 'KWDTRTFYLSRI', 'IWDPRTFYLSRI', 'IWDTGTFYLSRI',
    'VWDPRTFYLSRM', 'AWDPRTFYLSRI', 'VWDSRTFYLSRI', 'VWDPGTFYLSRI',
    'VWDPRTFYMSRI', 'VWDPRTFYLSRS', 'VWDPRTFYLSRV', 'WNPRTIFIKRV',
    'VRDPRTFYLSRI', 'VWDPKTFYLSRI', 'VWDPRTFYLSRN', 'FRFPFYIQRR'
    ]

TOP_PEPTIDES_KDS = {
    'VWDPRTFYLSRI' : '3',
    'WDANTIFIKRV' : '4',
    'WNPRTIFIKRA' : '>1000',
    'VWDPRTFYLSRT' : '3',
    'IWDTGTFYLSRT' : '7',
    'WWNTRSFYLSRI' : '12',
    'FWDPRTFYLSRI' : '4',
    'VWDPSTFYLSRI' : '3',
    'KWDTRTFYLSRY' : '5',
    'KWDTRTFYLSRI' : '6',
    'IWDPRTFYLSRI' : '1',
    'VWDPRTFYLSRM' : '4',
    'IWDTGTFYLSRI' : '>1000',
    'VWDPGTFYLSRI' :  '<1',
    'VWDSRTFYLSRI' : '3',
    'AWDPRTFYLSRI': '6',
    'VWDPRTFYLSRS' : '6',
    'VWDPRTFYMSRI' : '1',
    'VWDPRTFYLSRV' : '3',
    'WNPRTIFIKRV' : '>1000',
    'VRDPRTFYLSRI' : '>1000',
    'VWDPRTFYLSRN' : '>1000',
    'VWDPKTFYLSRI' : '14',
    'FRFPFYIQRR' : '>1000'
    }


def display_summaryReport(
    data_directory_path, base_cycle, n_top_peptides, file_name):
    
    today = TodaysDate() 
    
    display_summaryFileNameCSV = str(today) + 'display_summary' + file_name + '.csv'
    display_summaryReportFile = open(display_summaryFileNameCSV, 'w')
    
    display_summary = Completedisplay_summary(data_directory_path)
    SortedRoundsList = sorted(display_summary.keys())
    
    peptides_BY_Round = peptidesOccurrences_BY_Round(data_directory_path)
    Totalpeptides_BY_Round = TotalReads_BY_Round(data_directory_path)
    
    BaseRoundSortedpeptides = BaseRoundSortedpeptidesList(data_directory_path, base_cycle)
    #for i in range(len(BaseRoundSortedpeptides)):
     #   print ('>seq' + str(i + 1) + '\n' + BaseRoundSortedpeptides[i])
    #BaseRoundTopSortedpeptides = BaseRoundSortedpeptides[0 : (n_top_peptides)]
    BaseRoundTopSortedpeptides = ['VWDPRTFYLSRI', 'WDANTIFIKRV', 'WNPRTIFIKRA', 'VWDPRTFYLSRT',
                                'IWDTGTFYLSRT', 'WWNTRSFYLSRI', 'FWDPRTFYLSRI', 'VWDPSTFYLSRI',
                                'KWDTRTFYLSRY', 'KWDTRTFYLSRI', 'IWDPRTFYLSRI', 'IWDTGTFYLSRI',
                                'VWDPRTFYLSRM', 'AWDPRTFYLSRI', 'VWDSRTFYLSRI', 'VWDPGTFYLSRI',
                                'VWDPRTFYMSRI', 'VWDPRTFYLSRS', 'VWDPRTFYLSRV', 'WNPRTIFIKRV',
                                'VRDPRTFYLSRI', 'VWDPKTFYLSRI', 'VWDPRTFYLSRN', 'FRFPFYIQRR'
                                 ]
    BaseRoundpeptidesRank = peptidesRank_IN_BaseRound(data_directory_path, base_cycle)
    #print (BaseRoundSortedpeptides)
    
    Top24peptidesKDs = {'VWDPRTFYLSRI' : '3', 'WDANTIFIKRV' : '4', 'WNPRTIFIKRA' : '>1000', 'VWDPRTFYLSRT' : '3',
                        'IWDTGTFYLSRT' : '7', 'WWNTRSFYLSRI' : '12', 'FWDPRTFYLSRI' : '4', 'VWDPSTFYLSRI' : '3',
                        'KWDTRTFYLSRY' : '5', 'KWDTRTFYLSRI' : '6', 'IWDPRTFYLSRI' : '1', 'VWDPRTFYLSRM' : '4',
                        'IWDTGTFYLSRI' : '>1000', 'VWDPGTFYLSRI' :  '<1', 'VWDSRTFYLSRI' : '3', 'AWDPRTFYLSRI': '6',
                        'VWDPRTFYLSRS' : '6', 'VWDPRTFYMSRI' : '1', 'VWDPRTFYLSRV' : '3', 'WNPRTIFIKRV' : '>1000',
                        'VRDPRTFYLSRI' : '>1000', 'VWDPRTFYLSRN' : '>1000', 'VWDPKTFYLSRI' : '14', 'FRFPFYIQRR' : '>1000'
                       }
        
    display_summaryReportFile.write('peptide sequence' + ',' +
                                     'rank (#)' + ',' +
                                     'cDNA mutants' + ',')
    for Round in SortedRoundsList:
        display_summaryReportFile.write('C' +
                                         str(Round) +
                                         ' count (#) [frequency(%)]' + ',')
    display_summaryReportFile.write('\n')
    
    for peptide in BaseRoundTopSortedpeptides:
    #for peptide in Top24peptidesKDs:
        BaseRoundpeptideFraction = float((peptides_BY_Round[Round].get(peptide, 0)))/float(Totalpeptides_BY_Round[base_cycle])
        peptideRank = BaseRoundpeptidesRank[peptide]
        Formatedpeptide = HammingDistanceBasedFormating(BaseRoundTopSortedpeptides[0], peptide)
        peptidecDNAMutants = len(display_summary[base_cycle][peptide])
        display_summaryReportFile.write(Formatedpeptide + ',' +
                                         str(peptideRank) + ',' +
                                         str(peptidecDNAMutants) + ',')

            
        for Round in SortedRoundsList:
            peptideFraction = float((peptides_BY_Round[Round].get(peptide, 0)))/float(Totalpeptides_BY_Round[Round])

            BaseFraction = peptideFraction
            
            display_summaryReportFile.write(str(peptides_BY_Round[Round].get(peptide, 0)) +
                                             ' [' + '{:.1%}'.format(peptideFraction) + ']' + ',')
            
            BaseFraction = peptideFraction
        display_summaryReportFile.write('\n')
        
    display_summaryReportFile.write('total count (#)' + ',' + ',')
    for Round in SortedRoundsList:
        display_summaryReportFile.write(str(Totalpeptides_BY_Round[Round]) + ',')
    display_summaryReportFile.write('\n\n\n')
            
    display_summaryReportFile.close()
    
#-------------------------------------------------------------------------------
   
    # Create a figure of size 8x6 inches, 500 dots per inch
    plt.figure(
        figsize = (8, 6),
        dpi = 500)
    # Create 'ggplot' style
    plt.style.use('fivethirtyeight')
    # Create a new subplot from a grid of 1x1
    Graph = plt.subplot(1, 1, 1)
    
    Xs = []
    Ys = []
    
    Rank = 1
    peptideFractionInFinalRound = 0
    
    # Map colours onto lines
    
    cNorm  = matplotlib.colors.Normalize(vmin = 0,
                                         vmax = n_top_peptides - 1)
    scalarMap = matplotlib.cm.ScalarMappable(norm = cNorm,
                                             cmap = 'Paired')
    
    peptideLabels = []
    
    for peptide in BaseRoundTopSortedpeptides:
    #for peptide in Top24peptidesKDs:
        peptidesFractions_BY_Round = []
        for Round in SortedRoundsList:
            peptidesFractions_BY_Round += [float((peptides_BY_Round[Round].get(peptide, 0)))/float(Totalpeptides_BY_Round[Round])]
        
        x = SortedRoundsList
        y = peptidesFractions_BY_Round
        Xs += x
        Ys += y
        
        peptideColour = scalarMap.to_rgba(BaseRoundTopSortedpeptides.index(peptide))
        
        peptideRank = str(BaseRoundpeptidesRank[peptide])
#         peptideColour = scalarMap.to_rgba(peptideRank)
        peptideKD = Top24peptidesKDs[peptide]
        Formatedpeptide = HammingDistanceBasedFormating(BaseRoundTopSortedpeptides[0], peptide)
        
        peptideLabel =  Formatedpeptide + ' (' + peptideRank + ', ' + peptideKD +' nM)'
        
        #Set peptideLabel
        peptideLabels += [peptideLabel]
        
        plt.plot(x, y,
                 'o-',
                 c = peptideColour,
                 lw = 2.0,
                 ms = 4.0,
                 mew = 0.1,
                 mec = '#191919')

    
    XMin = min(Xs) - 0.05*(max(Xs) - min(Xs))
    XMax = max(Xs) + 0.05*(max(Xs) - min(Xs))
    YMin = min(Ys) - 0.05*(max(Ys) - min(Ys))
    YMax = max(Ys) + 0.05*(max(Ys) - min(Ys))
    
    plt.axis([XMin, XMax, YMin, YMax])
    
    plt.xticks(fontsize = 10)
    plt.yticks(fontsize = 10)
    
    plt.xlabel('Selection Cycle (#)',
               fontsize = 10)
    plt.ylabel('peptide Fraction (%)',
               fontsize = 10)
    
    legend = plt.legend(peptideLabels,
                        title = 'cyclic-peptide random region',
                        loc = 'upper center',
                        bbox_to_anchor = (0.5, -0.10),
                        fancybox = True,
                        shadow = False,
                        fontsize = 10,
                        ncol = 3)
    
    Graph.get_legend().get_title().set_size('small')
    
    display_summaryFileNamePNG = str(today) + 'display_summary' + file_name + '.png'
    
    plt.savefig(display_summaryFileNamePNG,
                bbox_extra_artists = [legend],
                bbox_inches = 'tight',
                dpi = 300)
    plt.show()
    plt.close()