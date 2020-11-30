def PeptideMasses(PeptideMassesSummaryFileName, PeptidesListFileLocation):
    # 'amino acid', monoisotopic residue mass, average residue mass
    AAResidueMassData = {'A':('Ala', 71.037114, 71.0779),
                'C':('Cys', 103.009185, 103.1429),
                'D':('Asp', 115.026943, 115.0874),
                'E':('Glu', 129.042593, 129.114),
                'F':('Phe', 147.068414, 147.1739),
                'G':('Gly', 57.021464, 57.0513),
                'H':('His', 137.058912, 137.1393),
                'I':('Ile', 113.084064, 113.1576),
                'K':('Lys', 128.094963, 128.1723),
                'L':('Leu',  113.084064, 113.1576),
                'M':('Met', 131.040485, 131.1961),
                'N':('Asn', 114.042927, 114.1026),
                'P':('Pro', 97.052764, 97.1152),
                'Q':('Gln', 128.058578, 128.1292),
                'R':('Arg', 156.101111, 156.1857),
                'S':('Ser', 87.032028, 87.0773),
                'T':('Thr', 101.047679, 101.1039),
                'V':('Val', 99.068414, 99.1311),
                'W':('Trp', 186.079313, 186.2099),
                'Y':('Tyr', 163.06332, 163.1733),
                'y':('D-Tyr', 163.06332, 163.1733),
                'X':('HONH-Glu',144.05349, 144.1300),
                'Z':('HONH-ASub',186.10044, 186.2110)
                }
            
    PeptidesListFile = open(PeptidesListFileLocation, 'r')
    Lines = PeptidesListFile.readlines()
    PeptidesListFile.close    
# print Lines

    PeptideMassesFile = open(PeptideMassesSummaryFileName, 'w')
    
    PeptideMassesFile.write('#' + ',' +
                            'peptide' + ',' +
                            'MI linear (Da)' + ',' +
                            'A linear (Da)' + ',' +
                            'MI cyclic (Da)' + ',' +
                            'A cyclic (Da)' + ',' +
                            'Ext 280nm (1/M*cm)''\n')
    
    PeptideNumber = 0
        
    for Line in Lines:
        Peptide = Line.strip('\n')
        
        PeptideNumber += 1
        LinearPeptideMIM = 17.02655
        LinearPeptideAM = 17.0310
        
        CysExt = 110.0
        TrpExt = 5630.0
        TyrExt = 1215.0
        #uses the data from this paper http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2143013/pdf/8563639.pdf
        #alternatively these data can be used http://www.rpgroup.caltech.edu/courses/aph162/2007/pdf/articles/Gill1989.pdf
        
        CysCount = 0
        TrpCount = 0
        TyrCount = 0
        
        for i in range(len(Peptide)):
            LinearPeptideMIM += AAResidueMassData[Peptide[i]][1]
            LinearPeptideAM += AAResidueMassData[Peptide[i]][2]
            if Peptide[i] == 'C':
                CysCount += 1
            elif Peptide[i] == 'W':
                TrpCount += 1
            elif Peptide[i] == 'Y' or Peptide[i] == 'y':
                TyrCount += 1
                
        CyclicPeptideMIM = LinearPeptideMIM + 42.01056
        #print MICyclic
        CyclicPeptideAM = LinearPeptideAM + 42.0370
        #print ACyclic
        PeptideExt = CysCount * CysExt + TrpCount *TrpExt + TyrCount * TyrExt
        #
        
        PeptideMassesFile.write(str(PeptideNumber) + ',' +
                                Peptide + ',' +
                                '{:.2f}'.format(LinearPeptideMIM) + ',' +
                                '{:.2f}'.format(LinearPeptideAM) + ',' +
                                '{:.2f}'.format(CyclicPeptideMIM) + ',' +
                                '{:.2f}'.format(CyclicPeptideAM) + ',' +
                                '{:.0f}'.format(PeptideExt) + '\n')
                                
    PeptideMassesFile.close

#_____________________________RUNNING THE FUNCTION_____________________________#
#___PeptideMassesSummaryFileName, PeptidesListFileLocation___
PeptideMasses('CloneMasses Result Test.csv', '/Users/NikitaLoik/CloneSynthesis.txt')