def AAsInPeptideListCount(PeptidesListFileLocation):
    PeptidesListFile = open(PeptidesListFileLocation, 'r')
    Lines = PeptidesListFile.readlines()
    PeptidesListFile.close
    
    AminoAcidsCount = {'A':0, 'C':0, 'D':0, 'E':0,
                        'F':0, 'G':0, 'H':0, 'I':0,
                        'K':0, 'L':0, 'M':0, 'N':0,
                        'P':0, 'Q':0, 'R':0, 'S':0,
                        'T':0, 'V':0, 'W':0, 'Y':0,
                        'y':0, 'X':0, 'Z':0
                    }
        
    # populate the dictionary, so that Peptides are the keys and
    for Line in Lines:
        Line = Line.strip('\n')
        for i in range(len(Line)):
           AminoAcidsCount[Line[i]] += 1
    return AminoAcidsCount
    
def AAQuantitiesForSYRO(AAQuantitiesForSYROFileName, PeptidesListFileLocation):
    AAData = {'A':('Ala','Fmoc-Ala-OH H2O', 311.34),
                'C':('Cys','Fmoc-Cys(Trt)-OH', 585.72),
                'D':('Asp','Fmoc-Asp(OtBu)-OH',411.46),
                'E':('Glu','Fmoc-Glu(OtBu)-OH',425.49),
                'F':('Phe','Fmoc-Phe-OH',387.40),
                'G':('Gly','Fmoc-Gly-OH',297.31),
                'H':('His','Fmoc-His(Trt)-OH',619.72),
                'I':('Ile','Fmoc-Ile-OH',353.42),
                'K':('Lys','Fmoc-Lys(Boc)-OH',468.55),
                'L':('Leu','Fmoc-Leu-OH',353.42),
                'M':('Met','Fmoc-Met-OH',371.45),
                'N':('Asn','Fmoc-Asn(Trt)-OH',596.68),
                'P':('Pro','Fmoc-Pro-OH',337.38),
                'Q':('Gln','Fmoc-Gln(Trt)-OH',610.72),
                'R':('Arg','Fmoc-Arg(Pbf)-OH',648.77),
                'S':('Ser','Fmoc-Ser(tBu)-OH',383.45),
                'T':('Thr','Fmoc-Thr(tBu)-OH',397.48),
                'V':('Val','Fmoc-Val-OH',339.39),
                'W':('Trp','Fmoc-Trp(Boc)-OH',526.59),
                'Y':('Tyr','Fmoc-Tyr(tBu)-OH',459.54),
                'y':('D-Tyr','Fmoc-D-Tyr(tBu)-OH',459.54),
                'X':('HONH-Glu','Fmoc-(tBu)ONH-Glu-OH',440.50),
                'Z':('HONH-ASub','Fmoc-(tBu)ONH-ASub-OH',482.50)
                }
        
    AAQuantitiesForSYRO = open(AAQuantitiesForSYROFileName, 'w')
    
    AAQuantitiesForSYRO.write('AA' + ',' +
                            'Name' + ',' +
                            'Formula' + ',' +
                            '#' + ',' +
                            'MW(g/mol)' + ',' +
                            'Q(mol)' + ',' +
                            'Q(g)' + ',' +
                            'V (mL)' + '\n')
    
    AAList = AAsInPeptideListCount(PeptidesListFileLocation)
    
    TotalNumberOfSteps = sum(AAList.values())
    
    for AminoAcid in AAList:
        AAName = AAData[AminoAcid][0]
        AAFormula = AAData[AminoAcid][1]
        AACount = AAList[AminoAcid]
        
        if AACount > 0:
            AADilutionVolume = 1.1 * (2 * (0.0006 + (AACount - 1) * 0.0003))
        else:
            AADilutionVolume = 0
                    
        AAMolecularWeight = AAData[AminoAcid][2]
        AAMoleQuantity = AADilutionVolume * 0.5
        AAGrQuantity = AAMoleQuantity * AAMolecularWeight
        
        AAQuantitiesForSYRO.write(AminoAcid + ',' +
                                AAName + ',' +
                                AAFormula + ',' +
                                str(AACount) + ',' +
                                '{:.2f}'.format(AAMolecularWeight) + ',' +
                                '{:.3f}'.format(AAMoleQuantity) + ',' +
                                '{:.3f}'.format(AAGrQuantity) + ',' +
                                '{:.3f}'.format(AADilutionVolume * 1000) + ',' + '\n')
                                
    MWofHBTU = 379.25
    MWofHOBt = 171.134
    MWofDIPEA = 129.25
    DofDIPEA = 0.742
                                            
    HBTUinDMFvolume = 2.2 * TotalNumberOfSteps * 0.000345
    QofHBTU = 0.43 * HBTUinDMFvolume * MWofHBTU
    QofHOBt = 0.43 * HBTUinDMFvolume * MWofHOBt
        
    DIPEAinNMPvolume = 2.2 * TotalNumberOfSteps * 0.000150
    QofDIPEA = 2.2 * DIPEAinNMPvolume * MWofDIPEA
    VofDIPEA = QofDIPEA/DofDIPEA
    
    VofPiperidineInDMF = 2.2 * TotalNumberOfSteps * 0.000900
    
    AAQuantitiesForSYRO.write('HBTU quantity (g)' + ',' + '{:.2f}'.format(QofHBTU)  + '\n' +
                            'HOBt.H2O quantity (g)' + ',' + '{:.2f}'.format(QofHOBt)  + '\n' +
                            'HBTU & HOBt.H2O in DMF (mL)' + ',' + '{:.2f}'.format(1000 * HBTUinDMFvolume)  + '\n' +
                            'DIPEA quantity (g)' + ',' + '{:.2f}'.format(QofDIPEA)  + '\n' +
                            'DIPEA volume (mL)' + ',' + '{:.2f}'.format(VofDIPEA)  + '\n' +
                            'DIPEA in NMP (mL)' + ',' + '{:.2f}'.format(1000 * DIPEAinNMPvolume)  + '\n' +
                            '40% piperidine in DMF (mL)' + ',' + '{:.2f}'.format(1000 * VofPiperidineInDMF)  + '\n')
       
    AAQuantitiesForSYRO.close

#_____________________________RUNNING THE FUNCTION_____________________________#
#___AAQuantitiesForSYROFileName, PeptidesListFileLocation___
AAQuantitiesForSYRO('CloneSynthesis Test.csv', '/Volumes/NIKITA 2GB/CloneSynthesis.txt')
