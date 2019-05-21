import os
import datetime
#rcParams['font.sans-serif'] = ['Tahoma']

import math
import numpy as np
import networkx as nx

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'monospace'

import globals

def todays_date():
    '''
    returns todays date in format DDMMYYYY
    '''
    return datetime.date.today().strftime('%d%b%Y')


def dna_coding_sequence (
        dna_sequence: str,
        quality_sequence: str,
        start_sequence: str,
        stop_sequence: str,
        cdna_min_length: int,
        cdna_max_length: int,
        quality_score,
        ):
    '''
    validates the dna coding sequence based on
    (1) its lenght AND
    (2) quality score of each nucleotide;
    uses ONLY ONE stop sequence;
    returns ONLY ONE coding sequence
    '''
    quality_score_string = globals.QUALITY_SCORE_STRING
    threshold_quality_string = quality_score_string[quality_score:]
            
    start_index = dna_sequence.find(start_sequence) + len(start_sequence)
    stop_index = dna_sequence.rfind(stop_sequence)
    coding_sequence =  dna_sequence[start_index:stop_index]
    if ((cdna_min_length <= len(coding_sequence))
    and (len(coding_sequence) <= cdna_max_length)
    and (len(coding_sequence)%3 == 0)):
        for character in quality_sequence[start_index:stop_index]:
            if character not in threshold_quality_string:
                return None
        return str(coding_sequence)


# 

def translation(
        coding_sequence: str,
        ):
    '''
    translates coding sequence into peptide;
    (peptide can be encoded as dna or as RNA)
    stop codons are as follows:
    UAA (ochre) — #
    UAG (amber) — *
    UGA (opal) — &
    '''
    translation_code = {
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
    transcription_code = {
        'A':'A','C':'C','G':'G','T':'U','U':'T'
        }
    
    # Convert dna to RNA
    rna_sequence = ''
    for nucleotide in coding_sequence:
        rna_sequence += transcription_code.get(nucleotide,'X')
        
    #rna_sequence-Test Print
    #print (rna_sequence)
        
    peptide = ''
    while len(rna_sequence) != 0:
        peptide += translation_code.get(
            rna_sequence[0:3],
            "Do not fuck with me!")
        rna_sequence = rna_sequence[3:]
    return peptide


# Define single_selection_cycle_summary function, which returns the occurence of coding_sequences, grouped by peptide
# single_selection_cycle_summary = {peptideY:    {coding_sequence_YZ:    occurence_YZ}}

def single_selection_cycle_summary(
        fastq_file_path: str,
        start_sequence: str,
        stop_sequence: str,
        cdna_min_length: int,
        cdna_max_length: int,
        quality_score: int,
        ):
    
    with open(fastq_file_path, 'r') as raw_data_file:
        lines = raw_data_file.readlines()

    # to store the results from a single cycle of selection,
    # create an empty single_selection_cycle_summary dictionary:
    # {peptideY:    {coding_sequence_YZ:    occurence_YZ}}
    single_selection_cycle_summary = {}

    # go through the NGS-data file (.fastq) line by line
    # with the results from a single cycle of selection,
    # populate single_selection_cycle_summary 
    for i, line in enumerate(lines):
        # Check whether the line contains a valid sequence
        if start_sequence in line and stop_sequence in line:
            coding_sequence = dna_coding_sequence(
                line,
                lines[i + 2],
                start_sequence,
                stop_sequence,
                cdna_min_length,
                cdna_max_length,
                quality_score)                      
            if coding_sequence != None:
                # translate dna or RNA sequence into peptide
                peptide_sequence = translation(coding_sequence)
                # add sequence to a single_selection_cycle_summary
                if peptide_sequence not in single_selection_cycle_summary:
                    single_selection_cycle_summary[str(peptide_sequence)] = {str(coding_sequence): 1}
                else:
                    if coding_sequence not in single_selection_cycle_summary[str(peptide_sequence)]:
                        single_selection_cycle_summary[str(peptide_sequence)][str(coding_sequence)] = 1
                    else:
                        single_selection_cycle_summary[str(peptide_sequence)][str(coding_sequence)] += 1

    return single_selection_cycle_summary


def hamming_distance(
        sequence_1: str,
        sequence_2: str,
        ):
    '''
    returns the number of mismatches between two sequences
    (i.e. hamming distance)
    '''
    
    if len(sequence_1) < len(sequence_2):
        sequence_1 = sequence_1 + (len(sequence_2) - len(sequence_1)) * '%'
    elif len(sequence_1) > len(sequence_2):
        sequence_2 = sequence_2 + (len(sequence_1) - len(sequence_2)) * '%'
    
    hamming_distance = 0
    for i in range(len(sequence_1)):
        if sequence_1[i] == sequence_2[i]:
            hamming_distance = hamming_distance
        else:
            hamming_distance = hamming_distance + 1
            
    return hamming_distance


def hamming_distance_based_formating(
        sequence_1: str,
        sequence_2: str,
        ):
    '''
    # returns formated sequence,
    # such that mismatches are capitalised,
    # and the rest of the formated sequence is lowercase 
    '''
    
    if len(sequence_1) < len(sequence_2):
        sequence_1 = sequence_1 + (len(sequence_2) - len(sequence_1)) * '-'
    elif len(sequence_1) > len(sequence_2):
        sequence_2 = sequence_2 + (len(sequence_1) - len(sequence_2)) * '-'
    
    hamming_distance = 0
    formated_sequence_2 = ''
    for i in range(len(sequence_1)):
        if sequence_1[i] == sequence_2[i]:
            formated_sequence_2 += sequence_2[i].lower()
            hamming_distance = hamming_distance
        else:
            formated_sequence_2 += sequence_2[i]
            hamming_distance = hamming_distance + 1            
    return formated_sequence_2



def complete_selection_summary(
        data_directory_path: str,
        start_sequence: str,
        stop_sequence: str,
        cdna_min_length: int,
        cdna_max_length: int,
        quality_score: int,
        ):
    '''
    returns the occurence of coding sequences
    grouped by peptide, and grouped by selection cycle:
    {Selectioncycle_X:    {peptideXY:    {Codingdna_XYZ:    occurence_XYZ}}}
    '''

    # Create an empty selection_summary dictionary to store the results from all the cycles of selection
    complete_selection_summary = {}
    
    for file in os.listdir(data_directory_path):
        file_path = os.path.join(data_directory_path, file)
        
        # A. get the cycle number from the file name
        # (file name should have two digit number before full stop — '00.')    
        if file.endswith('.fastq'):  # this condition is necessary; without it some shit appears in the beginning of the file list
            cycle_number_first_digit = file[file.find('.')-2]
            cycle_number_second_digit = file[file.find('.')-1]
            if cycle_number_first_digit == '0':
                cycle_number = int(cycle_number_second_digit)
                # cycle_number_Test Print
                #print cycle_number
            elif cycle_number_first_digit != '0':
                cycle_number = int(file[file.find('.')-2 : file.find('.')])
                # cycle_number_Test Print
                #print cycle_number
            
            # B. get single cycle results
            selection_cycle_summary = single_selection_cycle_summary(
                file_path,
                start_sequence,
                stop_sequence,
                cdna_min_length,
                cdna_max_length,
                quality_score)
             
            # C populate complete selection summary
            complete_selection_summary[cycle_number] = selection_cycle_summary
            
            # ConcatenatedResultsList-Test Print
            #print ConcatenatedResultsList
            
    return complete_selection_summary


def peptides_occurences_by_cycle(
        data_directory_path: str,
        cdna_min_length: int,
        cdna_max_length: int,
        quality_score: int,
        ):
    '''
    returns the occurences of peptides groupped by cycle:
    {cycle_X:    {peptideXY:    occurence_XY}}
    '''
    selection_summary = complete_selection_summary(
        data_directory_path,
        start_sequence,
        stop_sequence,
        cdna_min_length,
        cdna_max_length,
        quality_score)
    
    peptides_occurences_by_cycle = {}
    for cycle in selection_summary:
        peptides_occurences_in_cycle = {}
        for peptide in selection_summary[cycle]:
            peptides_occurences_in_cycle[peptide] = sum(
                selection_summary[cycle][peptide].values())
        peptides_occurences_by_cycle[cycle] = peptides_occurences_in_cycle
        
    return peptides_occurences_by_cycle


def dnas_occurences_by_cycle(
        data_directory_path: str,
        cdna_min_length: int,
        cdna_max_length: int,
        quality_score: int,
        ):
    '''
    returns the occurences of dnas groupped by cycle:
    {cycle_X:    {dna_XY:    occurence_XY}}
    '''
    selection_summary = complete_selection_summary(
        data_directory_path,
        start_sequence,
        stop_sequence,
        cdna_min_length,
        cdna_max_length,
        quality_score)
    
    dnas_occurences_by_cycle = {}
    for cycle in selection_summary:
        dnas_occurences_in_cycle = {}
        for peptide in selection_summary[cycle]:
            for dna in selection_summary[cycle][peptide]:
                dnas_occurences_in_cycle[dna] = selection_summary[cycle][peptide][dna]
        dnas_occurences_by_cycle[cycle] = dnas_occurences_in_cycle

    return dnas_occurences_by_cycle


def total_reads_by_cycle(
        data_directory_path: str,
        cdna_min_length: int,
        cdna_max_length: int,
        quality_score: int,
        ):
    '''
    returns number of reads by cycle:
    {cycle_X:    TotalReads_X}
    '''
    selection_summary = complete_selection_summary(
        data_directory_path,
        start_sequence,
        stop_sequence,
        cdna_min_length,
        cdna_max_length,
        quality_score)

    peptides_by_cycle = peptides_occurences_by_cycle(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        quality_score)
    
    total_reads_by_cycle = {}
    for cycle in selection_summary:
        total_reads_by_cycle[cycle] = sum(peptides_by_cycle[cycle].values())
        
    return total_reads_by_cycle


def base_cycle_sorted_peptides_list(
        data_directory_path: str,
        base_selection_cycle_number: int,
        ):
    '''
    returns list of peptides in base cycle sorted by their occurence:
    [peptide_1, ..., peptide_N]; occurence(peptide_1) > ... > occurence(peptide_N)
    '''
    peptides_by_cycle = peptides_occurences_by_cycle(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        quality_score)
            
    peptides_occurences_in_base_cycle = peptides_by_cycle[base_selection_cycle_number]
    base_cycle_sorted_peptides_list = sorted(
        peptides_occurences_in_base_cycle,
        key=peptides_occurences_in_base_cycle.get,
        reverse=True)
    return base_cycle_sorted_peptides_list


def base_cycle_sorted_dnas_list(
        data_directory_path: str,
        base_selection_cycle_number: int,
        ):
    '''
    returns list of dnas in base cycle sorted by their occurence:
    [dna_1, ..., dna_n]; occurence(dna_1) > ... > occurence(dna_n)
    '''
    dnas_by_cycle = dnas_occurences_by_cycle(
        data_directory_path: str,
        cdna_min_length: int,
        cdna_max_length: int,
        quality_score: int,
        )
            
    dnas_occurences_in_base_cycle = dnas_by_cycle[base_selection_cycle_number]
    base_cycle_sorted_dnas_list = sorted(
        dnas_occurences_in_base_cycle,
        key=dnas_occurences_in_base_cycle.get,
        reverse=True)
    
    return base_cycle_sorted_dnas_list


def peptides_rank_in_base_cycle(
        data_directory_path: str,
        base_selection_cycle_number: int,
        ):

    peptides_by_cycle = peptides_occurences_by_cycle(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        quality_score)

    base_cycle_sorted_peptides = base_cycle_sorted_peptides_list(
        data_directory_path,
        base_selection_cycle_number)
    
    base_peptide_count = 0
    peptide_rank = 1
    
    peptides_rank_in_base_cycle = {}
    
    for peptide in base_cycle_sorted_peptides:
        peptideCount = peptides_by_cycle[base_selection_cycle_number][peptide]
        if peptideCount < base_peptide_count:
            peptide_rank += 1
        
        peptides_rank_in_base_cycle[peptide] = peptide_rank
        base_peptide_count = peptideCount
        
    return peptides_rank_in_base_cycle


# occurences is a bad word-choice here
def dna_clones_occurences_by_cycle_by_peptide(
        data_directory_path: str,
        cdna_min_length: int,
        cdna_max_length: int,
        quality_score: int,
        ):
    '''
    returns number of clones for each peptide groupped by cycle:
    {Selectioncycle_X:    {peptideXY:    dnaClonesoccurences}}
    '''
    selection_summary = complete_selection_summary(
        data_directory_path,
        start_sequence,
        stop_sequence,
        cdna_min_length,
        cdna_max_length,
        quality_score)
    
    dna_clones_occurences_by_cycle_by_peptide = {}
    for cycle in selection_summary:
        dna_clones_occurences_by_peptide = {}
        for peptide in selection_summary[cycle]:
            dna_clones_occurences_by_peptide[peptide] = len(
                selection_summary[cycle][peptide])
        dna_clones_occurences_by_cycle_by_peptide[cycle] = dna_clones_occurences_by_peptide
        
    return dna_clones_occurences_by_cycle_by_peptide


def peptides_appearances_by_cycle(
        base_cycle_sorted_peptides_list: list,
        peptides_occurences_by_cycle):
    '''
    returns for each peptide in selection a list of cycles in which this peptide appears:
    {peptide_x:    [cycle_1, ..., cycle_n]}
    '''
    peptides_appearances_by_cycle = {}
    
    for peptide in base_cycle_sorted_peptides_list:
        peptides_appearances_by_cycle[peptide] = []
        for cycle in peptides_occurences_by_cycle:
            if peptide in peptides_occurences_by_cycle[cycle]:
                peptides_appearances_by_cycle[peptide] += [cycle]
    return peptides_appearances_by_cycle


# Define dnas_appearances_by_cycle, 
# dnas_appearances_by_cycle = 
def dnas_appearances_by_cycle(
        base_cycle_sorted_dnas_list: list,
        dnas_occurences_by_cycle):
    '''
    which returns for each dna in selection a list of cycles in which this dna appears:
    {dna_x:    [cycle_1, ..., cycle_N]}
    '''
    dnas_appearances_by_cycle = {}
    
    for dna in base_cycle_sorted_dnas_list:
        dnas_appearances_by_cycle[dna] = []
        for cycle in dnas_occurences_by_cycle:
            if dna in dnas_occurences_by_cycle[cycle]:
                dnas_appearances_by_cycle[dna] += [cycle]
    return dnas_appearances_by_cycle


# Define selection_summaryReport, which, for the TopNpeptidesNumber, returns a summary table (.txt) and provides a summary grpaph (.png).
def selection_summaryReport(data_directory_path,
                           base_selection_cycle_number,
                           TopNpeptidesNumber,
                           start_sequence,
                           stop_sequence,
                           cdna_min_length,
                           cdna_max_length,
                           quality_score,
                           FileName):
    
    today = todays_date() 
    
    selection_summaryFileNameCSV = str(today) + 'selection_summary' + FileName + '.csv'
    selection_summaryReportFile = open(selection_summaryFileNameCSV, 'w')
    
    selection_summary = complete_selection_summary(data_directory_path, start_sequence, stop_sequence, cdna_min_length, cdna_max_length, quality_score)
    SortedcyclesList = sorted(selection_summary.keys())
    
    peptides_by_cycle = peptides_occurences_by_cycle(data_directory_path, cdna_min_length, cdna_max_length, quality_score)
    Totalpeptides_by_cycle = total_reads_by_cycle(data_directory_path, cdna_min_length, cdna_max_length, quality_score)
    
    base_cycle_sorted_peptides = base_cycle_sorted_peptides_list(data_directory_path, base_selection_cycle_number)
    BasecycleTopSortedpeptides = base_cycle_sorted_peptides[0 : (TopNpeptidesNumber)]

    BasecyclepeptidesRank = peptides_rank_in_base_cycle(data_directory_path, base_selection_cycle_number)
        
    selection_summaryReportFile.write('peptide sequence' + ',' +
                                     'rank (#)' + ',' +
                                     'cdna mutants' + ',')
    for cycle in SortedcyclesList:
        selection_summaryReportFile.write('C' +
                                         str(cycle) +
                                         ' count (#) [frequency(%)]' + ',')
    selection_summaryReportFile.write('\n')
    
    for peptide in BasecycleTopSortedpeptides:
        BasecyclepeptideFraction = float((peptides_by_cycle[cycle].get(peptide, 0)))/float(Totalpeptides_by_cycle[base_selection_cycle_number])
        peptide_rank = BasecyclepeptidesRank[peptide]
        Formatedpeptide = hamming_distance_based_formating(BasecycleTopSortedpeptides[0], peptide)
        peptidecdnaMutants = len(selection_summary[base_selection_cycle_number][peptide])
        selection_summaryReportFile.write(Formatedpeptide + ',' +
                                         str(peptide_rank) + ',' +
                                         str(peptidecdnaMutants) + ',')
            
        for cycle in SortedcyclesList:
            peptideFraction = float((peptides_by_cycle[cycle].get(peptide, 0)))/float(Totalpeptides_by_cycle[cycle])

            BaseFraction = peptideFraction
            
            selection_summaryReportFile.write(str(peptides_by_cycle[cycle].get(peptide, 0)) +
                                             ' [' + '{:.1%}'.format(peptideFraction) + ']' + ',')
            
            BaseFraction = peptideFraction
        selection_summaryReportFile.write('\n')
        
    selection_summaryReportFile.write('total count (#)' + ',' + ',')
    for cycle in SortedcyclesList:
        selection_summaryReportFile.write(str(Totalpeptides_by_cycle[cycle]) + ',')
    selection_summaryReportFile.write('\n\n\n')
            
    selection_summaryReportFile.close()
    
#-------------------------------------------------------------------------------
   
    # Create a figure of size 8x6 inches, 500 dots per inch
    plt.figure(figsize = (8, 6),
               dpi = 500)
    # Create 'ggplot' style
    plt.style.use('fivethirtyeight')
    # Create a new subplot from a grid of 1x1
    Graph = plt.subplot(1, 1, 1)
    
    Xs = []
    Ys = []

    # Map colours onto lines  
    cNorm  = matplotlib.colors.Normalize(
        vmin = 0,
        vmax = TopNpeptidesNumber - 1)
    scalarMap = matplotlib.cm.ScalarMappable(
        norm = cNorm,
        cmap = 'gist_rainbow')
    
    peptideLabels = []
    
    for peptide in BasecycleTopSortedpeptides:
    #for peptide in Top24peptidesKDs:
        peptidesFractions_BY_cycle = []
        for cycle in SortedcyclesList:
            peptidesFractions_BY_cycle += [float((peptides_by_cycle[cycle].get(peptide, 0)))/float(Totalpeptides_by_cycle[cycle])]
        
        x = SortedcyclesList
        y = peptidesFractions_BY_cycle
        Xs += x
        Ys += y
        
        peptide_rank = BasecyclepeptidesRank[peptide]
        #peptideColour = scalarMap.to_rgba(peptide_rank)
        peptideColour = scalarMap.to_rgba(BasecycleTopSortedpeptides.index(peptide))
        Formatedpeptide = hamming_distance_based_formating(BasecycleTopSortedpeptides[0], peptide)
        
        peptideLabel =  Formatedpeptide + ' (' + str(peptide_rank) +')'
        
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
    
    plt.xlabel("Selection Cycle (#)",
               fontsize = 10)
    plt.ylabel("peptide Fraction",
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
    
    selection_summaryFileNamePNG = str(today) + 'selection_summary' + FileName + '.png'
    
    plt.savefig(selection_summaryFileNamePNG,
                bbox_extra_artists = [legend],
                bbox_inches = 'tight',
                dpi = 300)
    plt.show()
    plt.close()


def dnaMutantsAnalysis(data_directory_path,
                       base_selection_cycle_number,
                       TopNpeptidesNumber,
                       start_sequence,
                       stop_sequence,
                       FileName):
    
    today = todays_date() 
    
    dnaMutantsAnalysisFileNameCSV =  str(today) + 'dnasMutantsAnalysis' + FileName + '.csv'
    dnaMutantsAnalysisFile = open(dnaMutantsAnalysisFileNameCSV, 'w')
    
    selection_summary = complete_selection_summary(data_directory_path, start_sequence, stop_sequence, cdna_min_length, cdna_max_length, quality_score)
    SortedcyclesList = sorted(selection_summary.keys())
    
    peptides_by_cycle = peptides_occurences_by_cycle(data_directory_path, cdna_min_length, cdna_max_length, quality_score)
    Totalpeptides_by_cycle = total_reads_by_cycle(data_directory_path, cdna_min_length, cdna_max_length, quality_score)
    
    base_cycle_sorted_peptides = base_cycle_sorted_peptides_list(data_directory_path, base_selection_cycle_number)
    BasecycleTopSortedpeptides = base_cycle_sorted_peptides[0 : (TopNpeptidesNumber)]
    
    dnaClones_BY_cycle_BY_peptide = dna_clones_occurences_by_cycle_by_peptide(data_directory_path, cdna_min_length, cdna_max_length, quality_score)
    
    dnaMutantsAnalysisFile.write('peptide sequence' + ',')
    for cycle in SortedcyclesList:
        dnaMutantsAnalysisFile.write('cycle # ' + str(cycle) + ' dna clones (#)' + ',')
    dnaMutantsAnalysisFile.write('\n')
    
    for peptide in BasecycleTopSortedpeptides:
        dnaMutantsAnalysisFile.write(peptide + ',')
        for cycle in SortedcyclesList:
            dnaMutantsAnalysisFile.write(str(dnaClones_BY_cycle_BY_peptide[cycle].get(peptide, 0)) + ',')
        dnaMutantsAnalysisFile.write('\n')
    dnaMutantsAnalysisFile.close()
    
#-------------------------------------------------------------------------------        

    # Create a figure of size 8x6 inches, 500 dots per inch
    plt.figure(figsize = (8, 6),
               dpi = 500)
    # Create 'ggplot' style
    plt.style.use('fivethirtyeight')
    # Create a new subplot from a grid of 1x1
    Graph = plt.subplot(1, 1, 1)
    
#    peptidednaClonesNumber_IN_Basecycle = []
#    peptideoccurence_IN_Basecycle = []

    # Map colours onto lines
    cNorm  = matplotlib.colors.Normalize(
        vmin = 0,
        vmax = len(base_cycle_sorted_peptides) - 1)
    scalarMap = matplotlib.cm.ScalarMappable(
        norm = cNorm,
        cmap = 'gist_rainbow')

    cycleIndex = base_selection_cycle_number

    Xs = []
    Ys = []        
    for peptide in dnaClones_BY_cycle_BY_peptide[cycleIndex]:
        peptidednaClonesNumber_IN_Basecycle = math.log(dnaClones_BY_cycle_BY_peptide[cycleIndex].get(peptide, 0), 2)
        peptideoccurence_IN_Basecycle = math.log(peptides_by_cycle[cycleIndex].get(peptide, 0), 2)
        
        peptideColour = scalarMap.to_rgba(base_cycle_sorted_peptides.index(peptide))
    
        x = peptidednaClonesNumber_IN_Basecycle
        y = peptideoccurence_IN_Basecycle
        Xs += [x]
        Ys += [y]
        
        plt.plot(x, y,
                'o',
                c = peptideColour,
                ms = 5.0,
                mew = 0.1,
                mec = '#191919')

    XMin = min(Xs) - 0.05*(max(Xs) - min(Xs))
    XMax = max(Xs) + 0.05*(max(Xs) - min(Xs))
    YMin = min(Ys) - 0.05*(max(Ys) - min(Ys))
    YMax = max(Ys) + 0.05*(max(Ys) - min(Ys))
    
    plt.axis([XMin, XMax, YMin, YMax])
        
    XLabel = 'log$2$ (dna Clones #)' #$_$ makes subscript possible
    plt.xlabel(XLabel, fontsize = 14)
    YLabel = 'log$2$ (peptide occurence)'
    plt.ylabel(YLabel, fontsize = 14)
    
    legend = plt.legend(base_cycle_sorted_peptides,
                        loc = 'upper center',
                        bbox_to_anchor = (0.5, -0.15),
                        fancybox = True,
                        shadow = False,
                        ncol = 4)
    
    dnaClonesAnalysisFileNamePNG = str(today) + 'dnasMutantsAnalysisRegression' + 'R' + str(cycle) + FileName + '.png'
    
    plt.savefig(dnaClonesAnalysisFileNamePNG, bbox_extra_artists=[legend], bbox_inches='tight', dpi = 300)
    plt.show()
    plt.close()


def peptidesRelatednessAnalysis(data_directory_path,
                                base_selection_cycle_number,
                                TopNpeptidesNumber,
                                start_sequence,
                                stop_sequence,
                                cdna_min_length,
                                cdna_max_length,
                                quality_score,
                                FileName):
    
    # to extract todays_date
    today = todays_date()
    
    # to collect dnas-based summary information By_cycle
    dnas_by_cycle = dnas_occurences_by_cycle(data_directory_path, cdna_min_length, cdna_max_length, quality_score)
    Totaldnas_by_cycle = total_reads_by_cycle(data_directory_path, cdna_min_length, cdna_max_length, quality_score)
    BasecycleSorteddnas = base_cycle_sorted_dnas_list(data_directory_path, base_selection_cycle_number)
    dnasAppearances = dnas_appearances_by_cycle(BasecycleSorteddnas, dnas_by_cycle)
    
    # to collect peptides-based summary information By_cycle
    peptides_by_cycle = peptides_occurences_by_cycle(data_directory_path, cdna_min_length, cdna_max_length, quality_score)
    Totalpeptides_by_cycle = total_reads_by_cycle(data_directory_path, cdna_min_length, cdna_max_length, quality_score)
    base_cycle_sorted_peptides = base_cycle_sorted_peptides_list(data_directory_path, base_selection_cycle_number)
    peptidesAppearances = peptides_appearances_by_cycle(base_cycle_sorted_peptides, peptides_by_cycle)
    
    selection_summary = complete_selection_summary(data_directory_path, start_sequence, stop_sequence, cdna_min_length, cdna_max_length, quality_score)
    SortedcyclesList = sorted(selection_summary.keys())
    
    # to create a disjoint graph (Forest), based on dnas in the Basecycle (joint subgraphs are the Trees and the unique dna sequences are the Leaves)
    BasecyclednasForest = nx.Graph()
    # to add nodes (Leaves, unique dna sequences) to the BasecyclednasForest disjoint graph
    BasecyclednasForest.add_nodes_from(BasecycleSorteddnas)
    # to add edges (Twigs, dna-to-dna connections based on the hamming distance between unique dna sequences) to the BasecyclednasForest so that disjoint graphs (stand alone Trees) can be identified
    UsedNodes = []
    for dna1 in BasecycleSorteddnas:
        UsedNodes += [dna1]
        for dna2 in BasecycleSorteddnas:
            if dna2 not in UsedNodes and hamming_distance(dna1, dna2) == 1:
                BasecyclednasForest.add_edge(dna1,dna2,
                                                MutationsNumber = 1)
    # to extract individual joint subgraphs (stand alone Trees) from the disjoint graph (Forest)
    BasecyclednasTrees = list(nx.connected_component_subgraphs(BasecyclednasForest, copy = True))
    
    # to create a peptideSummarydnaPerspectiveCSV file
    peptidesSummaryFileNameCSV =  str(today) + 'peptideFamiliesSummary' + FileName + '.csv'
    peptidesSummaryFile = open(peptidesSummaryFileNameCSV, 'w')
    
    # to convert list of dnas Trees into a list of peptides Trees Leaves
    peptidesTreesLeaves = []
    for dnasTree in BasecyclednasTrees:
        peptideLeaves = []
        for dna in dnasTree:
            peptide = translation(dna)
            if peptide not in peptideLeaves:
                peptideLeaves += [peptide]
        peptidesTreesLeaves += [peptideLeaves]
    # to sort the resulting list of lists from the largest to smallest
    peptidesTreesLeaves.sort(key = len, reverse = True)
    
    # to fix the coordinates of the origin of the graph
    Positions = {}
    X_0_Coordinate = 1
    Y_0_Coordinate = 0
    Y_X0_Coordinate = 0
    
    TreesXCoordinates = []
    
    peptideGraphFigure = plt.figure()
    peptideGraph = peptideGraphFigure.add_subplot(1, 1, 1)
    
    #to introduce peptideFamilyCounter
    MultiplepeptideFamilyCounter = 0
    #to introduce SinglepeptideFamilyCounter
    SinglepeptideFamilyCounter = 0
    #to introduce peptideFamilySize
    peptideTreeSize = []
    
    
    # to create a tree for each set of peptides Trees Leaves 
    for peptideLeaves in peptidesTreesLeaves:
        
        peptideTree = nx.Graph()
        # to convert each peptide (Leave) into a node of a peptide graph (peptideLeave on a peptideTree)
        
        for peptide in peptideLeaves:
            peptideTree.add_node(peptide,
                                        occurence = peptides_by_cycle[base_selection_cycle_number][peptide],
                                        FirstAppearance = min(peptidesAppearances[peptide]))
        # to join the peptide nodes of a graph by edges (Twigs)
        for peptide1 in peptideLeaves:
            for peptide2 in peptideLeaves:
                if hamming_distance(peptide1, peptide2) == 1:
                    peptideTree.add_edge(peptide1,peptide2)
        
        # to identify the Rootpeptide of a peptideTree graph
        Treepeptidesoccurences = nx.get_node_attributes(peptideTree, 'occurence')              
        Rootpeptide = max(Treepeptidesoccurences, key=Treepeptidesoccurences.get)
        
        # to create a dictionary holder for peptide and their properties (Predecessor and Occurrence)
        Treepeptides = {}
        Treepeptides[Rootpeptide] = [0, '', 0, peptideTree.node[Rootpeptide]['occurence'], peptideTree.node[Rootpeptide]['FirstAppearance']]
        TreepeptidesList = list(peptideTree.nodes())
        TreepeptidesList.remove(Rootpeptide)

        for peptide in TreepeptidesList:
            peptidePredecessor = nx.shortest_path(peptideTree, source = peptide, target = Rootpeptide, weight = None)[1]
            # Predecessoroccurence can be used to sort the peptides, but does not seem to be useful
            Predecessoroccurence = peptideTree.node[peptidePredecessor]['occurence']
            peptideoccurence = peptideTree.node[peptide]['occurence']

            Treepeptides[peptide] = [peptidePredecessor, Predecessoroccurence, peptideoccurence]
            
        
        # to sort peptides in a peptideTree by their distance to the Rootpeptide
        peptides_BY_DistanceToTheRoot = {}
        for peptide in peptideTree.nodes():
            DistanceToTheRoot = nx.shortest_path_length(peptideTree, source = peptide, target = Rootpeptide, weight = None)
            if DistanceToTheRoot not in peptides_BY_DistanceToTheRoot:
                peptides_BY_DistanceToTheRoot[DistanceToTheRoot] = [peptide]
            else:
                peptides_BY_DistanceToTheRoot[DistanceToTheRoot] += [peptide]
        
        # to identify the largest group of equidistanced peptides 
        MaxpeptidesNumber = max(map(lambda k: len(peptides_BY_DistanceToTheRoot[k]), peptides_BY_DistanceToTheRoot))

        Sortedpeptides_BY_DistanceToTheRoot = {}
        # to sort peptides by their distance to the Rootpeptide
        for DistanceToTheRoot in peptides_BY_DistanceToTheRoot:

            Equidistantpeptides = peptides_BY_DistanceToTheRoot[DistanceToTheRoot]

            Equidistantpeptides = sorted(Equidistantpeptides, key = lambda peptide: (Treepeptides[peptide][2]), reverse = True)
            # Predecessoroccurence can be used to sort the peptides, but does not seem to be useful
            # Equidistantpeptides = sorted(Equidistantpeptides, key = lambda peptide: (Treepeptides[peptide][1]), reverse = True)
            Equidistantpeptides = sorted(Equidistantpeptides, key = lambda peptide: (Treepeptides[peptide][0]), reverse = False)

            AdditionalElements = MaxpeptidesNumber - len(Equidistantpeptides)
            Sortedpeptides_BY_DistanceToTheRoot[DistanceToTheRoot] = Equidistantpeptides + AdditionalElements * ['']

            if len(peptideTree.nodes()) > 1:
                for peptide in Equidistantpeptides:
                    XCoordinate = X_0_Coordinate + DistanceToTheRoot
                    YCoordinate = Y_0_Coordinate - Equidistantpeptides.index(peptide)
                    Positions[peptide] = (XCoordinate, YCoordinate)
                    
                                    
            elif len(peptideTree.nodes()) == 1:
                for peptide in Equidistantpeptides:
                    XCoordinate = 0
                    YCoordinate = Y_X0_Coordinate
                    Positions[peptide] = (XCoordinate, YCoordinate)
                    

        #BasecyclepeptidesGraph = nx.Graph()    
        #BasecyclepeptidesGraph.add_nodes_from(base_cycle_sorted_peptides)

        Sizes = []
        for peptide in peptideTree.nodes():
            Sizes.append(math.log(peptides_by_cycle[base_selection_cycle_number][peptide], 2) + 5)

        Colours = []
        for peptide in peptideTree.nodes():
            Colours.append(min(peptidesAppearances[peptide]))
        
        XSpan = max(map(lambda peptide: Positions[peptide][0], Positions)) - min(map(lambda peptide: Positions[peptide][0], Positions))
        YSpan = max(map(lambda peptide: Positions[peptide][1], Positions)) - min(map(lambda peptide: Positions[peptide][1], Positions))
                            
        XMin = min(map(lambda peptide: Positions[peptide][0], Positions)) - 0.01 * XSpan
        XMax = max(map(lambda peptide: Positions[peptide][0], Positions)) + 0.01 * XSpan
        YMin = min(map(lambda peptide: Positions[peptide][1], Positions)) - 0.02 * YSpan
        YMax = max(map(lambda peptide: Positions[peptide][1], Positions)) + 0.02 * YSpan
        
        
        
        NumberOfColours = len(SortedcyclesList)
        
        ColourMap = plt.get_cmap('Paired', NumberOfColours)
        
        nx.draw_networkx(peptideTree,
                        pos = Positions,
                        node_size = Sizes,
                        node_color = Colours,
                        cmap = ColourMap,
                        linewidths = 0.2,
                        width = 0.2,
                        with_labels = False,
                        #font_size = 6,
                        vmin = min(SortedcyclesList),
                        vmax = max(SortedcyclesList))
                    
        if len(peptideTree.nodes()) > 1:
            for DistanceToTheRoot in Sortedpeptides_BY_DistanceToTheRoot:
                peptidesSummaryFile.write(str(DistanceToTheRoot) + ' mutations' + ',' + 'frequency' + ',' + 'rank' + ',')
            peptidesSummaryFile.write('\n')

        for i in range(MaxpeptidesNumber):
            for MutationsNumber in Sortedpeptides_BY_DistanceToTheRoot:                        
                peptide = Sortedpeptides_BY_DistanceToTheRoot[MutationsNumber][i]

                if peptide != '':
                    Formatedpeptide = hamming_distance_based_formating(Rootpeptide, peptide)
                    peptide_rank = str(base_cycle_sorted_peptides.index(peptide) + 1)
                    #ClonesNumber = str(len(peptideTree.neighbors(peptide)))
                    peptideFraction = ('{:.2%}'.format(float((peptides_by_cycle[base_selection_cycle_number].get(peptide, 0)))/float(Totalpeptides_by_cycle[base_selection_cycle_number])))
                else:
                    Formatedpeptide = ''
                    #ClonesNumber = ''
                    peptide_rank = ''
                    peptideFraction = ''

                peptidesSummaryFile.write(Formatedpeptide + ',' +
                            peptideFraction + ',' +
                            peptide_rank + ',')
                            #ClonesNumber + ',')
            peptidesSummaryFile.write('\n')
        
        
        if len(peptideTree.nodes()) > 1:
            TreesXCoordinates += [X_0_Coordinate]
            X_0_Coordinate += max(peptides_BY_DistanceToTheRoot.keys()) + 1
            MultiplepeptideFamilyCounter += 1
            peptideTreeSize += [len(peptideTree.nodes())]
            
            

        if len(peptideTree.nodes()) == 1:
            Y_X0_Coordinate -= 1
            SinglepeptideFamilyCounter += 1

        peptidesSummaryFile.write('\n')
                    
    peptidesSummaryFile.close()
    
    #plt.axis('off')
    plt.axis([XMin, XMax, YMin, YMax])

    peptideLegendColour = peptideGraphFigure.add_subplot(1, 1, 1)
        
    ColourMap = plt.get_cmap('Paired', NumberOfColours)
    peptideLegendColours = SortedcyclesList
    
    LegendDotsX = XMax - 0.3 * XMax
    YIncrement = - 0.03 * YMin
    #print (YIncrement)
    
    LegendColourDotsX = np.array([LegendDotsX] * NumberOfColours)
    #print (LegendColourDotsX)
    FirstYColours = YMin + 12 * YIncrement
    #print (FirstYColours)
    LastYColours = YMin + (12 + NumberOfColours) * YIncrement
    #print (LastYColours)
    LegendColourDotsY = np.linspace(FirstYColours, LastYColours, NumberOfColours, endpoint = False)
    #print (LegendColourDotsY)
    
    peptideLegendColour.scatter(x = LegendColourDotsX,
                                y = LegendColourDotsY,
                                s = 15,
                                c = peptideLegendColours,
                                cmap = ColourMap,
                                linewidths = 0.2)
    
    ColourLabels = SortedcyclesList
#     this way of setting the colours seems to be redundant
#     ColourLabels = ['{0}'.format(i) for i in range(NumberOfColours)]

    for label, x, y in zip(ColourLabels, LegendColourDotsX, LegendColourDotsY):
        plt.annotate(label, xy = (x, y), xytext = (5, 0),
                     textcoords = 'offset points',
                     fontsize = 5,
                     ha = 'left', va = 'center')
    plt.text(x = LegendDotsX, y = (max(LegendColourDotsY) + YIncrement),
             s = 'first-appearance cycle #',
             fontsize = 5)
    #plt.axis('off')

    peptideLegendSize = peptideGraphFigure.add_subplot(1, 1, 1)

    Size = []
    for i in [1, 10, 100, 1000, 10000]:
        Size.append(math.log(i, 2) + 5)

    LegendSizeDotsX = np.array([LegendDotsX] * 5)
    FirstYSizes = YMin + 5 * YIncrement
    LastYSizez = YMin + 10 * YIncrement
    LegendSizeDotsY = np.linspace(FirstYSizes, LastYSizez, 5, endpoint = False)
    peptideLegendSize.scatter(x = LegendSizeDotsX,
                              y = LegendSizeDotsY,
                              s = Size,
                              c = 'w',
                              linewidths = 0.2)

    SizeLabels = ['{0}'.format(i) for i in [1, 10, 100, 1000, 10000]]

    for label, x, y in zip(SizeLabels, LegendSizeDotsX, LegendSizeDotsY):
        plt.annotate(label, xy = (x, y), xytext = (5, 0),
                     textcoords = 'offset points',
                     fontsize = 5,
                     ha = 'left', va = 'center')
    plt.text(x = LegendDotsX, y = (max(LegendSizeDotsY) - 0.03 * YMin),
             s = 'frequency in the last cycle',
             fontsize = 5)
    
    
    
    #FamilySizeLabels = ['{0}'.format(i) for i in peptideTreeSize]

    #for label, x, y in zip(FamilySizeLabels, LegendSizeDotsX, LegendSizeDotsY):
    #    plt.annotate(label, xy = (x, y), xytext = (5, 0),
    #                 textcoords = 'offset points',
    #                 fontsize = 5,
    #                 ha = 'left', va = 'center')
    for i in range(len(peptideTreeSize)):
        plt.text(x = TreesXCoordinates[i], y = YIncrement,
                 s = peptideTreeSize[i],
                 fontsize = 5)
    
    
    
    
    
    plt.text(x = LegendDotsX, y = YMin + 3 * YIncrement,
             s = ('total # unique peptide sequence ' + str(len(base_cycle_sorted_peptides))),
             fontsize = 5)
    plt.text(x = LegendDotsX, y = YMin + 2 * YIncrement,
             s = 'single-member peptide family # ' + str(SinglepeptideFamilyCounter),
             fontsize = 5)
    plt.text(x = LegendDotsX, y = YMin + 1 * YIncrement,
             s = 'multiple-member peptide family # ' + str(MultiplepeptideFamilyCounter),
             fontsize = 5)
    
    plt.axis('off')

    peptidesSummaryFileNamePNG = str(today) + 'peptideFamiliesSummary' + FileName + '.png'
    plt.savefig(peptidesSummaryFileNamePNG, dpi = 500)
    
    
    fig = plt.gcf()
    SizeInches = fig.get_size_inches()*fig.dpi
    SizeDots = fig.get_size_inches()
    
    #print (SizeInches)
    #print (SizeDots)
    
    
    #print (XMin)
    #print (XMax)
    #print (YMin)
    #print (YMax)
    
    #print (peptideTreeSize)
    #print (len(peptideTreeSize))
    #print (TreesXCoordinates)
    #print (len(TreesXCoordinates))
    
    plt.show()
    plt.close()