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

import global_parameters as gs

# =============================================================================
def todays_date():
    '''
    returns todays date in format YYYYMMDD
    '''
    return datetime.date.today().strftime('%Y%m%d')


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
    (1) its length AND
    (2) quality score of each nucleotide;
    uses ONLY ONE stop sequence;
    returns ONLY ONE coding sequence
    '''
    quality_score_string = gs.QUALITY_SCORE_STRING
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
    translation_code = gs.TRANSLATION_CODE
    transcription_code = gs.TRANSCRIPTION_CODE
    
    # Convert DNA to RNA.
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


# Define single_selection_cycle_summary function,
# which returns the occurrence of coding_sequences, grouped by peptide
# {peptideY: {coding_sequence_YZ: occurrence_YZ}}

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
    # {peptideY:    {coding_sequence_YZ:    occurrence_YZ}}
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
    returns the occurrence of coding sequences
    grouped by peptide, and grouped by selection cycle:
    {Selectioncycle_X:    {peptideXY:    {Codingdna_XYZ:    occurrence_XYZ}}}
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


def peptides_occurrences_by_cycle(
        data_directory_path: str,
        cdna_min_length: int,
        cdna_max_length: int,
        quality_score: int,
        start_sequence: str,
        stop_sequence: str,
        ):
    '''
    returns the occurrences of peptides groupped by cycle:
    {cycle_X:    {peptideXY:    occurrence_XY}}
    '''
    selection_summary = complete_selection_summary(
        data_directory_path,
        start_sequence,
        stop_sequence,
        cdna_min_length,
        cdna_max_length,
        quality_score)
    
    peptides_occurrences_by_cycle = {}
    for cycle in selection_summary:
        peptides_occurrences_in_cycle = {}
        for peptide in selection_summary[cycle]:
            peptides_occurrences_in_cycle[peptide] = sum(
                selection_summary[cycle][peptide].values())
        peptides_occurrences_by_cycle[cycle] = peptides_occurrences_in_cycle
        
    return peptides_occurrences_by_cycle


def dnas_occurrences_by_cycle(
        data_directory_path: str,
        cdna_min_length: int,
        cdna_max_length: int,
        quality_score: int,
        start_sequence: str,
        stop_sequence: str,
        ):
    '''
    returns the occurrences of dnas groupped by cycle:
    {cycle_X:    {dna_XY:    occurrence_XY}}
    '''
    selection_summary = complete_selection_summary(
        data_directory_path,
        start_sequence,
        stop_sequence,
        cdna_min_length,
        cdna_max_length,
        quality_score)
    
    dnas_occurrences_by_cycle = {}
    for cycle in selection_summary:
        dnas_occurrences_in_cycle = {}
        for peptide in selection_summary[cycle]:
            for dna in selection_summary[cycle][peptide]:
                dnas_occurrences_in_cycle[dna] = selection_summary[cycle][peptide][dna]
        dnas_occurrences_by_cycle[cycle] = dnas_occurrences_in_cycle

    return dnas_occurrences_by_cycle


def total_reads_by_cycle(
        data_directory_path: str,
        cdna_min_length: int,
        cdna_max_length: int,
        quality_score: int,
        start_sequence: str,
        stop_sequence: str,
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

    peptides_by_cycle = peptides_occurrences_by_cycle(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        quality_score,
        start_sequence,
        stop_sequence,
        )
    
    total_reads_by_cycle = {}
    for cycle in selection_summary:
        total_reads_by_cycle[cycle] = sum(peptides_by_cycle[cycle].values())
        
    return total_reads_by_cycle


def base_cycle_sorted_peptides_list(
        data_directory_path: str,
        base_selection_cycle_number: int,
        cdna_min_length: int,
        cdna_max_length: int,
        quality_score: int,
        start_sequence: str,
        stop_sequence: str,
        ):
    '''
    returns list of peptides in base cycle sorted by their occurrence:
    [peptide_1, ..., peptide_N];
    occurrence(peptide_1) > ... > occurrence(peptide_N)
    '''
    peptides_by_cycle = peptides_occurrences_by_cycle(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        quality_score,
        start_sequence,
        stop_sequence,
        )
            
    peptides_occurrences_in_base_cycle = peptides_by_cycle[
        base_selection_cycle_number]
    base_cycle_sorted_peptides_list = sorted(
        peptides_occurrences_in_base_cycle,
        key=peptides_occurrences_in_base_cycle.get,
        reverse=True)
    return base_cycle_sorted_peptides_list


def base_cycle_sorted_dnas_list(
        data_directory_path: str,
        base_selection_cycle_number: int,
        cdna_min_length: int,
        cdna_max_length: int,
        quality_score: int,
        start_sequence: str,
        stop_sequence: str,
        ):
    '''
    returns list of dnas in base cycle sorted by their occurrence:
    [dna_1, ..., dna_n]; occurrence(dna_1) > ... > occurrence(dna_n)
    '''
    dnas_by_cycle = dnas_occurrences_by_cycle(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        quality_score,
        start_sequence,
        stop_sequence,
        )
            
    dnas_occurrences_in_base_cycle = dnas_by_cycle[base_selection_cycle_number]
    base_cycle_sorted_dnas_list = sorted(
        dnas_occurrences_in_base_cycle,
        key=dnas_occurrences_in_base_cycle.get,
        reverse=True)
    
    return base_cycle_sorted_dnas_list


def peptides_rank_in_base_cycle(
        data_directory_path: str,
        base_selection_cycle_number: int,
        cdna_min_length: int,
        cdna_max_length: int,
        quality_score: int,
        start_sequence: str,
        stop_sequence: str,
        ):

    peptides_by_cycle = peptides_occurrences_by_cycle(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        quality_score,
        start_sequence,
        stop_sequence,
        )

    base_cycle_sorted_peptides = base_cycle_sorted_peptides_list(
        data_directory_path,
        base_selection_cycle_number,
        cdna_min_length,
        cdna_max_length,
        quality_score,
        start_sequence,
        stop_sequence,
        )
    
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


# occurrences is a bad word-choice here
def dna_clones_occurrences_by_cycle_by_peptide(
        data_directory_path: str,
        cdna_min_length: int,
        cdna_max_length: int,
        quality_score: int,
        start_sequence: str,
        stop_sequence: str,
        ):
    '''
    returns number of clones for each peptide groupped by cycle:
    {Selectioncycle_X:    {peptideXY:    dnaClonesoccurrences}}
    '''
    selection_summary = complete_selection_summary(
        data_directory_path,
        start_sequence,
        stop_sequence,
        cdna_min_length,
        cdna_max_length,
        quality_score)
    
    dna_clones_occurrences_by_cycle_by_peptide = {}
    for cycle in selection_summary:
        dna_clones_occurrences_by_peptide = {}
        for peptide in selection_summary[cycle]:
            dna_clones_occurrences_by_peptide[peptide] = len(
                selection_summary[cycle][peptide])
        dna_clones_occurrences_by_cycle_by_peptide[cycle] = dna_clones_occurrences_by_peptide
        
    return dna_clones_occurrences_by_cycle_by_peptide


def peptides_appearances_by_cycle(
        base_cycle_sorted_peptides_list: list,
        peptides_occurrences_by_cycle):
    '''
    returns for each peptide in selection a list of cycles in which this peptide appears:
    {peptide_x:    [cycle_1, ..., cycle_n]}
    '''
    peptides_appearances_by_cycle = {}
    
    for peptide in base_cycle_sorted_peptides_list:
        peptides_appearances_by_cycle[peptide] = []
        for cycle in peptides_occurrences_by_cycle:
            if peptide in peptides_occurrences_by_cycle[cycle]:
                peptides_appearances_by_cycle[peptide] += [cycle]
    return peptides_appearances_by_cycle


# Define dnas_appearances_by_cycle, 
# dnas_appearances_by_cycle = 
def dnas_appearances_by_cycle(
        base_cycle_sorted_dnas_list: list,
        dnas_occurrences_by_cycle):
    '''
    which returns for each dna in selection a list of cycles in which this dna appears:
    {dna_x:    [cycle_1, ..., cycle_N]}
    '''
    dnas_appearances_by_cycle = {}
    
    for dna in base_cycle_sorted_dnas_list:
        dnas_appearances_by_cycle[dna] = []
        for cycle in dnas_occurrences_by_cycle:
            if dna in dnas_occurrences_by_cycle[cycle]:
                dnas_appearances_by_cycle[dna] += [cycle]
    return dnas_appearances_by_cycle


# Define selection_summary_report, which, for the top_n_peptides,
# returns a summary table (.txt) and provides a summary grpaph (.png).
def selection_summary_report(
        data_directory_path: str,
        base_selection_cycle_number: int,
        top_n_peptides: int,
        start_sequence: str,
        stop_sequence: str,
        cdna_min_length: int,
        cdna_max_length: int,
        quality_score: int,
        file_name: str
        ):
    
    today = todays_date() 
    
    selection_summary_csv = f"{today}_selection_summary_{file_name}.csv"
    selection_summary_report_file = open(selection_summary_csv, 'w')
    
    selection_summary = complete_selection_summary(
        data_directory_path,
        start_sequence,
        stop_sequence,
        cdna_min_length,
        cdna_max_length,
        quality_score)

    sorted_cycles_list = sorted(selection_summary.keys())
    
    peptides_by_cycle = peptides_occurrences_by_cycle(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        quality_score,
        start_sequence,
        stop_sequence,
        )

    total_peptides_by_cycle = total_reads_by_cycle(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        quality_score,
        start_sequence,
        stop_sequence,
        )
    
    base_cycle_sorted_peptides = base_cycle_sorted_peptides_list(
        data_directory_path,
        base_selection_cycle_number,
        cdna_min_length,
        cdna_max_length,
        quality_score,
        start_sequence,
        stop_sequence,
        )
    base_cycle_top_sorted_peptides = base_cycle_sorted_peptides[
        0 : (top_n_peptides)]

    base_cycle_peptides_rank = peptides_rank_in_base_cycle(
        data_directory_path,
        base_selection_cycle_number,
        cdna_min_length,
        cdna_max_length,
        quality_score,
        start_sequence,
        stop_sequence,
        )
        
    selection_summary_report_file.write(
        f"peptide sequence,rank (#),cdna mutants,")
    for cycle in sorted_cycles_list:
        selection_summary_report_file.write(
            f"C{cycle}count (#) [frequency(%)],")
    selection_summary_report_file.write(f"\n")
    
    for peptide in base_cycle_top_sorted_peptides:
        base_cycle_peptide_fraction = float(
            (peptides_by_cycle[cycle].get(peptide, 0)))/float(
                total_peptides_by_cycle[base_selection_cycle_number])
        peptide_rank = base_cycle_peptides_rank[peptide]
        formated_peptide = hamming_distance_based_formating(
            base_cycle_top_sorted_peptides[0], peptide)
        peptide_cdna_mutants = len(
            selection_summary[base_selection_cycle_number][peptide])
        selection_summary_report_file.write(
            f"{formated_peptide},{peptide_rank},{peptide_cdna_mutants},")
            
        for cycle in sorted_cycles_list:
            peptide_fraction = float(
                (peptides_by_cycle[cycle].get(peptide, 0)))/float(
                    total_peptides_by_cycle[cycle])

            base_fraction = peptide_fraction
            
            selection_summary_report_file.write(
                f"{str(peptides_by_cycle[cycle].get(peptide, 0))}"
                f" [{peptide_fraction:.1%}],")
            
            base_fraction = peptide_fraction
        selection_summary_report_file.write(f"\n")
        
    selection_summary_report_file.write(
        f"total count (#),,")
    for cycle in sorted_cycles_list:
        selection_summary_report_file.write(
            f"{str(total_peptides_by_cycle[cycle])},")
    selection_summary_report_file.write('\n\n\n')
            
    selection_summary_report_file.close()
    
# =============================================================================
   
    # Create a figure 8x6 inches, 500 dots per inch
    plt.figure(figsize = (8, 6),
               dpi = 500)
    # Create 'ggplot' style
    plt.style.use('fivethirtyeight')
    # Create a new subplot from a grid of 1x1
    Graph = plt.subplot(1, 1, 1)
    
    Xs = []
    Ys = []

    # Map colours onto lines  
    c_norm  = matplotlib.colors.Normalize(
        vmin = 0,
        vmax = top_n_peptides - 1)
    scalarMap = matplotlib.cm.ScalarMappable(
        norm = c_norm,
        cmap = 'gist_rainbow')
    
    peptide_labels = []
    
    for peptide in base_cycle_top_sorted_peptides:
    #for peptide in Top24peptidesKDs:
        peptidesFractions_BY_cycle = []
        for cycle in sorted_cycles_list:
            peptidesFractions_BY_cycle += [float(
                (peptides_by_cycle[cycle].get(peptide, 0)))/float(
                    total_peptides_by_cycle[cycle])]
        
        x = sorted_cycles_list
        y = peptidesFractions_BY_cycle
        Xs += x
        Ys += y
        
        peptide_rank = base_cycle_peptides_rank[peptide]
        #peptideColour = scalarMap.to_rgba(peptide_rank)
        peptideColour = scalarMap.to_rgba(
            base_cycle_top_sorted_peptides.index(peptide))
        formated_peptide = hamming_distance_based_formating
        (base_cycle_top_sorted_peptides[0], peptide)
        
        peptide_label =  formated_peptide + ' (' + str(peptide_rank) +')'
        
        #Set peptide_label
        peptide_labels += [peptide_label]
        
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
    
    legend = plt.legend(peptide_labels,
                        title = 'cyclic-peptide random region',
                        loc = 'upper center',
                        bbox_to_anchor = (0.5, -0.10),
                        fancybox = True,
                        shadow = False,
                        fontsize = 10,
                        ncol = 3)
    
    Graph.get_legend().get_title().set_size('small')
    
    selection_summary_png_path = f"{today}_selection_summary_{file_name}.png"
    
    plt.savefig(selection_summary_png_path,
                bbox_extra_artists = [legend],
                bbox_inches = 'tight',
                dpi = 300)
    plt.show()
    plt.close()


def dna_mutants_analysis(
        data_directory_path: str,
        base_selection_cycle_number: int,
        top_n_peptides: int,
        start_sequence: str,
        stop_sequence: str,
        file_name: str,
        cdna_min_length: int,
        cdna_max_length: int,
        quality_score: int,
        ):
        
    
    today = todays_date() 
    
    dna_mutants_analysis_csv =  f"{today}_dnas_mutants_analysis_{file_name}.csv"
    dna_mutants_analysis_file = open(dna_mutants_analysis_csv, 'w')
    
    selection_summary = complete_selection_summary(
        data_directory_path,
        start_sequence,
        stop_sequence,
        cdna_min_length,
        cdna_max_length,
        quality_score)
    sorted_cycles_list = sorted(selection_summary.keys())
    
    peptides_by_cycle = peptides_occurrences_by_cycle(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        quality_score,
        start_sequence,
        stop_sequence,
        )
    total_peptides_by_cycle = total_reads_by_cycle(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        quality_score,
        start_sequence,
        stop_sequence,
        )
    
    base_cycle_sorted_peptides = base_cycle_sorted_peptides_list(
        data_directory_path,
        base_selection_cycle_number,
        cdna_min_length,
        cdna_max_length,
        quality_score,
        start_sequence,
        stop_sequence,
        )
    base_cycle_top_sorted_peptides = base_cycle_sorted_peptides[0 : (top_n_peptides)]
    
    dna_clones_by_cycle_by_peptide = dna_clones_occurrences_by_cycle_by_peptide(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        quality_score,
        start_sequence,
        stop_sequence
        )
    
    dna_mutants_analysis_file.write("peptide sequence,")
    for cycle in sorted_cycles_list:
        dna_mutants_analysis_file.write(f"cycle # {cycle} dna clones (#),")
    dna_mutants_analysis_file.write("\n")
    
    for peptide in base_cycle_top_sorted_peptides:
        dna_mutants_analysis_file.write(f"{peptide}," )
        for cycle in sorted_cycles_list:
            dna_mutants_analysis_file.write(
                f"{str(dna_clones_by_cycle_by_peptide[cycle].get(peptide, 0)),}")
        dna_mutants_analysis_file.write("\n")
    dna_mutants_analysis_file.close()
    
#-------------------------------------------------------------------------------        

    # Create a figure of size 8x6 inches, 500 dots per inch
    plt.figure(figsize = (8, 6),
               dpi = 500)
    # Create 'ggplot' style
    plt.style.use('fivethirtyeight')
    # Create a new subplot from a grid of 1x1
    Graph = plt.subplot(1, 1, 1)
    
#    peptide_dna_clones_number_in_base_cycle = []
#    peptide_occurrence_in_base_cycle = []

    # Map colours onto lines
    c_norm  = matplotlib.colors.Normalize(
        vmin = 0,
        vmax = len(base_cycle_sorted_peptides) - 1)
    scalarMap = matplotlib.cm.ScalarMappable(
        norm = c_norm,
        cmap = 'gist_rainbow')

    cycle_index = base_selection_cycle_number

    Xs = []
    Ys = []        
    for peptide in dna_clones_by_cycle_by_peptide[cycle_index]:
        peptide_dna_clones_number_in_base_cycle = math.log(
            dna_clones_by_cycle_by_peptide[cycle_index].get(peptide, 0), 2)
        peptide_occurrence_in_base_cycle = math.log(
            peptides_by_cycle[cycle_index].get(peptide, 0), 2)
        
        peptideColour = scalarMap.to_rgba(
            base_cycle_sorted_peptides.index(peptide))
    
        x = peptide_dna_clones_number_in_base_cycle
        y = peptide_occurrence_in_base_cycle
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
    YLabel = 'log$2$ (peptide occurrence)'
    plt.ylabel(YLabel, fontsize = 14)
    
    legend = plt.legend(base_cycle_sorted_peptides,
                        loc = 'upper center',
                        bbox_to_anchor = (0.5, -0.15),
                        fancybox = True,
                        shadow = False,
                        ncol = 4)
    
    dna_clones_analysis_png_path = (
        f"{today}dnas_mutants_analysis_regression_C{cycle}_{file_name}.png")
    
    plt.savefig(
        dna_clones_analysis_png_path,
        bbox_extra_artists=[legend],
        bbox_inches='tight',
        dpi = 300)
    plt.show()
    plt.close()


def peptides_relatedness_analysis(
        data_directory_path,
        base_selection_cycle_number,
        top_n_peptides,
        start_sequence,
        stop_sequence,
        cdna_min_length,
        cdna_max_length,
        quality_score,
        file_name
        ):
    
    # to extract todays_date
    today = todays_date()
    
    # to collect dnas-based summary information By_cycle
    dnas_by_cycle = dnas_occurrences_by_cycle(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        quality_score,
        start_sequence,
        stop_sequence,
        )
    total_dnas_by_cycle = total_reads_by_cycle(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        quality_score,
        start_sequence,
        stop_sequence,
        )
    base_cycle_sorted_dnas = base_cycle_sorted_dnas_list(
        data_directory_path,
        base_selection_cycle_number,
        cdna_min_length,
        cdna_max_length,
        quality_score,
        start_sequence,
        stop_sequence,
        )
    dnas_appearances = dnas_appearances_by_cycle(
        base_cycle_sorted_dnas,
        dnas_by_cycle)
    
    # to collect peptides-based summary information By_cycle
    peptides_by_cycle = peptides_occurrences_by_cycle(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        quality_score,
        start_sequence,
        stop_sequence,
        )
    total_peptides_by_cycle = total_reads_by_cycle(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        quality_score,
        start_sequence,
        stop_sequence,
        )
    base_cycle_sorted_peptides = base_cycle_sorted_peptides_list(
        data_directory_path,
        base_selection_cycle_number,
        cdna_min_length,
        cdna_max_length,
        quality_score,
        start_sequence,
        stop_sequence,
        )
    peptides_appearances = peptides_appearances_by_cycle(
        base_cycle_sorted_peptides,
        peptides_by_cycle)
    
    selection_summary = complete_selection_summary(
        data_directory_path,
        start_sequence,
        stop_sequence,
        cdna_min_length,
        cdna_max_length,
        quality_score)
    sorted_cycles_list = sorted(selection_summary.keys())
    
    # to create a disjoint graph (forest),
    # based on dnas in the base_cycle
    # (joint subgraphs are the trees and the unique dna sequences are the leaves)
    base_cycle_dnas_forest = nx.Graph()
    # to add nodes (leaves, unique dna sequences)
    # to the base_cycle_dnas_forest disjoint graph
    base_cycle_dnas_forest.add_nodes_from(base_cycle_sorted_dnas)
    # to add edges
    # (twigs, dna-to-dna connections based on the hamming
    # distance between unique dna sequences) to the base_cycle_dnas_forest
    # so that disjoint graphs (stand alone trees) can be identified
    used_nodes = []
    for dna1 in base_cycle_sorted_dnas:
        used_nodes += [dna1]
        for dna2 in base_cycle_sorted_dnas:
            if dna2 not in used_nodes and hamming_distance(dna1, dna2) == 1:
                base_cycle_dnas_forest.add_edge(dna1,dna2,
                                                mutations_number = 1)
    # to extract individual joint subgraphs (stand alone trees) from the disjoint graph (forest)
    base_cycle_dnas_trees = list(
        nx.connected_component_subgraphs(base_cycle_dnas_forest, copy = True))
    
    # to create a peptideSummarydnaPerspectiveCSV file
    peptides_summary_csv = f"{today}_peptide_families_summary_{file_name}.csv"
    peptides_summary_file = open(peptides_summary_csv, 'w')
    
    # to convert list of dnas trees into a list of peptides trees leaves
    peptides_trees_leaves = []
    for dnas_tree in base_cycle_dnas_trees:
        peptide_leaves = []
        for dna in dnas_tree:
            peptide = translation(dna)
            if peptide not in peptide_leaves:
                peptide_leaves += [peptide]
        peptides_trees_leaves += [peptide_leaves]
    # to sort the resulting list of lists from the largest to smallest
    peptides_trees_leaves.sort(key = len, reverse = True)
    
    # to fix the coordinates of the origin of the graph
    positions = {}
    X_0_Coordinate = 1
    Y_0_Coordinate = 0
    Y_X0_Coordinate = 0
    
    treesXCoordinates = []
    
    peptideGraphFigure = plt.figure()
    peptideGraph = peptideGraphFigure.add_subplot(1, 1, 1)
    
    #to introduce peptideFamilyCounter
    MultiplepeptideFamilyCounter = 0
    #to introduce SinglepeptideFamilyCounter
    SinglepeptideFamilyCounter = 0
    #to introduce peptideFamilySize
    peptide_treeSize = []
    
    
    # to create a tree for each set of peptides trees leaves 
    for peptide_leaves in peptides_trees_leaves:
        
        peptide_tree = nx.Graph()
        # to convert each peptide (Leave) into a node of a peptide graph (peptideLeave on a peptide_tree)
        
        for peptide in peptide_leaves:
            peptide_tree.add_node(peptide,
                                        occurrence = peptides_by_cycle[base_selection_cycle_number][peptide],
                                        first_appearance = min(peptides_appearances[peptide]))
        # to join the peptide nodes of a graph by edges (twigs)
        for peptide1 in peptide_leaves:
            for peptide2 in peptide_leaves:
                if hamming_distance(peptide1, peptide2) == 1:
                    peptide_tree.add_edge(peptide1,peptide2)
        
        # to identify the root_peptide of a peptide_tree graph
        tree_peptides_occurrences = nx.get_node_attributes(
            peptide_tree,
            'occurrence')              
        root_peptide = max(
            tree_peptides_occurrences,
            key=tree_peptides_occurrences.get)
        
        # to create a dictionary holder for peptide and their properties (Predecessor and Occurrence)
        tree_peptides = {}
        tree_peptides[root_peptide] = [
            0,
            '',
            0,
            peptide_tree.node[root_peptide]['occurrence'],
            peptide_tree.node[root_peptide]['first_appearance']
            ]
        tree_peptides_list = list(peptide_tree.nodes())
        tree_peptides_list.remove(root_peptide)

        for peptide in tree_peptides_list:
            peptide_predecessor = nx.shortest_path(
                peptide_tree,
                source=peptide,
                target=root_peptide,
                weight=None)[1]
            # predecessor_occurrence can be used to sort the peptides, but does not seem to be useful
            predecessor_occurrence = peptide_tree.node[peptide_predecessor]['occurrence']
            peptide_occurrence = peptide_tree.node[peptide]['occurrence']

            tree_peptides[peptide] = [
                peptide_predecessor,
                predecessor_occurrence,
                peptide_occurrence]
            
        
        # to sort peptides in a peptide_tree by their distance to the root_peptide
        peptides_by_distance_to_the_root = {}
        for peptide in peptide_tree.nodes():
            distance_to_the_root = nx.shortest_path_length(
                peptide_tree,
                source=peptide,
                target=root_peptide,
                weight=None)
            if distance_to_the_root not in peptides_by_distance_to_the_root:
                peptides_by_distance_to_the_root[distance_to_the_root] = [peptide]
            else:
                peptides_by_distance_to_the_root[distance_to_the_root] += [peptide]
        
        # to identify the largest group of equidistanced peptides 
        max_peptides_number = max(
            map(lambda k:
            len(peptides_by_distance_to_the_root[k]),
            peptides_by_distance_to_the_root))

        sorted_peptides_by_distance_to_the_root = {}
        # to sort peptides by their distance to the root_peptide
        for distance_to_the_root in peptides_by_distance_to_the_root:

            equidistant_peptides = peptides_by_distance_to_the_root[
                distance_to_the_root]

            equidistant_peptides = sorted(
                equidistant_peptides,
                key=lambda peptide: (tree_peptides[peptide][2]),
                reverse=True)
            # predecessor_occurrence can be used to sort the peptides, but does not seem to be useful
            # equidistant_peptides = sorted(equidistant_peptides, key = lambda peptide: (tree_peptides[peptide][1]), reverse = True)
            equidistant_peptides = sorted(
                equidistant_peptides,
                key=lambda peptide: (tree_peptides[peptide][0]),
                reverse=False)

            additional_elements = max_peptides_number - len(equidistant_peptides)
            sorted_peptides_by_distance_to_the_root[
                distance_to_the_root] = (
                    equidistant_peptides + additional_elements * [''])

            if len(peptide_tree.nodes()) > 1:
                for peptide in equidistant_peptides:
                    XCoordinate = X_0_Coordinate + distance_to_the_root
                    YCoordinate = Y_0_Coordinate - equidistant_peptides.index(peptide)
                    positions[peptide] = (XCoordinate, YCoordinate)
                    
                                    
            elif len(peptide_tree.nodes()) == 1:
                for peptide in equidistant_peptides:
                    XCoordinate = 0
                    YCoordinate = Y_X0_Coordinate
                    positions[peptide] = (XCoordinate, YCoordinate)
                    

        #BasecyclepeptidesGraph = nx.Graph()    
        #BasecyclepeptidesGraph.add_nodes_from(base_cycle_sorted_peptides)

        sizes = []
        for peptide in peptide_tree.nodes():
            sizes.append(math.log(peptides_by_cycle[base_selection_cycle_number][peptide], 2) + 5)

        colours = []
        for peptide in peptide_tree.nodes():
            colours.append(min(peptides_appearances[peptide]))
        
        XSpan = max(map(lambda peptide: positions[peptide][0], positions)) - min(map(lambda peptide: positions[peptide][0], positions))
        YSpan = max(map(lambda peptide: positions[peptide][1], positions)) - min(map(lambda peptide: positions[peptide][1], positions))
                            
        XMin = min(map(lambda peptide: positions[peptide][0], positions)) - 0.01 * XSpan
        XMax = max(map(lambda peptide: positions[peptide][0], positions)) + 0.01 * XSpan
        YMin = min(map(lambda peptide: positions[peptide][1], positions)) - 0.02 * YSpan
        YMax = max(map(lambda peptide: positions[peptide][1], positions)) + 0.02 * YSpan
        
        
        
        number_of_colours = len(sorted_cycles_list)
        
        colour_map = plt.get_cmap('Paired', number_of_colours)
        
        nx.draw_networkx(
            peptide_tree,
            pos = positions,
            node_size = sizes,
            node_color = colours,
            cmap = colour_map,
            linewidths = 0.2,
            width = 0.2,
            with_labels = False,
            #font_size = 6,
            vmin = min(sorted_cycles_list),
            vmax = max(sorted_cycles_list))
                    
        if len(peptide_tree.nodes()) > 1:
            for distance_to_the_root in sorted_peptides_by_distance_to_the_root:
                peptides_summary_file.write(
                    f"{distance_to_the_root} mutations,frequency,rank,")
            peptides_summary_file.write("\n")

        for i in range(max_peptides_number):
            for mutations_number in sorted_peptides_by_distance_to_the_root:                        
                peptide = sorted_peptides_by_distance_to_the_root[mutations_number][i]

                if peptide != '':
                    formated_peptide = hamming_distance_based_formating(root_peptide, peptide)
                    peptide_rank = str(base_cycle_sorted_peptides.index(peptide) + 1)
                    #ClonesNumber = str(len(peptide_tree.neighbors(peptide)))
                    peptide_fraction = float(
                        (peptides_by_cycle[base_selection_cycle_number].get(peptide, 0)))/float(
                            total_peptides_by_cycle[base_selection_cycle_number])
                else:
                    formated_peptide = ''
                    #ClonesNumber = ''
                    peptide_rank = ''
                    peptide_fraction = ''

                peptides_summary_file.write(
                    f"{formated_peptide},{peptide_fraction:.2%},{peptide_rank},")
            peptides_summary_file.write("\n")
        
        
        if len(peptide_tree.nodes()) > 1:
            treesXCoordinates += [X_0_Coordinate]
            X_0_Coordinate += max(peptides_by_distance_to_the_root.keys()) + 1
            MultiplepeptideFamilyCounter += 1
            peptide_treeSize += [len(peptide_tree.nodes())]
            
            

        if len(peptide_tree.nodes()) == 1:
            Y_X0_Coordinate -= 1
            SinglepeptideFamilyCounter += 1

        peptides_summary_file.write('\n')
                    
    peptides_summary_file.close()
    
    #plt.axis('off')
    plt.axis([XMin, XMax, YMin, YMax])

    peptideLegendColour = peptideGraphFigure.add_subplot(1, 1, 1)
        
    colour_map = plt.get_cmap('Paired', number_of_colours)
    peptideLegendcolours = sorted_cycles_list
    
    LegendDotsX = XMax - 0.3 * XMax
    YIncrement = - 0.03 * YMin
    #print (YIncrement)
    
    LegendColourDotsX = np.array([LegendDotsX] * number_of_colours)
    #print (LegendColourDotsX)
    FirstYcolours = YMin + 12 * YIncrement
    #print (FirstYcolours)
    LastYcolours = YMin + (12 + number_of_colours) * YIncrement
    #print (LastYcolours)
    LegendColourDotsY = np.linspace(FirstYcolours, LastYcolours, number_of_colours, endpoint = False)
    #print (LegendColourDotsY)
    
    peptideLegendColour.scatter(x = LegendColourDotsX,
                                y = LegendColourDotsY,
                                s = 15,
                                c = peptideLegendcolours,
                                cmap = colour_map,
                                linewidths = 0.2)
    
    ColourLabels = sorted_cycles_list
#     this way of setting the colours seems to be redundant
#     ColourLabels = ['{0}'.format(i) for i in range(number_of_colours)]

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
    FirstYsizes = YMin + 5 * YIncrement
    LastYSizez = YMin + 10 * YIncrement
    LegendSizeDotsY = np.linspace(FirstYsizes, LastYSizez, 5, endpoint = False)
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
    
    
    
    #FamilySizeLabels = ['{0}'.format(i) for i in peptide_treeSize]

    #for label, x, y in zip(FamilySizeLabels, LegendSizeDotsX, LegendSizeDotsY):
    #    plt.annotate(label, xy = (x, y), xytext = (5, 0),
    #                 textcoords = 'offset points',
    #                 fontsize = 5,
    #                 ha = 'left', va = 'center')
    for i in range(len(peptide_treeSize)):
        plt.text(x = treesXCoordinates[i], y = YIncrement,
                 s = peptide_treeSize[i],
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

    peptidesSummaryfile_namePNG = str(today) + 'peptideFamiliesSummary' + file_name + '.png'
    plt.savefig(peptidesSummaryfile_namePNG, dpi = 500)
    
    
    fig = plt.gcf()
    SizeInches = fig.get_size_inches()*fig.dpi
    SizeDots = fig.get_size_inches()
    
    #print (SizeInches)
    #print (SizeDots)
    
    
    #print (XMin)
    #print (XMax)
    #print (YMin)
    #print (YMax)
    
    #print (peptide_treeSize)
    #print (len(peptide_treeSize))
    #print (treesXCoordinates)
    #print (len(treesXCoordinates))
    
    plt.show()
    plt.close()