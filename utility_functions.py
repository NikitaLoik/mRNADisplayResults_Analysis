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

import global_parameters as gp

# =============================================================================
def get_todays_date():
    '''
    returns todays date in format YYYYMMDD
    '''
    return datetime.date.today().strftime('%Y%m%d')


def dna_coding_sequence(
        dna_sequence: str,
        quality_sequence: str,
        cdna_min_length: int = gp.CDNA_MIN_LENGTH,
        cdna_max_length: int = gp.CDNA_MAX_LENGTH,
        start_sequence: str = gp.START_SEQUENCE,
        stop_sequence: str = gp.STOP_SEQUENCE,
        quality_score: int = gp.QUALITY_SCORE,
        ):
    '''
    validates the dna coding sequence based on
    (1) its length AND
    (2) quality score of each nucleotide;
    uses ONLY ONE stop sequence;
    returns ONLY ONE coding sequence
    '''
    quality_score_string = gp.QUALITY_SCORE_STRING
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
    translation_code = gp.TRANSLATION_CODE
    transcription_code = gp.TRANSCRIPTION_CODE
    
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


# Get the occurrence of coding_sequences, grouped by peptide
# {peptide_Y: {coding_sequence_YZ: occurrence_YZ}}

def get_single_cycle_summary(
        fastq_file_path: str,
        cdna_min_length: int = gp.CDNA_MIN_LENGTH,
        cdna_max_length: int = gp.CDNA_MAX_LENGTH,
        start_sequence: str = gp.START_SEQUENCE,
        stop_sequence: str = gp.STOP_SEQUENCE,
        quality_score: int = gp.QUALITY_SCORE,
        ):
    
    with open(fastq_file_path, 'r') as f_in:
        lines = f_in.readlines()

    # to store the results from a single cycle of selection,
    # create an empty single_cycle_summary dictionary:
    # {peptideY:    {coding_sequence_YZ:    occurrence_YZ}}
    single_cycle_summary = {}

    # go through the NGS-data file (.fastq) line by line
    # with the results from a single cycle of selection,
    # populate single_cycle_summary 
    for i, line in enumerate(lines):
        # Check whether the line contains a valid sequence
        if start_sequence in line and stop_sequence in line:
            coding_sequence = dna_coding_sequence(
                line,
                lines[i + 2],
                cdna_min_length,
                cdna_max_length,
                start_sequence,
                stop_sequence,
                quality_score,)                      
            if coding_sequence != None:
                # translate dna or RNA sequence into peptide
                peptide_sequence = translation(coding_sequence)
                # add sequence to a single_cycle_summary
                if peptide_sequence not in single_cycle_summary:
                    single_cycle_summary[str(peptide_sequence)] = {str(coding_sequence): 1}
                else:
                    if coding_sequence not in single_cycle_summary[str(peptide_sequence)]:
                        single_cycle_summary[str(peptide_sequence)][str(coding_sequence)] = 1
                    else:
                        single_cycle_summary[str(peptide_sequence)][str(coding_sequence)] += 1

    return single_cycle_summary


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
        ) -> str:
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



def get_complete_selection_summary(
        data_directory_path: str,
        cdna_min_length: int = gp.CDNA_MIN_LENGTH,
        cdna_max_length: int = gp.CDNA_MAX_LENGTH,
        start_sequence: str = gp.START_SEQUENCE,
        stop_sequence: str = gp.STOP_SEQUENCE,
        quality_score: int = gp.QUALITY_SCORE,
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
            single_cycle_summary = get_single_cycle_summary(
                file_path,
                cdna_min_length,
                cdna_max_length,
                start_sequence,
                stop_sequence,
                quality_score,)
             
            # C populate complete selection summary
            complete_selection_summary[cycle_number] = single_cycle_summary
            
            # ConcatenatedResultsList-Test Print
            #print ConcatenatedResultsList
            
    return complete_selection_summary


def peptides_occurrences_by_cycle(
        data_directory_path: str,
        cdna_min_length: int = gp.CDNA_MIN_LENGTH,
        cdna_max_length: int = gp.CDNA_MAX_LENGTH,
        start_sequence: str = gp.START_SEQUENCE,
        stop_sequence: str = gp.STOP_SEQUENCE,
        quality_score: int = gp.QUALITY_SCORE,
        ):
    '''
    returns the occurrences of peptides groupped by cycle:
    {cycle_X:    {peptideXY:    occurrence_XY}}
    '''
    selection_summary = get_complete_selection_summary(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        start_sequence,
        stop_sequence,
        quality_score,)
    
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
        cdna_min_length: int = gp.CDNA_MIN_LENGTH,
        cdna_max_length: int = gp.CDNA_MAX_LENGTH,
        start_sequence: str = gp.START_SEQUENCE,
        stop_sequence: str = gp.STOP_SEQUENCE,
        quality_score: int = gp.QUALITY_SCORE,
        ):
    '''
    returns the occurrences of dnas groupped by cycle:
    {cycle_X:    {dna_XY:    occurrence_XY}}
    '''
    selection_summary = get_complete_selection_summary(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        start_sequence,
        stop_sequence,
        quality_score,)
    
    dnas_occurrences_by_cycle = {}
    for cycle in selection_summary:
        dnas_occurrences_in_cycle = {}
        for peptide in selection_summary[cycle]:
            for dna in selection_summary[cycle][peptide]:
                dnas_occurrences_in_cycle[dna] = selection_summary[cycle][peptide][dna]
        dnas_occurrences_by_cycle[cycle] = dnas_occurrences_in_cycle

    return dnas_occurrences_by_cycle


def get_total_reads_per_cycle(
        data_directory_path: str,
        cdna_min_length: int = gp.CDNA_MIN_LENGTH,
        cdna_max_length: int = gp.CDNA_MAX_LENGTH,
        start_sequence: str = gp.START_SEQUENCE,
        stop_sequence: str = gp.STOP_SEQUENCE,
        quality_score: int = gp.QUALITY_SCORE,
        ):
    '''
    returns number of reads per cycle:
    {cycle_X:    total_reads_X}
    '''
    selection_summary = get_complete_selection_summary(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        start_sequence,
        stop_sequence,
        quality_score)

    peptides_by_cycle = peptides_occurrences_by_cycle(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        start_sequence,
        stop_sequence,
        quality_score,
        )
    
    total_reads_per_cycle = {}
    for cycle in selection_summary:
        total_reads_per_cycle[cycle] = sum(peptides_by_cycle[cycle].values())
        
    return total_reads_per_cycle


def base_cycle_sorted_peptides_list(
        data_directory_path: str,
        base_cycle: int,
        cdna_min_length: int = gp.CDNA_MIN_LENGTH,
        cdna_max_length: int = gp.CDNA_MAX_LENGTH,
        start_sequence: str = gp.START_SEQUENCE,
        stop_sequence: str = gp.STOP_SEQUENCE,
        quality_score: int = gp.QUALITY_SCORE,
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
        start_sequence,
        stop_sequence,
        quality_score,
        )
            
    peptides_occurrences_in_base_cycle = peptides_by_cycle[
        base_cycle]
    base_cycle_sorted_peptides_list = sorted(
        peptides_occurrences_in_base_cycle,
        key=peptides_occurrences_in_base_cycle.get,
        reverse=True)
    return base_cycle_sorted_peptides_list


def base_cycle_sorted_dnas_list(
        data_directory_path: str,
        base_cycle: int,
        cdna_min_length: int = gp.CDNA_MIN_LENGTH,
        cdna_max_length: int = gp.CDNA_MAX_LENGTH,
        start_sequence: str = gp.START_SEQUENCE,
        stop_sequence: str = gp.STOP_SEQUENCE,
        quality_score: int = gp.QUALITY_SCORE,
        ):
    '''
    returns list of dnas in base cycle sorted by their occurrence:
    [dna_1, ..., dna_n]; occurrence(dna_1) > ... > occurrence(dna_n)
    '''
    dnas_by_cycle = dnas_occurrences_by_cycle(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        start_sequence,
        stop_sequence,
        quality_score,
        )
            
    dnas_occurrences_in_base_cycle = dnas_by_cycle[base_cycle]
    base_cycle_sorted_dnas_list = sorted(
        dnas_occurrences_in_base_cycle,
        key=dnas_occurrences_in_base_cycle.get,
        reverse=True)
    
    return base_cycle_sorted_dnas_list


def get_peptides_rank_in_base_cycle(
        data_directory_path: str,
        base_cycle: int,
        cdna_min_length: int = gp.CDNA_MIN_LENGTH,
        cdna_max_length: int = gp.CDNA_MAX_LENGTH,
        start_sequence: str = gp.START_SEQUENCE,
        stop_sequence: str = gp.STOP_SEQUENCE,
        quality_score: int = gp.QUALITY_SCORE,
        ):

    peptides_by_cycle = peptides_occurrences_by_cycle(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        start_sequence,
        stop_sequence,
        quality_score,
        )

    base_cycle_sorted_peptides = base_cycle_sorted_peptides_list(
        data_directory_path,
        base_cycle,
        cdna_min_length,
        cdna_max_length,
        start_sequence,
        stop_sequence,
        quality_score,
        )
    
    base_peptide_count = 0
    peptide_rank = 1
    
    peptides_rank_in_base_cycle = {}
    
    for peptide in base_cycle_sorted_peptides:
        peptideCount = peptides_by_cycle[base_cycle][peptide]
        if peptideCount < base_peptide_count:
            peptide_rank += 1
        
        peptides_rank_in_base_cycle[peptide] = peptide_rank
        base_peptide_count = peptideCount
        
    return peptides_rank_in_base_cycle


# occurrences is a bad word-choice here
def dna_clones_occurrences_by_cycle_by_peptide(
        data_directory_path: str,
        cdna_min_length: int = gp.CDNA_MIN_LENGTH,
        cdna_max_length: int = gp.CDNA_MAX_LENGTH,
        start_sequence: str = gp.START_SEQUENCE,
        stop_sequence: str = gp.STOP_SEQUENCE,
        quality_score: int = gp.QUALITY_SCORE,
        ):
    '''
    returns number of clones for each peptide groupped by cycle:
    {Selectioncycle_X:    {peptideXY:    dnaClonesoccurrences}}
    '''
    selection_summary = get_complete_selection_summary(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        start_sequence,
        stop_sequence,
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

# =============================================================================
# Get selection_summary_report, which, for the n_top_peptides,
# returns a summary table (.txt) and provides a summary grpaph (.png).
def get_selection_summary(
        data_directory_path: str,
        base_cycle: int,
        n_top_peptides: int,
        file_name: str,
        cdna_min_length: int = gp.CDNA_MIN_LENGTH,
        cdna_max_length: int = gp.CDNA_MAX_LENGTH,
        start_sequence: str = gp.START_SEQUENCE,
        stop_sequence: str = gp.STOP_SEQUENCE,
        quality_score: int = gp.QUALITY_SCORE,
        ):
    
    today = get_todays_date() 
    
    selection_summary_csv_name = f"{today}_selection_summary_{file_name}.csv"
    selection_summary_file = open(selection_summary_csv_name, 'w')
    
    selection_summary = get_complete_selection_summary(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        start_sequence,
        stop_sequence,
        quality_score)

    sorted_cycles_list = sorted(selection_summary.keys())
    
    peptides_by_cycle = peptides_occurrences_by_cycle(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        start_sequence,
        stop_sequence,
        quality_score,
        )

    total_peptides_by_cycle = get_total_reads_per_cycle(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        start_sequence,
        stop_sequence,
        quality_score,
        )
    
    base_cycle_sorted_peptides = base_cycle_sorted_peptides_list(
        data_directory_path,
        base_cycle,
        cdna_min_length,
        cdna_max_length,
        start_sequence,
        stop_sequence,
        quality_score,
        )
    
    base_cycle_top_sorted_peptides = base_cycle_sorted_peptides[
        0 : (n_top_peptides)]

    base_cycle_peptides_rank = get_peptides_rank_in_base_cycle(
        data_directory_path,
        base_cycle,
        cdna_min_length,
        cdna_max_length,
        start_sequence,
        stop_sequence,
        quality_score,
        )
        
    selection_summary_file.write(
        f"peptide sequence,rank (#),cdna mutants,")
    for cycle in sorted_cycles_list:
        selection_summary_file.write(
            f"C{cycle}count (#) [frequency(%)],")
    selection_summary_file.write(f"\n")
    
    for peptide in base_cycle_top_sorted_peptides:
        # base_cycle_peptide_fraction = float(
        #     (peptides_by_cycle[cycle].get(peptide, 0)))/float(
        #         total_peptides_by_cycle[base_cycle])
        peptide_rank = base_cycle_peptides_rank[peptide]
        formated_peptide = hamming_distance_based_formating(
            base_cycle_top_sorted_peptides[0],
            peptide)
        peptide_cdna_mutants = len(
            selection_summary[base_cycle][peptide])
        selection_summary_file.write(
            f"{formated_peptide},{peptide_rank},{peptide_cdna_mutants},")
            
        for cycle in sorted_cycles_list:
            peptide_fraction = float(
                (peptides_by_cycle[cycle].get(peptide, 0)))/float(
                    total_peptides_by_cycle[cycle])

            # base_fraction = peptide_fraction
            
            selection_summary_file.write(
                f"{str(peptides_by_cycle[cycle].get(peptide, 0))}"
                f" [{peptide_fraction:.1%}],")
            
            # base_fraction = peptide_fraction
        selection_summary_file.write(f"\n")
        
    selection_summary_file.write(
        f"total count (#),,")
    for cycle in sorted_cycles_list:
        selection_summary_file.write(
            f"{str(total_peptides_by_cycle[cycle])},")
    selection_summary_file.write('\n\n\n')
            
    selection_summary_file.close()
    
# =============================================================================
    # Use 'ggplot' style
    plt.style.use('fivethirtyeight')
    # Create a figure 8 x 6 inches, 300 dots per inch.
    fig, ax = plt.subplots(
        1, 1,
        figsize = (8, 6),
        dpi = 300)
    Xs = []
    Ys = []
    # Map colors onto lines  
    c_norm  = matplotlib.colors.Normalize(
        vmin = 0,
        vmax = n_top_peptides - 1)
    scalar_map = matplotlib.cm.ScalarMappable(
        norm = c_norm,
        cmap = 'gist_rainbow')
    
    peptide_labels = []
    for peptide in base_cycle_top_sorted_peptides:
    #for peptide in Top24peptidesKDs:
        peptides_fractions_by_cycle = []
        for cycle in sorted_cycles_list:
            peptides_fractions_by_cycle.append(
                float((peptides_by_cycle[cycle].get(peptide, 0)))
                / float(total_peptides_by_cycle[cycle]))
        
        x = sorted_cycles_list
        y = peptides_fractions_by_cycle
        Xs += x
        Ys += y
        # print(Ys)
        
        peptide_rank = base_cycle_peptides_rank[peptide]
        #peptide_color = scalar_map.to_rgba(peptide_rank)
        peptide_color = scalar_map.to_rgba(
            base_cycle_top_sorted_peptides.index(peptide))
        formated_peptide = hamming_distance_based_formating(
            base_cycle_top_sorted_peptides[0],
            peptide)
        
        peptide_label =  f"{formated_peptide} ({peptide_rank})"
        
        #Set peptide_label
        peptide_labels += [peptide_label]
        
        ax.plot(
            x, y,
            'o-',
            c = peptide_color,
            lw = 2.0,
            ms = 4.0,
            mew = 0.1,
            mec = '#191919')

    x_margin = 0.05 * (max(Xs) - min(Xs))
    x_min = min(Xs) - x_margin
    x_max = max(Xs) + x_margin
    y_margin = 0.05 * (max(Ys) - min(Ys))
    y_min = min(Ys) - y_margin
    y_max = max(Ys) + y_margin
    
    ax.axis([x_min, x_max, y_min, y_max])
    
    ax.tick_params(labelsize = 10)
    ax.tick_params(labelsize = 10)
    
    ax.set_xlabel(
        "Selection Cycle #",
        fontsize = 10)
    ax.set_ylabel(
        "Peptide Fraction",
        fontsize = 10)
    
    legend = plt.legend(
        peptide_labels,
        title = 'cyclic-peptide random region',
        loc = 'upper center',
        bbox_to_anchor = (0.5, -0.10),
        fancybox = True,
        shadow = False,
        fontsize = 10,
        ncol = 3)
    
    ax.get_legend().get_title().set_size('small')
    
    selection_summary_png_path = f"{today}_selection_summary_{file_name}.png"
    
    fig.savefig(
        selection_summary_png_path,
        bbox_extra_artists = [legend],
        bbox_inches = 'tight',
        dpi = 300)
    plt.show()
    plt.close()


# =============================================================================
def dna_mutants_analysis(
        data_directory_path: str,
        base_cycle: int,
        n_top_peptides: int,
        file_name: str,
        cdna_min_length: int = gp.CDNA_MIN_LENGTH,
        cdna_max_length: int = gp.CDNA_MAX_LENGTH,
        start_sequence: str = gp.START_SEQUENCE,
        stop_sequence: str = gp.STOP_SEQUENCE,
        quality_score: int = gp.QUALITY_SCORE,
        ):
        
    
    today = get_todays_date() 
    
    dna_mutants_analysis_csv =  f"{today}_dnas_mutants_analysis_{file_name}.csv"
    dna_mutants_analysis_file = open(dna_mutants_analysis_csv, 'w')
    
    selection_summary = get_complete_selection_summary(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        start_sequence,
        stop_sequence,
        quality_score)
    sorted_cycles_list = sorted(selection_summary.keys())
    
    peptides_by_cycle = peptides_occurrences_by_cycle(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        start_sequence,
        stop_sequence,
        quality_score,
        )
    # total_peptides_by_cycle = get_total_reads_per_cycle(
    #     data_directory_path,
    #     cdna_min_length,
    #     cdna_max_length,
    #     start_sequence,
    #     stop_sequence,
    #     quality_score,
    #     )
    
    base_cycle_sorted_peptides = base_cycle_sorted_peptides_list(
        data_directory_path,
        base_cycle,
        cdna_min_length,
        cdna_max_length,
        start_sequence,
        stop_sequence,
        quality_score,
        )
    base_cycle_top_sorted_peptides = base_cycle_sorted_peptides[0 : (n_top_peptides)]
    
    dna_clones_by_cycle_by_peptide = dna_clones_occurrences_by_cycle_by_peptide(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        start_sequence,
        stop_sequence,
        quality_score,
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
    
# =============================================================================        

    # # Create a figure of size 8x6 inches, 500 dots per inch
    # plt.figure(figsize = (8, 6),
    #            dpi = 500)
    # # Create 'ggplot' style
    # plt.style.use('fivethirtyeight')
    # # Create a new subplot from a grid of 1x1
    # Graph = plt.subplot(1, 1, 1)
    # Use 'ggplot' style
    plt.style.use('fivethirtyeight')
    # Create a figure 8 x 6 inches, 300 dots per inch.
    fig, ax = plt.subplots(
        1, 1,
        figsize = (8, 6),
        dpi = 300)
    
#    peptide_dna_clones_number_in_base_cycle = []
#    peptide_occurrence_in_base_cycle = []

    # Map colors onto lines
    c_norm  = matplotlib.colors.Normalize(
        vmin = 0,
        vmax = len(base_cycle_sorted_peptides) - 1)
    scalar_map = matplotlib.cm.ScalarMappable(
        norm = c_norm,
        cmap = 'gist_rainbow')

    cycle_index = base_cycle

    Xs = []
    Ys = []        
    for peptide in dna_clones_by_cycle_by_peptide[cycle_index]:
        peptide_dna_clones_number_in_base_cycle = math.log(
            dna_clones_by_cycle_by_peptide[cycle_index].get(peptide, 0), 2)
        peptide_occurrence_in_base_cycle = math.log(
            peptides_by_cycle[cycle_index].get(peptide, 0), 2)
        
        peptide_color = scalar_map.to_rgba(
            base_cycle_sorted_peptides.index(peptide))
    
        x = peptide_dna_clones_number_in_base_cycle
        y = peptide_occurrence_in_base_cycle
        Xs += [x]
        Ys += [y]
        
        plt.plot(x, y,
                'o',
                c = peptide_color,
                ms = 5.0,
                mew = 0.1,
                mec = '#191919')

    x_min = min(Xs) - 0.05*(max(Xs) - min(Xs))
    x_max = max(Xs) + 0.05*(max(Xs) - min(Xs))
    y_min = min(Ys) - 0.05*(max(Ys) - min(Ys))
    y_max = max(Ys) + 0.05*(max(Ys) - min(Ys))
    
    ax.axis([x_min, x_max, y_min, y_max])
        
    # x_label = 
    ax.set_xlabel(
        'log$2$ (# dna clones)',  # $_$ makes subscript possible
        fontsize = 14)
    # YLabel = 
    ax.set_ylabel(
        'log$2$ (peptide occurrence)',  # $_$ makes subscript possible
        fontsize = 14)
    
    legend = plt.legend(base_cycle_sorted_peptides,
                        loc = 'upper center',
                        bbox_to_anchor = (0.5, -0.15),
                        fancybox = True,
                        shadow = False,
                        ncol = 4)
    
    dna_clones_analysis_png_path = (
        f"{today}_dnas_mutants_analysis_regression_C{cycle}_{file_name}.png")
    
    fig.savefig(
        dna_clones_analysis_png_path,
        bbox_extra_artists=[legend],
        bbox_inches='tight',
        dpi = 300)
    plt.show()
    plt.close()

# =============================================================================
def peptides_relatedness_analysis(
        data_directory_path,
        base_cycle,
        n_top_peptides,
        file_name,
        cdna_min_length: int = gp.CDNA_MIN_LENGTH,
        cdna_max_length: int = gp.CDNA_MAX_LENGTH,
        start_sequence: str = gp.START_SEQUENCE,
        stop_sequence: str = gp.STOP_SEQUENCE,
        quality_score: int = gp.QUALITY_SCORE,
        ):
    
    # Get todays date/
    today = get_todays_date()
    
    # Get DNAs-based summary by cycle
    # dnas_by_cycle = dnas_occurrences_by_cycle(
    #     data_directory_path,
    #     cdna_min_length,
    #     cdna_max_length,
    #     start_sequence,
    #     stop_sequence,
    #     quality_score,
    #     )
    # total_dnas_by_cycle = get_total_reads_per_cycle(
    #     data_directory_path,
    #     cdna_min_length,
    #     cdna_max_length,
    #     start_sequence,
    #     stop_sequence,
    #     quality_score,
    #     )
    base_cycle_sorted_dnas = base_cycle_sorted_dnas_list(
        data_directory_path,
        base_cycle,
        cdna_min_length,
        cdna_max_length,
        start_sequence,
        stop_sequence,
        quality_score,
        )
    # dnas_appearances = dnas_appearances_by_cycle(
    #     base_cycle_sorted_dnas,
    #     dnas_by_cycle)
    
    # Get peptides-based summary by cycle
    peptides_by_cycle = peptides_occurrences_by_cycle(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        start_sequence,
        stop_sequence,
        quality_score,
        )
    total_peptides_by_cycle = get_total_reads_per_cycle(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        start_sequence,
        stop_sequence,
        quality_score,
        )
    base_cycle_sorted_peptides = base_cycle_sorted_peptides_list(
        data_directory_path,
        base_cycle,
        cdna_min_length,
        cdna_max_length,
        start_sequence,
        stop_sequence,
        quality_score,
        )

    peptides_appearances = peptides_appearances_by_cycle(
        base_cycle_sorted_peptides,
        peptides_by_cycle)
    
    selection_summary = get_complete_selection_summary(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        start_sequence,
        stop_sequence,
        quality_score)

    sorted_cycles_list = sorted(selection_summary.keys())
    
    # Get a disjoint graph (forest),
    # based on DNAs in the base cycle
    # (joint subgraphs are the trees and the unique DNA-sequences are the leaves)
    base_cycle_dnas_forest = nx.Graph()
    # to add nodes (leaves, unique dna sequences)
    # to the base_cycle_dnas_forest disjoint graph
    base_cycle_dnas_forest.add_nodes_from(base_cycle_sorted_dnas)
    # Add edges between DNA sequences
    # (if hamming distance between two DNA sequences = 1)
    # to the base_cycle_dnas_forest
    # so that disjoint graphs (stand alone trees) can be identified
    used_nodes = []
    for dna1 in base_cycle_sorted_dnas:
        used_nodes += [dna1]
        for dna2 in base_cycle_sorted_dnas:
            if ((dna2 not in used_nodes)
                and (hamming_distance(dna1, dna2) == 1)):
                base_cycle_dnas_forest.add_edge(
                    dna1,
                    dna2,
                    n_mutations = 1)
    # Get individual joint subgraphs (stand alone trees)
    # from the disjoint graph (forest).
    base_cycle_dnas_trees = list(
        nx.connected_component_subgraphs(
            base_cycle_dnas_forest,
            copy=True))
    
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
    peptides_trees_leaves.sort(
        key=len,
        reverse=True)
    
    # Fix the coordinates of the origin of the graph.
    positions = {}
    x_0_coordinate = 1
    y_0_coordinate = 0
    y_x_0_coordinate = 0
    
    trees_x_coordinates = []
    
    # Introduce peptideFamilyCounter
    many_peptides_families_counter = 0
    # Introduce one_peptide_families_counter
    one_peptide_families_counter = 0
    # Introduce peptideFamilysize
    peptide_tree_size = []
    
    # PLOT SETUP START ========================================================
    # Set color map.
    n_colors = len(sorted_cycles_list)
    color_map = plt.get_cmap('Paired', n_colors)
    # Set figure and axes.
    fig, ax = plt.subplots(
        1, 1,
        # figsize = (8, 6),
        dpi = 300)
    # PLOT SETUP END ==========================================================

    # Get a tree for each set of peptides trees leaves.
    for peptide_leaves in peptides_trees_leaves:
        peptide_tree = nx.Graph()
        # to convert each peptide (Leave) into a node of a peptide graph (peptideLeave on a peptide_tree)
        
        for peptide in peptide_leaves:
            peptide_tree.add_node(
                peptide,
                occurrence = peptides_by_cycle[base_cycle][peptide],
                first_appearance = min(peptides_appearances[peptide]))
        # Join peptide nodes by edges, if hamming distance = 1
        for peptide1 in peptide_leaves:
            for peptide2 in peptide_leaves:
                if hamming_distance(peptide1, peptide2) == 1:
                    peptide_tree.add_edge(peptide1,peptide2)
        
        # Get root-peptide of a peptide-tree.
        tree_peptides_occurrences = nx.get_node_attributes(
            peptide_tree,
            'occurrence')              
        root_peptide = max(
            tree_peptides_occurrences,
            key=tree_peptides_occurrences.get)
        
        # Make a dictionary holder for peptides and their properties
        # (predecessor and occurrence)
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
            # Predecessor occurrence can be used to sort the peptides.
            # However, it does not seem to be useful.
            predecessor_occurrence = peptide_tree.node[peptide_predecessor]['occurrence']
            peptide_occurrence = peptide_tree.node[peptide]['occurrence']

            tree_peptides[peptide] = [
                peptide_predecessor,
                predecessor_occurrence,
                peptide_occurrence]
            
        # Sort peptides in a peptide_tree by their distance to the root_peptide.
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
        
        # Find the largest group of equidistanced peptides.
        n_peptides_max = max(
            map(lambda k:
            len(peptides_by_distance_to_the_root[k]),
            peptides_by_distance_to_the_root))

        sorted_peptides_by_distance_to_the_root = {}
        # Sort peptides by their distance to the root_peptide.
        for distance_to_the_root in peptides_by_distance_to_the_root:

            equidistant_peptides = peptides_by_distance_to_the_root[
                distance_to_the_root]

            equidistant_peptides = sorted(
                equidistant_peptides,
                key=lambda peptide: (tree_peptides[peptide][2]),
                reverse=True)
            # Predecessor occurrence can be used to sort the peptides.
            # However, it does not seem to be useful.
            # equidistant_peptides = sorted(
            #     equidistant_peptides,
            #     key = lambda peptide: (tree_peptides[peptide][1]),
            #     reverse = True)
            equidistant_peptides = sorted(
                equidistant_peptides,
                key=lambda peptide: (tree_peptides[peptide][0]),
                reverse=False)

            additional_elements = n_peptides_max - len(equidistant_peptides)
            sorted_peptides_by_distance_to_the_root[
                distance_to_the_root] = (
                    equidistant_peptides + additional_elements * [''])

            if len(peptide_tree.nodes()) > 1:
                for peptide in equidistant_peptides:
                    XCoordinate = x_0_coordinate + distance_to_the_root
                    YCoordinate = y_0_coordinate - equidistant_peptides.index(peptide)
                    positions[peptide] = (XCoordinate, YCoordinate)
                    
                                    
            elif len(peptide_tree.nodes()) == 1:
                for peptide in equidistant_peptides:
                    XCoordinate = 0
                    YCoordinate = y_x_0_coordinate
                    positions[peptide] = (XCoordinate, YCoordinate)
                    
        # Get marker size, proportional to peptides occurence in a base cycle.
        sizes = []
        for peptide in peptide_tree.nodes():
            sizes.append(
                math.log(peptides_by_cycle[base_cycle][peptide], 2) + 5)
        # Get marker color based on the peptides first appearance.
        colors = []
        for peptide in peptide_tree.nodes():
            colors.append(min(peptides_appearances[peptide]))
        
        x_span = (
            max(map(lambda peptide: positions[peptide][0], positions))
            - min(map(lambda peptide: positions[peptide][0], positions)))
        y_span = (
            max(map(lambda peptide: positions[peptide][1], positions))
            - min(map(lambda peptide: positions[peptide][1], positions)))
                            
        x_min = min(map(lambda peptide: positions[peptide][0], positions)) - 0.01 * x_span
        x_max = max(map(lambda peptide: positions[peptide][0], positions)) + 0.01 * x_span
        y_min = min(map(lambda peptide: positions[peptide][1], positions)) - 0.02 * y_span
        y_max = max(map(lambda peptide: positions[peptide][1], positions)) + 0.02 * y_span
        
        # Plot peptide tree.
        nx.draw_networkx(
            peptide_tree,
            pos = positions,
            node_size = sizes,
            node_color = colors,
            cmap = color_map,
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

        for i in range(n_peptides_max):
            for n_mutations in sorted_peptides_by_distance_to_the_root:                        
                peptide = sorted_peptides_by_distance_to_the_root[n_mutations][i]

                if peptide != '':
                    formated_peptide = hamming_distance_based_formating(
                        root_peptide,
                        peptide)
                    peptide_rank = str(base_cycle_sorted_peptides.index(peptide) + 1)
                    #n_clones = str(len(peptide_tree.neighbors(peptide)))
                    peptide_fraction = (float(
                        (peptides_by_cycle[base_cycle].get(peptide, 0)))
                        / float(
                            total_peptides_by_cycle[base_cycle]))
                else:
                    formated_peptide = ''
                    #n_clones = 0
                    peptide_rank = ''
                    peptide_fraction = 0.

                peptides_summary_file.write(
                    f"{formated_peptide},{peptide_fraction:.1%},{peptide_rank},")
            peptides_summary_file.write("\n")
        
        
        if len(peptide_tree.nodes()) > 1:
            trees_x_coordinates += [x_0_coordinate]
            x_0_coordinate += max(peptides_by_distance_to_the_root.keys()) + 1
            many_peptides_families_counter += 1
            peptide_tree_size += [len(peptide_tree.nodes())]
            
        if len(peptide_tree.nodes()) == 1:
            y_x_0_coordinate -= 1
            one_peptide_families_counter += 1

        peptides_summary_file.write("\n")
                    
    peptides_summary_file.close()


# =============================================================================
    ax.axis([x_min, x_max, y_min, y_max])
    
    legend_dots_x = x_max - 0.3 * x_max
    y_increment = - 0.03 * y_min
    #print (y_increment)
    
    legend_color_dots_X = np.array([legend_dots_x] * n_colors)
    #print (legend_color_dots_X)
    first_y_colors = y_min + 12 * y_increment
    #print (first_y_colors)
    last_y_colors = y_min + (12 + n_colors) * y_increment
    #print (last_y_colors)
    legend_color_dots_Y = np.linspace(
        first_y_colors,
        last_y_colors,
        n_colors,
        endpoint=False)
    #print (legend_color_dots_Y)
    peptide_legend_colors = sorted_cycles_list
    
    ax.scatter(
        x = legend_color_dots_X,
        y = legend_color_dots_Y,
        s = 15,
        c = peptide_legend_colors,
        cmap = color_map,
        linewidths = 0.2)
    
    color_labels = sorted_cycles_list
#     this way of setting the colors seems to be redundant
#     color_labels = ['{0}'.format(i) for i in range(n_colors)]

    for label, x, y in zip(color_labels, legend_color_dots_X, legend_color_dots_Y):
        ax.annotate(
            label,
            xy = (x, y),
            xytext = (5, 0),
            textcoords = 'offset points',
            fontsize = 5,
            ha = 'left', va = 'center')
    ax.text(
        x = legend_dots_x,
        y = (max(legend_color_dots_Y) + y_increment),
        s = 'first-appearance cycle #',
        fontsize = 5)

    size = []
    for i in [1, 10, 100, 1000, 10000]:
        size.append(math.log(i, 2) + 5)

    legend_size_dots_x = np.array([legend_dots_x] * 5)
    first_y_sizes = y_min + 5 * y_increment
    last_y_sizes = y_min + 10 * y_increment
    legend_size_dots_y = np.linspace(
        first_y_sizes,
        last_y_sizes,
        5,
        endpoint = False)
    ax.scatter(
        x = legend_size_dots_x,
        y = legend_size_dots_y,
        s = size,
        c = 'w',
        linewidths = 0.2)

    size_labels = [f"{i}" for i in [1, 10, 100, 1000, 10000]]

    for label, x, y in zip(size_labels, legend_size_dots_x, legend_size_dots_y):
        ax.annotate(
            label,
            xy = (x, y),
            xytext = (5, 0),
            textcoords = 'offset points',
            fontsize = 5,
            ha = 'left',
            va = 'center')
    ax.text(
        x = legend_dots_x,
        y = (max(legend_size_dots_y) - 0.03 * y_min),
        s = "frequency in the last cycle",
        fontsize = 5)
    
    
    
    #Familysize_labels = ['{0}'.format(i) for i in peptide_tree_size]

    #for label, x, y in zip(Familysize_labels, legend_size_dots_x, legend_size_dots_y):
    #    plt.annotate(label, xy = (x, y), xytext = (5, 0),
    #                 textcoords = 'offset points',
    #                 fontsize = 5,
    #                 ha = 'left', va = 'center')
    for i in range(len(peptide_tree_size)):
        ax.text(x = trees_x_coordinates[i], y = y_increment,
                 s = peptide_tree_size[i],
                 fontsize = 5)
    
    ax.text(
        x=legend_dots_x,
        y=y_min + 3 * y_increment,
        s = f"total # unique peptide sequence {len(base_cycle_sorted_peptides)}",
        fontsize = 5)
    ax.text(
        x = legend_dots_x,
        y = y_min + 2 * y_increment,
        s = f"single-member peptide family # {one_peptide_families_counter}",
        fontsize = 5)
    ax.text(
        x = legend_dots_x,
        y = y_min + 1 * y_increment,
        s = f"multiple-member peptide family # {many_peptides_families_counter}",
        fontsize = 5)
    
    ax.axis('off')

    peptides_summary_png = f"{today}_peptide_families_summary_{file_name}.png"
    fig.savefig(peptides_summary_png, dpi = 500)
    
    
    # fig = plt.gcf()
    # size_inches = fig.get_size_inches()*fig.dpi
    # size_dots = fig.get_size_inches()
    
    #print (size_inches)
    #print (size_dots)
    
    #print (x_min)
    #print (x_max)
    #print (y_min)
    #print (y_max)
    
    #print (peptide_tree_size)
    #print (len(peptide_tree_size))
    #print (trees_x_coordinates)
    #print (len(trees_x_coordinates))
    
    plt.show()
    plt.close()