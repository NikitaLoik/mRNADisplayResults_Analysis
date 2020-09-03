import os
import datetime

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


def get_coding_sequence(
        dna_sequence: str,
        quality_sequence: str,
        cdna_min_length: int = gp.CDNA_MIN_LENGTH,
        cdna_max_length: int = gp.CDNA_MAX_LENGTH,
        start_sequence: str = gp.START_SEQUENCE,
        stop_sequence: str = gp.STOP_SEQUENCE,
        quality_score: int = gp.QUALITY_SCORE,
        ):
    '''
    Validate the DNA coding sequence based on
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


def translate(
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


# Get the count of coding_sequences, grouped by peptide
# {peptide_Y: {coding_sequence_YZ: count_YZ}}

def get_single_cycle_data(
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
    # create an empty single_cycle_data dictionary:
    # {peptideY:    {coding_sequence_YZ:    count_YZ}}
    single_cycle_data = {}

    # go through the NGS-data file (.fastq) line by line
    # with the results from a single cycle of selection,
    # populate single_cycle_data 
    for i, line in enumerate(lines):
        # Check whether the line contains a valid sequence
        if ((start_sequence in line)
            and (stop_sequence in line)):
            coding_sequence = get_coding_sequence(
                line,  # line with DNA sequence
                lines[i + 2],  # line with quality score
                cdna_min_length,
                cdna_max_length,
                start_sequence,
                stop_sequence,
                quality_score,)                      
            if coding_sequence != None:
                # translate dna or RNA sequence into peptide
                peptide = translate(coding_sequence)
                # add sequence to a single_cycle_data
                if peptide not in single_cycle_data:
                    single_cycle_data[str(peptide)] = {str(coding_sequence): 1}
                else:
                    if coding_sequence not in single_cycle_data[str(peptide)]:
                        single_cycle_data[str(peptide)][str(coding_sequence)] = 1
                    else:
                        single_cycle_data[str(peptide)][str(coding_sequence)] += 1

    return single_cycle_data


def get_hamming_distance(
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


def format_sequence_based_on_mismatches(
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
    
    formated_sequence_2 = ''
    for i in range(len(sequence_1)):
        if sequence_1[i] == sequence_2[i]:
            formated_sequence_2 += sequence_2[i].lower()
        else:
            formated_sequence_2 += sequence_2[i]
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
    returns the count of coding sequences
    grouped by peptide, and grouped by selection cycle:
    {Selectioncycle_X:    {peptideXY:    {Codingdna_XYZ:    count_XYZ}}}
    '''

    # Create an empty selection_summary dictionary to store the results from all the cycles of selection
    complete_selection_summary = {}
    
    for file in os.listdir(data_directory_path):
        file_path = os.path.join(data_directory_path, file)
        
        # A. get the cycle number from the file name
        # (file name should have two digit number before full stop — '00.')
        # without the following condition some shit appears in the beginning of the file list  
        if file.endswith('.fastq'):
            cycle_number_first_digit = file[file.find('.')-2]
            cycle_number_second_digit = file[file.find('.')-1]
            if cycle_number_first_digit == '0':
                cycle_number = int(cycle_number_second_digit)
            elif cycle_number_first_digit != '0':
                cycle_number = int(file[file.find('.')-2 : file.find('.')])
            
            # B. get single cycle results
            single_cycle_data = get_single_cycle_data(
                file_path,
                cdna_min_length,
                cdna_max_length,
                start_sequence,
                stop_sequence,
                quality_score,)
             
            # C populate complete selection summary
            complete_selection_summary[cycle_number] = single_cycle_data
            
    return complete_selection_summary


def get_peptides_counts_by_cycle(
        data_directory_path: str,
        cdna_min_length: int = gp.CDNA_MIN_LENGTH,
        cdna_max_length: int = gp.CDNA_MAX_LENGTH,
        start_sequence: str = gp.START_SEQUENCE,
        stop_sequence: str = gp.STOP_SEQUENCE,
        quality_score: int = gp.QUALITY_SCORE,
        ):
    '''
    returns the counts of peptides groupped by cycle:
    {cycle_X:    {peptideXY:    count_XY}}
    '''
    selection_summary = get_complete_selection_summary(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        start_sequence,
        stop_sequence,
        quality_score,)
    
    peptides_counts_by_cycle = {}
    for cycle in selection_summary:
        peptides_counts_in_cycle = {}
        for peptide in selection_summary[cycle]:
            peptides_counts_in_cycle[peptide] = sum(
                selection_summary[cycle][peptide].values())
        peptides_counts_by_cycle[cycle] = peptides_counts_in_cycle
        
    return peptides_counts_by_cycle


def dnas_counts_by_cycle(
        data_directory_path: str,
        cdna_min_length: int = gp.CDNA_MIN_LENGTH,
        cdna_max_length: int = gp.CDNA_MAX_LENGTH,
        start_sequence: str = gp.START_SEQUENCE,
        stop_sequence: str = gp.STOP_SEQUENCE,
        quality_score: int = gp.QUALITY_SCORE,
        ):
    '''
    returns the counts of dnas groupped by cycle:
    {cycle_X:    {dna_XY:    count_XY}}
    '''
    selection_summary = get_complete_selection_summary(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        start_sequence,
        stop_sequence,
        quality_score,)
    
    dnas_counts_by_cycle = {}
    for cycle in selection_summary:
        dnas_counts_in_cycle = {}
        for peptide in selection_summary[cycle]:
            for dna in selection_summary[cycle][peptide]:
                dnas_counts_in_cycle[dna] = selection_summary[cycle][peptide][dna]
        dnas_counts_by_cycle[cycle] = dnas_counts_in_cycle

    return dnas_counts_by_cycle


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

    peptides_by_cycle = get_peptides_counts_by_cycle(
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


def get_base_cycle_sorted_peptides_list(
        data_directory_path: str,
        base_cycle: int,
        cdna_min_length: int = gp.CDNA_MIN_LENGTH,
        cdna_max_length: int = gp.CDNA_MAX_LENGTH,
        start_sequence: str = gp.START_SEQUENCE,
        stop_sequence: str = gp.STOP_SEQUENCE,
        quality_score: int = gp.QUALITY_SCORE,
        ):
    '''
    returns list of peptides in base cycle sorted by their count:
    [peptide_1, ..., peptide_N];
    count(peptide_1) > ... > count(peptide_N)
    '''
    peptides_by_cycle = get_peptides_counts_by_cycle(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        start_sequence,
        stop_sequence,
        quality_score,
        )
            
    peptides_counts_in_base_cycle = peptides_by_cycle[
        base_cycle]
    base_cycle_sorted_peptides_list = sorted(
        peptides_counts_in_base_cycle,
        key=peptides_counts_in_base_cycle.get,
        reverse=True)
    return base_cycle_sorted_peptides_list


def get_base_cycle_sorted_dnas_list(
        data_directory_path: str,
        base_cycle: int,
        cdna_min_length: int = gp.CDNA_MIN_LENGTH,
        cdna_max_length: int = gp.CDNA_MAX_LENGTH,
        start_sequence: str = gp.START_SEQUENCE,
        stop_sequence: str = gp.STOP_SEQUENCE,
        quality_score: int = gp.QUALITY_SCORE,
        ):
    '''
    returns list of dnas in base cycle sorted by their count:
    [dna_1, ..., dna_n]; count(dna_1) > ... > count(dna_n)
    '''
    dnas_by_cycle = dnas_counts_by_cycle(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        start_sequence,
        stop_sequence,
        quality_score,
        )
            
    dnas_counts_in_base_cycle = dnas_by_cycle[base_cycle]
    base_cycle_sorted_dnas_list = sorted(
        dnas_counts_in_base_cycle,
        key=dnas_counts_in_base_cycle.get,
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

    peptides_by_cycle = get_peptides_counts_by_cycle(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        start_sequence,
        stop_sequence,
        quality_score,
        )

    base_cycle_sorted_peptides = get_base_cycle_sorted_peptides_list(
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


def get_dna_clones_counts_by_cycle_by_peptide(
        data_directory_path: str,
        cdna_min_length: int = gp.CDNA_MIN_LENGTH,
        cdna_max_length: int = gp.CDNA_MAX_LENGTH,
        start_sequence: str = gp.START_SEQUENCE,
        stop_sequence: str = gp.STOP_SEQUENCE,
        quality_score: int = gp.QUALITY_SCORE,
        ):
    '''
    returns number of clones for each peptide groupped by cycle:
    {Selectioncycle_X:    {peptideXY:    dnaClonescounts}}
    '''
    selection_summary = get_complete_selection_summary(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        start_sequence,
        stop_sequence,
        quality_score)
    
    dna_clones_counts_by_cycle_by_peptide = {}
    for cycle in selection_summary:
        dna_clones_counts_by_peptide = {}
        for peptide in selection_summary[cycle]:
            dna_clones_counts_by_peptide[peptide] = len(
                selection_summary[cycle][peptide])
        dna_clones_counts_by_cycle_by_peptide[cycle] = dna_clones_counts_by_peptide
        
    return dna_clones_counts_by_cycle_by_peptide


def get_peptides_appearances_by_cycle(
        base_cycle_sorted_peptides_list: list,
        peptides_counts_by_cycle):
    '''
    returns for each peptide in selection a list of cycles in which this peptide appears:
    {peptide_x:    [cycle_1, ..., cycle_n]}
    '''
    peptides_appearances_by_cycle = {}
    
    for peptide in base_cycle_sorted_peptides_list:
        peptides_appearances_by_cycle[peptide] = []
        for cycle in peptides_counts_by_cycle:
            if peptide in peptides_counts_by_cycle[cycle]:
                peptides_appearances_by_cycle[peptide] += [cycle]
    return peptides_appearances_by_cycle


def dnas_appearances_by_cycle(
        base_cycle_sorted_dnas_list: list,
        dnas_counts_by_cycle):
    '''
    which returns for each dna in selection a list of cycles in which this dna appears:
    {dna_x:    [cycle_1, ..., cycle_N]}
    '''
    dnas_appearances_by_cycle = {}
    
    for dna in base_cycle_sorted_dnas_list:
        dnas_appearances_by_cycle[dna] = []
        for cycle in dnas_counts_by_cycle:
            if dna in dnas_counts_by_cycle[cycle]:
                dnas_appearances_by_cycle[dna] += [cycle]
    return dnas_appearances_by_cycle
