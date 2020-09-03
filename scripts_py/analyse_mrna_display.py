import os, sys, inspect
import datetime

import argparse

import math
import numpy as np
import networkx as nx

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'monospace'

working_dir = os.path.dirname(
    os.path.abspath(inspect.getfile(inspect.currentframe())))
parent_dir = os.path.dirname(working_dir)
sys.path.insert(0, parent_dir)

from scripts_py import global_parameters as gp
from scripts_py import utility_functions as uf

# SET LOGGER ==================================================================
import logging
logging.basicConfig(
    level=logging.INFO,
    format='%(levelname)s: %(asctime)s: %(filename)s: %(lineno)s:\n%(message)s')
logger = logging.getLogger(__name__)


# DEFINE SELECTION CLASS ======================================================
class Selection:
    def __init__(
            self,
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
        self.data_directory_path = data_directory_path
        self.base_cycle = base_cycle
        self.n_top_peptides = n_top_peptides
        self.file_name = file_name
        self.cdna_min_length = cdna_min_length
        self.cdna_max_length = cdna_max_length
        self.start_sequence = start_sequence
        self.stop_sequence = stop_sequence
        self.quality_score = quality_score

        self.selection_summary = uf.get_complete_selection_summary(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        start_sequence,
        stop_sequence,
        quality_score)

        self.sorted_cycles_list = sorted(self.selection_summary.keys())
        
        self.peptides_by_cycle = uf.get_peptides_counts_by_cycle(
            data_directory_path,
            cdna_min_length,
            cdna_max_length,
            start_sequence,
            stop_sequence,
            quality_score,
            )

        self.total_peptides_by_cycle = uf.get_total_reads_per_cycle(
            data_directory_path,
            cdna_min_length,
            cdna_max_length,
            start_sequence,
            stop_sequence,
            quality_score,
            )
        
        self.base_cycle_sorted_peptides = uf.get_base_cycle_sorted_peptides_list(
            data_directory_path,
            base_cycle,
            cdna_min_length,
            cdna_max_length,
            start_sequence,
            stop_sequence,
            quality_score,
            )
        
        self.base_cycle_top_sorted_peptides = self.base_cycle_sorted_peptides[
            0 : (n_top_peptides)]

        self.base_cycle_peptides_rank = uf.get_peptides_rank_in_base_cycle(
            data_directory_path,
            base_cycle,
            cdna_min_length,
            cdna_max_length,
            start_sequence,
            stop_sequence,
            quality_score,
            )

        self.dna_clones_by_cycle_by_peptide = uf.get_dna_clones_counts_by_cycle_by_peptide(
        data_directory_path,
        cdna_min_length,
        cdna_max_length,
        start_sequence,
        stop_sequence,
        quality_score,
        )

        self.base_cycle_sorted_dnas = uf.get_base_cycle_sorted_dnas_list(
        data_directory_path,
        base_cycle,
        cdna_min_length,
        cdna_max_length,
        start_sequence,
        stop_sequence,
        quality_score,
        )
        self.peptides_appearances = uf.get_peptides_appearances_by_cycle(
            self.base_cycle_sorted_peptides,
            self.peptides_by_cycle)


    def get_summary(self):
    # Get selection_summary_report, which, for the n_top_peptides,
    # returns a summary table (.txt) and provides a summary grpaph (.png).
        # A. Get summary .csv and .png path.
        today = uf.get_todays_date()
        selection_summary_csv_path = f"{today}_selection_summary_{self.file_name}.csv"
        selection_summary_png_path = f"{today}_selection_summary_{self.file_name}.png"

        # B. Populate summary .csv.
        summary_file = open(selection_summary_csv_path, 'w')
        summary_file.write(
            f"peptide sequence,rank (#),cdna mutants,")
        for cycle in self.sorted_cycles_list:
            summary_file.write(
                f"C{cycle}count (#) [frequency(%)],")
        summary_file.write(f"\n")
        
        for peptide in self.base_cycle_top_sorted_peptides:
            peptide_rank = self.base_cycle_peptides_rank[peptide]
            formated_peptide = uf.format_sequence_based_on_mismatches(
                self.base_cycle_top_sorted_peptides[0],
                peptide)
            peptide_cdna_mutants = len(
                self.selection_summary[self.base_cycle][peptide])
            summary_file.write(
                f"{formated_peptide},{peptide_rank},{peptide_cdna_mutants},")
                
            for cycle in self.sorted_cycles_list:
                peptide_fraction = float(
                    (self.peptides_by_cycle[cycle].get(peptide, 0)))/float(
                        self.total_peptides_by_cycle[cycle])
                
                summary_file.write(
                    f"{str(self.peptides_by_cycle[cycle].get(peptide, 0))}"
                    f" [{peptide_fraction:.1%}],")
                
            summary_file.write(f"\n")
            
        summary_file.write(
            f"total count (#),,")
        for cycle in self.sorted_cycles_list:
            summary_file.write(
                f"{str(self.total_peptides_by_cycle[cycle])},")
        summary_file.write('\n\n\n')
                
        summary_file.close()
        
        # C. Get summary .png plot.
        # Use 'ggplot' style
        plt.style.use('fivethirtyeight')
        fig, ax = plt.subplots(
            1, 1,
            figsize = (8, 6),
            dpi = 300)
        Xs = []
        Ys = []
        # Map colors onto lines  
        c_norm  = matplotlib.colors.Normalize(
            vmin = 0,
            vmax = self.n_top_peptides - 1)
        scalar_map = matplotlib.cm.ScalarMappable(
            norm = c_norm,
            cmap = 'gist_rainbow')

        # Set peptide_labels
        peptide_labels = []
        for peptide in self.base_cycle_top_sorted_peptides:
            peptides_fractions_by_cycle = []
            for cycle in self.sorted_cycles_list:
                peptides_fractions_by_cycle.append(
                    float((self.peptides_by_cycle[cycle].get(peptide, 0)))
                    / float(self.total_peptides_by_cycle[cycle]))
            
            x = self.sorted_cycles_list
            y = peptides_fractions_by_cycle
            Xs += x
            Ys += y
            
            peptide_rank = self.base_cycle_peptides_rank[peptide]
            peptide_color = scalar_map.to_rgba(
                self.base_cycle_top_sorted_peptides.index(peptide))
            formated_peptide = uf.format_sequence_based_on_mismatches(
                self.base_cycle_top_sorted_peptides[0],
                peptide)
            
            peptide_label =  f"{formated_peptide} ({peptide_rank})"
            
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
        # ax.tick_params(labelsize = 10)
        
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
        
        fig.savefig(
            selection_summary_png_path,
            bbox_extra_artists = [legend],
            bbox_inches = 'tight',
            dpi = 300)
        plt.show()
        plt.close()


# =============================================================================
    def get_dna_mutants_summary(self):
        # A. Get mutants .csv and .png path.
        today = uf.get_todays_date()
        dna_mutants_analysis_csv =  f"{today}_dnas_mutants_analysis_{self.file_name}.csv"
        dna_clones_analysis_png_path = (
            f"{today}_dnas_mutants_analysis_regression_C{self.base_cycle}_{self.file_name}.png")
        
        # B. Populate mutants .csv.
        dna_mutants_analysis_file = open(dna_mutants_analysis_csv, 'w')
        dna_mutants_analysis_file.write("peptide sequence,")
        for cycle in self.sorted_cycles_list:
            dna_mutants_analysis_file.write(f"cycle # {cycle} dna clones (#),")
        dna_mutants_analysis_file.write("\n")
        
        for peptide in self.base_cycle_top_sorted_peptides:
            dna_mutants_analysis_file.write(f"{peptide}," )
            for cycle in self.sorted_cycles_list:
                dna_mutants_analysis_file.write(
                    f"{str(self.dna_clones_by_cycle_by_peptide[cycle].get(peptide, 0)),}")
            dna_mutants_analysis_file.write("\n")
        dna_mutants_analysis_file.close()
        
    # =============================================================================        
        # C. Get mutants .png plot.
        plt.style.use('fivethirtyeight')
        # Create a figure 8 x 6 inches, 300 dots per inch.
        fig, ax = plt.subplots(
            1, 1,
            figsize = (8, 6),
            dpi = 300)

        # Map colors onto lines
        c_norm  = matplotlib.colors.Normalize(
            vmin = 0,
            vmax = len(self.base_cycle_sorted_peptides) - 1)
        scalar_map = matplotlib.cm.ScalarMappable(
            norm = c_norm,
            cmap = 'gist_rainbow')

        cycle_index = self.base_cycle

        Xs = []
        Ys = []        
        for peptide in self.dna_clones_by_cycle_by_peptide[cycle_index]:
            peptide_dna_clones_number_in_base_cycle = math.log(
                self.dna_clones_by_cycle_by_peptide[cycle_index].get(peptide, 0), 2)
            peptide_count_in_base_cycle = math.log(
                self.peptides_by_cycle[cycle_index].get(peptide, 0), 2)
            
            peptide_color = scalar_map.to_rgba(
                self.base_cycle_sorted_peptides.index(peptide))
        
            x = peptide_dna_clones_number_in_base_cycle
            y = peptide_count_in_base_cycle
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
            
        ax.set_xlabel(
            'log$_2$ (# dna clones)',  # $_X$ makes subscript possible
            fontsize = 14)
        ax.set_ylabel(
            'log$_2$ (peptide count)',  # $_X$ makes subscript possible
            fontsize = 14)
        
        legend = plt.legend(self.base_cycle_sorted_peptides,
                            loc = 'upper center',
                            bbox_to_anchor = (0.5, -0.15),
                            fancybox = True,
                            shadow = False,
                            ncol = 4)
        
        fig.savefig(
            dna_clones_analysis_png_path,
            bbox_extra_artists=[legend],
            bbox_inches='tight',
            dpi = 300)
        plt.show()
        plt.close()


# =============================================================================
    def get_peptides_relatedness_summary(self):
        # Get a disjoint graph (forest),
        # based on DNAs in the base cycle
        # (joint subgraphs are the trees and the unique DNA-sequences are the leaves)
        base_cycle_dnas_forest = nx.Graph()
        # to add nodes (leaves, unique dna sequences)
        # to the base_cycle_dnas_forest disjoint graph
        base_cycle_dnas_forest.add_nodes_from(self.base_cycle_sorted_dnas)
        # Add edges between DNA sequences
        # (if hamming distance between two DNA sequences = 1)
        # to the base_cycle_dnas_forest
        # so that disjoint graphs (stand alone trees) can be identified
        used_nodes = []
        for dna1 in self.base_cycle_sorted_dnas:
            used_nodes += [dna1]
            for dna2 in self.base_cycle_sorted_dnas:
                if ((dna2 not in used_nodes)
                    and (uf.get_hamming_distance(dna1, dna2) == 1)):
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
        
        # to convert list of dnas trees into a list of peptides trees leaves
        peptides_trees_leaves = []
        for dnas_tree in base_cycle_dnas_trees:
            peptide_leaves = []
            for dna in dnas_tree:
                peptide = uf.translate(dna)
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
        

        # B. Get mutants .csv and .png path
        today = uf.get_todays_date()
        peptides_summary_csv = f"{today}_peptide_families_summary_{self.file_name}.csv"
        peptides_summary_png = f"{today}_peptide_families_summary_{self.file_name}.png"
        peptides_summary_file = open(peptides_summary_csv, 'w')

        # PLOT SETUP START ========================================================
        # Set color map.
        n_colors = len(self.sorted_cycles_list)
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
                    count = self.peptides_by_cycle[self.base_cycle][peptide],
                    first_appearance = min(self.peptides_appearances[peptide]))
            # Join peptide nodes by edges, if hamming distance = 1
            for peptide1 in peptide_leaves:
                for peptide2 in peptide_leaves:
                    if uf.get_hamming_distance(peptide1, peptide2) == 1:
                        peptide_tree.add_edge(peptide1,peptide2)
            
            # Get root-peptide of a peptide-tree.
            tree_peptides_counts = nx.get_node_attributes(
                peptide_tree,
                'count')              
            root_peptide = max(
                tree_peptides_counts,
                key=tree_peptides_counts.get)
            
            # Make a dictionary holder for peptides and their properties
            # (predecessor and count)
            tree_peptides = {}
            tree_peptides[root_peptide] = [
                0,
                '',
                0,
                peptide_tree.node[root_peptide]['count'],
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
                # Predecessor count can be used to sort the peptides.
                # However, it does not seem to be useful.
                predecessor_count = peptide_tree.node[peptide_predecessor]['count']
                peptide_count = peptide_tree.node[peptide]['count']

                tree_peptides[peptide] = [
                    peptide_predecessor,
                    predecessor_count,
                    peptide_count]
                
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
                # Predecessor count can be used to sort the peptides.
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
                    math.log(self.peptides_by_cycle[self.base_cycle][peptide], 2) + 5)
            # Get marker color based on the peptides first appearance.
            colors = []
            for peptide in peptide_tree.nodes():
                colors.append(min(self.peptides_appearances[peptide]))
            
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
                vmin = min(self.sorted_cycles_list),
                vmax = max(self.sorted_cycles_list))
                        
            if len(peptide_tree.nodes()) > 1:
                for distance_to_the_root in sorted_peptides_by_distance_to_the_root:
                    peptides_summary_file.write(
                        f"{distance_to_the_root} mutations,frequency,rank,")
                peptides_summary_file.write("\n")

            for i in range(n_peptides_max):
                for n_mutations in sorted_peptides_by_distance_to_the_root:                        
                    peptide = sorted_peptides_by_distance_to_the_root[n_mutations][i]

                    if peptide != '':
                        formated_peptide = uf.format_sequence_based_on_mismatches(
                            root_peptide,
                            peptide)
                        peptide_rank = str(self.base_cycle_sorted_peptides.index(peptide) + 1)
                        #n_clones = str(len(peptide_tree.neighbors(peptide)))
                        peptide_fraction = (float(
                            (self.peptides_by_cycle[self.base_cycle].get(peptide, 0)))
                            / float(
                                self.total_peptides_by_cycle[self.base_cycle]))
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
        peptide_legend_colors = self.sorted_cycles_list
        
        ax.scatter(
            x = legend_color_dots_X,
            y = legend_color_dots_Y,
            s = 15,
            c = peptide_legend_colors,
            cmap = color_map,
            linewidths = 0.2)
        
        color_labels = self.sorted_cycles_list
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
        
        for i in range(len(peptide_tree_size)):
            ax.text(x = trees_x_coordinates[i], y = y_increment,
                    s = peptide_tree_size[i],
                    fontsize = 5)
        
        ax.text(
            x=legend_dots_x,
            y=y_min + 3 * y_increment,
            s = f"total # unique peptide sequence {len(self.base_cycle_sorted_peptides)}",
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

        fig.savefig(peptides_summary_png, dpi = 500)
        
        plt.show()
        plt.close()



# RUN THE SCRIPT ==============================================================
if __name__ == "__main__":

    # SET PARSER ==============================================================
    parser = argparse.ArgumentParser(
        description="""
        This script combines the .fastq of mRNA display cycles;
        it produces a summary report, analyses peptide relatedness,
        and DNA mutants.
        """
        )
    parser.add_argument(
        "--data_path",
        required=True,
        type=str,
        # action='store_true',
        help="""Provide path to directory containing .fastq files by round.
        Files have to end in two digits denoting selection cycle number."""
        )
    parser.add_argument(
        "--base_cycle",
        required=True,
        type=int,
        help="""Indicate selection base cycle number.
        This cycle will be used as a reference point for comparison."""
        )
    parser.add_argument(
        "--n_top_peptides",
        required=True,
        type=int,
        help="""Indicate the number of top peptides
        (most abundant peptides in base cycle) to include into graphs."""
        )
    parser.add_argument(
        "--file_name",
        required=True,
        type=str,
        # action='store_true',
        help="""Indicate selection name.
        This name is included in all generated files names."""
        )
    parser.add_argument(
        "--cdna_min_length",
        required=False,
        default=gp.CDNA_MIN_LENGTH,
        type=int,
        help="""Indicate minimum expected cDNA length."""
        )
    parser.add_argument(
        "--cdna_max_length",
        required=False,
        default=gp.CDNA_MAX_LENGTH,
        type=int,
        help="""Indicate maximum expected cDNA length."""
        )
    parser.add_argument(
        "--start_sequence",
        required=False,
        default=gp.START_SEQUENCE,
        type=str,
        help="""Indicate start sequence."""
        )
    parser.add_argument(
        "--stop_sequence",
        required=False,
        default=gp.STOP_SEQUENCE,
        type=str,
        help="""Indicate stop sequence."""
        )
    parser.add_argument(
        "--quality_score",
        required=False,
        default=gp.QUALITY_SCORE,
        type=int,
        help="""Indicate threshold quality score."""
        )

    parsed_args = vars(parser.parse_args())
    data_path = parsed_args['data_path']
    base_cycle = parsed_args['base_cycle']
    n_top_peptides = parsed_args['n_top_peptides']
    file_name = parsed_args['file_name']
    cdna_min_length = parsed_args['cdna_min_length']
    cdna_max_length = parsed_args['cdna_max_length']
    start_sequence = parsed_args['start_sequence']
    stop_sequence = parsed_args['stop_sequence']
    quality_score = parsed_args['quality_score']

    new_selection = Selection(
        data_path,
        base_cycle,
        n_top_peptides,
        file_name,
        cdna_min_length,
        cdna_max_length,
        start_sequence,
        stop_sequence,
        quality_score,
        )
    
    new_selection.get_summary()
    new_selection.get_dna_mutants_summary()
    new_selection.get_peptides_relatedness_summary()