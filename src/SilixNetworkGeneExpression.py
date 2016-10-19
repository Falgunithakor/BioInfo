import glob

import math
import os
import numpy as np
import re
import matplotlib.pyplot as plt
from scipy.stats import stats

from src.GeneExpressionManager import GeneExpressionManager

__author__ = 'FalguniT'


class SilixNetworkGeneExpression(object):
    def __init__(self, gene_expression_manager):
        self.network_data = []
        self.multiple_networks = False
        self.gene_expression_manager = gene_expression_manager
        self.silix_nw_file_path = ""
        self.silix_nw_exp_data_folder_path = ""
        self.silix_nw_exp_data_filename = "Network_Gene_Expression.txt"
        self.network_gene_expressions = []
        self.network_ehux_gene_expression = []
        self.network_geph_gene_expression = []

    def load_network_data(self):
        self.network_data = np.genfromtxt(self.silix_nw_file_path, delimiter='\n', dtype=str)
        if len(self.network_data.shape) == 0:
            self.network_data = np.genfromtxt(self.silix_nw_file_path, delimiter=',', dtype=str)
        print "Networks file data shape:", type(self.network_data), self.network_data.shape
        if ',' in self.network_data[0]:
            print "Multiple Networks in file"
            self.multiple_networks = True
        else:
            print "one network in file"

    def map_network_to_gene_expression(self):
        f_handle = open(self.silix_nw_exp_data_folder_path + self.silix_nw_exp_data_filename, 'w')
        for i, rowdata in enumerate(self.network_data):
            if self.multiple_networks:
                f_handle = open(
                    str(self.silix_nw_exp_data_folder_path + self.silix_nw_exp_data_filename).replace('Gene_Expression',
                                                                                                      'Gene_Expression_' + str(
                                                                                                          i)), 'w')
            for node in rowdata.split(","):
                gene_expression_row = self.gene_expression_manager.map_sequence_to_gene_expression(node)
                f_handle.write(str(gene_expression_row.tolist()).strip('[]').replace("'", "") + '\n')
        print "Generating Gene Expression data for networks completed"

    def generate_sequence_gene_expression_statistics(self, show_species_charts=True, show_chart=True):
        i = -1
        if self.multiple_networks:
            for nw_ge_file in glob.glob(self.silix_nw_exp_data_folder_path + '/*.txt'):
                i += 1
                mapping_data = np.genfromtxt(nw_ge_file, delimiter=',', dtype=str)
                if len(mapping_data) > 0:
                    # print 'Network: ', nw_ge_file, i, mapping_data.shape
                    x = np.array(mapping_data[:, 2], dtype=float)
                    y = np.array(mapping_data[:, 3], dtype=float)
                    ttest = stats.ttest_ind(x, y, equal_var=False)
                    nw_number = (int)(re.findall(r'\d+', nw_ge_file)[0])
                    nw_statistics = (
                        [nw_number, x.mean(), x.var(), x.std(), y.mean(), y.var(), y.std(), ttest.statistic,
                         ttest.pvalue])
                    self.network_gene_expressions.append(nw_statistics)
        else:
            mapping_data = np.genfromtxt(self.silix_nw_exp_data_folder_path + self.silix_nw_exp_data_filename,
                                         delimiter=',', dtype=str)
            if len(mapping_data) > 0:
                print 'Network: ', mapping_data.shape
                x = np.array(mapping_data[:, 2], dtype=float)
                y = np.array(mapping_data[:, 3], dtype=float)
                ttest = stats.ttest_ind(x, y, equal_var=False)
                nw_statistics = (
                    [0, x.mean(), x.var(), x.std(), y.mean(), y.var(), y.std(), ttest.statistic,
                     ttest.pvalue])
                self.network_gene_expressions.append(nw_statistics)
        # convert list into array
        self.network_gene_expressions = np.asarray(self.network_gene_expressions)

        # Save network gene expression statistics to csv file
        gene_expression_statistics_file = self.silix_nw_exp_data_folder_path + 'gene_expression_statistics.csv'
        if not os.path.isfile(gene_expression_statistics_file):
            with open(gene_expression_statistics_file, 'a') as f_handle:
                np.savetxt(f_handle, self.network_gene_expressions, delimiter=',')

        if show_species_charts:
            self.generate_species_wise_gene_expression_statistics()

        if show_chart:
            self.plot_chart()

    def generate_species_wise_gene_expression_statistics(self):
        ehux_data = []
        geph_data = []
        if self.multiple_networks:
            for nw_ge_file in glob.glob(self.silix_nw_exp_data_folder_path + '/*.txt'):
                mapping_data = np.genfromtxt(nw_ge_file, delimiter=',', dtype=str)
                if len(mapping_data) > 0:
                    nw_number = (int)(re.findall(r'\d+', nw_ge_file)[0])
                    if len([x for x in mapping_data if str(x[0]).isdigit()]) > 0:
                        ehux_data = np.vstack([x for x in mapping_data if str(x[0]).isdigit()])
                    if len([x for x in mapping_data if re.match(r'^evm.model.Contig', x[0])]) > 0:
                        geph_data = np.vstack([x for x in mapping_data if re.match(r'^evm.model.Contig', x[0])])
                    if len(ehux_data) > 0:
                        ehux_ca = np.array(ehux_data[:, 2], dtype=float)
                        ehux_spike = np.array(ehux_data[:, 3], dtype=float)
                        ehux_ttest = stats.ttest_ind(ehux_ca, ehux_spike, equal_var=False)
                        ehux_nw_statistics = (
                            [nw_number, ehux_ca.mean(), ehux_ca.var(), ehux_ca.std(), ehux_spike.mean(),
                             ehux_spike.var(),
                             ehux_spike.std(), ehux_ttest.statistic, ehux_ttest.pvalue])
                        self.network_ehux_gene_expression.append(ehux_nw_statistics)
                    if len(geph_data) > 0:
                        geph_ca = np.array(geph_data[:, 2], dtype=float)
                        geph_spike = np.array(geph_data[:, 3], dtype=float)
                        geph_ttest = stats.ttest_ind(geph_ca, geph_spike, equal_var=False)
                        geph_nw_statistics = (
                            [nw_number, geph_ca.mean(), geph_ca.var(), geph_ca.std(), geph_spike.mean(),
                             geph_spike.var(),
                             geph_spike.std(), geph_ttest.statistic, geph_ttest.pvalue])
                        self.network_geph_gene_expression.append(geph_nw_statistics)
        self.network_ehux_gene_expression = np.asarray(self.network_ehux_gene_expression)
        self.network_geph_gene_expression = np.asarray(self.network_geph_gene_expression)
        #print self.network_ehux_gene_expression

    def plot_chart(self):
        # PLOT Network Vs 9mM Ca factor Mean and Standard deviation  and Spike factor Mean and Standard Deviation
        chart1_9mmsd = self.network_gene_expressions[:, 3]
        chart1_spikesd = self.network_gene_expressions[:, 6]
        fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True, sharey=True)
        ax0.errorbar(self.network_gene_expressions[:, 0], self.network_gene_expressions[:, 1], yerr=chart1_9mmsd,
                     linestyle='None', fmt='-o')
        ax0.set_title('9mM Ca factor Mean and Standard deviation')
        ax0.set_ylabel(r'9mM Ca factor')
        ax0.axhline(0)
        ax1.errorbar(self.network_gene_expressions[:, 0], self.network_gene_expressions[:, 4], yerr=chart1_spikesd,
                     linestyle='None', fmt='o')
        ax1.set_title('Spike factor Mean and Standard deviation')
        ax1.set_ylabel(r'Spike factor')
        ax1.axhline(0)
        plt.xlim(-1, len(self.network_gene_expressions))
        plt.xlabel("Networks", fontsize=16)

        # PLOT 9mM Ca factor Mean vs Spike factor Mean with Standard Deviations
        chart2_data = self.network_gene_expressions
        chart2_9mmsd = chart2_data[:, 3]
        chart2_spikesd = chart2_data[:, 6]
        fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True, sharey=True)
        ax0.errorbar(chart2_data[:, 1], chart2_data[:, 4], yerr=chart2_9mmsd, linestyle='None', fmt='-o')
        ax0.set_title('9mM Ca factor Standard deviation')
        ax0.axhline(0)
        ax1.errorbar(chart2_data[:, 1], chart2_data[:, 4], yerr=chart2_spikesd, linestyle='None', fmt='o')
        ax1.set_title('Spike factor Standard deviation')
        ax1.axhline(0)
        plt.xlabel('9mM Ca factor', fontsize=16)
        plt.ylabel(r'Spike factor', fontsize=16)

        # PLOT Network vs Species mean and Standard Deviation
        if len(self.network_ehux_gene_expression) > 0 and len(self.network_geph_gene_expression) > 0:
            chart3_ehux9mmsd = self.network_ehux_gene_expression[:, 3]
            chart3_ehuxspikesd = self.network_ehux_gene_expression[:, 6]
            chart3_geph9mmsd = self.network_geph_gene_expression[:, 3]
            chart3_gephspikesd = self.network_geph_gene_expression[:, 6]
            fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True, sharey=True)
            ax0.errorbar(self.network_ehux_gene_expression[:, 0], self.network_ehux_gene_expression[:, 1],
                         yerr=chart3_ehux9mmsd,
                         linestyle='None', fmt='-o', label="Ehux")
            ax0.errorbar(self.network_geph_gene_expression[:, 0], self.network_geph_gene_expression[:, 1],
                         yerr=chart3_geph9mmsd,
                         linestyle='None', fmt='s', label="Geph")
            ax0.set_title('9mM Ca factor Mean and Standard deviation')
            ax0.set_ylabel(r'9mM Ca factor')
            ax0.legend()
            ax0.axhline(0)
            ax1.errorbar(self.network_ehux_gene_expression[:, 0], self.network_ehux_gene_expression[:, 4],
                         yerr=chart3_ehuxspikesd,
                         linestyle='None', fmt='o', label="Ehux")
            ax1.errorbar(self.network_geph_gene_expression[:, 0], self.network_geph_gene_expression[:, 4],
                         yerr=chart3_gephspikesd,
                         linestyle='None', fmt='s', label="Geph")
            ax1.set_title('Spike factor Mean and Standard deviation')
            ax1.set_ylabel(r'Spike factor')
            ax1.legend()
            ax1.axhline(0)
            plt.xlim(-1, len(self.network_ehux_gene_expression))
            plt.xlabel("Networks", fontsize=16)
        # show all plots
        plt.show()


if __name__ == "__main__":
    # Instantiate Gene Expression Manager object
    gene_expression_manager = GeneExpressionManager()
    # LOAD Gene Expression data
    gene_expression_manager.load_gene_expression_data()

    # Create Object of Class
    objSilixNw = SilixNetworkGeneExpression(gene_expression_manager)
    # Provide Network File path and output folder path
    objSilixNw.silix_nw_file_path = "../output/silix/more_than_six_nodes/node_data_geph_ehux.txt"
    objSilixNw.silix_nw_exp_data_folder_path = "../output/silix/more_than_six_nodes/gene_expr/"

    # LOAD Network Cluster data
    objSilixNw.load_network_data()

    # MAP Network sequences to gene expression sequences
    # objSilixNw.map_network_to_gene_expression()

    # GENERATE scatter plot and statistics
    objSilixNw.generate_sequence_gene_expression_statistics(show_chart=True)


