import glob

import math
import os
import numpy as np
import re
import matplotlib.pyplot as plt
from matplotlib import gridspec
import scipy.stats as stats

from src.GeneExpressionManager import GeneExpressionManager

__author__ = 'FalguniT'


class SilixNetworkGeneExpression(object):
    def __init__(self, gene_expression_manager):
        self.network_data = []
        self.multiple_networks = False
        self.gene_expression_manager = gene_expression_manager
        self.silix_nw_file_path = ""
        self.output_silix_nw_exp_data_folder_path = ""
        self.silix_nw_exp_data_filename = "Network_Gene_Expression.txt"
        self.network_gene_expressions = []
        self.network_ehux_gene_expression = []
        self.network_geph_gene_expression = []
        self.network_iso_gene_expression = []

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
        if not os.path.exists(self.output_silix_nw_exp_data_folder_path):
            os.makedirs(self.output_silix_nw_exp_data_folder_path)

    def map_network_to_gene_expression(self):
        f_handle = open(self.output_silix_nw_exp_data_folder_path + self.silix_nw_exp_data_filename, 'w')
        for i, rowdata in enumerate(self.network_data):
            if self.multiple_networks:
                f_handle = open(
                    str(self.output_silix_nw_exp_data_folder_path + self.silix_nw_exp_data_filename).replace(
                        'Gene_Expression',
                        'Gene_Expression_' + str(
                            i)), 'w')
            for node in rowdata.split(","):
                gene_expression_row = self.gene_expression_manager.map_sequence_to_gene_expression(node)
                f_handle.write(str(gene_expression_row.tolist()).strip('[]').replace("'", "") + '\n')
        print "Generating Gene Expression data for networks completed"

    def generate_sequence_gene_expression_statistics(self, show_species_charts=True, show_chart=True):
        i = -1
        if self.multiple_networks:
            for nw_ge_file in glob.glob(self.output_silix_nw_exp_data_folder_path + '/*.txt'):
                i += 1
                mapping_data = np.genfromtxt(nw_ge_file, delimiter=',', dtype=str)
                if len(mapping_data) > 0:
                    print 'Network: ', i, mapping_data.shape
                    x = np.array(mapping_data[:, 2], dtype=float)
                    y = np.array(mapping_data[:, 3], dtype=float)
                    ca_stat = ca_pvalue = spike_stat = spike_pvalue = ind_stat = ind_pvalue = 0

                    if not np.all(x == 0):
                        ca_stat, ca_pvalue = stats.ttest_1samp(x[x != 0], 0)
                        spike_stat, spike_pvalue = stats.ttest_1samp(y[y != 0], 0)
                        ind_stat, ind_pvalue = stats.ttest_ind(x[x != 0], y[y != 0], equal_var=False)
                    nw_number = (int)(re.findall(r'\d+', nw_ge_file)[0])
                    nw_statistics = (
                        [nw_number, x[x != 0].mean(), x[x != 0].var(), x[x != 0].std(), y[y != 0].mean(),
                         y[y != 0].var(), y[y != 0].std(), ca_stat, ca_pvalue,
                         spike_stat, spike_pvalue, ind_stat, ind_pvalue])
                    self.network_gene_expressions.append(nw_statistics)
        else:
            mapping_data = np.genfromtxt(self.output_silix_nw_exp_data_folder_path + self.silix_nw_exp_data_filename,
                                         delimiter=',', dtype=str)
            if len(mapping_data) > 0:
                print 'Network: ', mapping_data.shape
                x = np.array(mapping_data[:, 2], dtype=float)
                y = np.array(mapping_data[:, 3], dtype=float)
                ca_stat = ca_pvalue = spike_stat = spike_pvalue = ind_stat = ind_pvalue = 0
                if not np.all(x == 0):
                    ca_stat, ca_pvalue = stats.ttest_1samp(x[x != 0], 0)
                    spike_stat, spike_pvalue = stats.ttest_1samp(y[y != 0], 0)
                    ind_stat, ind_pvalue = stats.ttest_ind(x[x != 0], y[y != 0], equal_var=False)
                nw_statistics = (
                    [0, x[x != 0].mean(), x[x != 0].var(), x[x != 0].std(), y[y != 0].mean(),
                     y[y != 0].var(), y[y != 0].std(), ca_stat, ca_pvalue,
                     spike_stat, spike_pvalue, ind_stat, ind_pvalue])
                self.network_gene_expressions.append(nw_statistics)
        # convert list into array
        self.network_gene_expressions = np.asarray(self.network_gene_expressions)

        # Save network gene expression statistics to csv file
        gene_expression_statistics_file = self.output_silix_nw_exp_data_folder_path + 'gene_expression_statistics.csv'
        with open(gene_expression_statistics_file, 'w') as f_handle:
            f_handle.write(
                'Network, 9mM CA Mean, 9mM CA Var, 9mM CA SD, Spike Mean, Spike Var, Spike SD, 9mM CA ttest-stat, 9mM CA ttest-pvalue, Spike ttest-stat, Spike ttest-pvalue, Ind ttest-stat, Ind ttest-pvalue  \n')
            np.savetxt(f_handle, self.network_gene_expressions, delimiter=',')

        if show_species_charts:
            self.generate_species_wise_gene_expression_statistics()

        if self.multiple_networks and show_chart:
            self.plot_all_nw_gene_expr_stats_chart()
        elif show_chart and not self.multiple_networks:
            self.plot_single_network_gene_expr_stats_chart()

    def generate_species_wise_gene_expression_statistics(self):
        ehux_data = []
        geph_data = []
        iso_data = []
        for nw_ge_file in glob.glob(self.output_silix_nw_exp_data_folder_path + '/*.txt'):
            mapping_data = np.genfromtxt(nw_ge_file, delimiter=',', dtype=str)
            if len(mapping_data) > 0:
                if self.multiple_networks:
                    nw_number = (int)(re.findall(r'\d+', nw_ge_file)[0])
                else:
                    nw_number = 0
                if len([x for x in mapping_data if str(x[0]).isdigit()]) > 0:
                    ehux_data = np.vstack([x for x in mapping_data if str(x[0]).isdigit()])
                if len([x for x in mapping_data if re.match(r'^evm.model.Contig', x[0])]) > 0:
                    geph_data = np.vstack([x for x in mapping_data if re.match(r'^evm.model.Contig', x[0])])
                if len([x for x in mapping_data if re.match(r'^evm.model.scaffold', x[0])]) > 0:
                    iso_data = np.vstack([x for x in mapping_data if re.match(r'^evm.model.scaffold', x[0])])
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
                if len(iso_data) > 0:
                    iso_ca = np.array(iso_data[:, 2], dtype=float)
                    iso_spike = np.array(iso_data[:, 3], dtype=float)
                    iso_ttest = stats.ttest_ind(iso_ca, iso_spike, equal_var=False)
                    iso_nw_statistics = (
                        [nw_number, iso_ca.mean(), iso_ca.var(), iso_ca.std(), iso_spike.mean(),
                         iso_spike.var(),
                         iso_spike.std(), iso_ttest.statistic, iso_ttest.pvalue])
                    self.network_iso_gene_expression.append(iso_nw_statistics)
        self.network_ehux_gene_expression = np.asarray(self.network_ehux_gene_expression)
        self.network_geph_gene_expression = np.asarray(self.network_geph_gene_expression)
        self.network_iso_gene_expression = np.asarray(self.network_iso_gene_expression)
        # print self.network_ehux_gene_expression

    def plot_all_nw_gene_expr_stats_chart(self):
        # PLOT Network Vs 9mM Ca factor Mean and Standard deviation  and Spike factor Mean and Standard Deviation
        chart1_9mmsd = self.network_gene_expressions[:, 3]
        chart1_spikesd = self.network_gene_expressions[:, 6]
        fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True, sharey=True)
        ax0.errorbar(self.network_gene_expressions[:, 0], self.network_gene_expressions[:, 1], yerr=chart1_9mmsd,
                     linestyle='None', fmt='-o')
        self.apply_network_plot_style(ax0, title='9mM Ca factor Mean and Standard deviation', ylabel='9mM Ca factor (mean)')
        ax1.errorbar(self.network_gene_expressions[:, 0], self.network_gene_expressions[:, 4], yerr=chart1_spikesd,
                     linestyle='None', fmt='o')
        self.apply_network_plot_style(ax1, title='Spike factor Mean and Standard deviation', ylabel='Spike factor (mean)',
                                      xlabel='Networks')

        # PLOT 9mM Ca factor Mean vs Spike factor Mean with Standard Deviations
        chart2_data = self.network_gene_expressions
        chart2_9mmsd = chart2_data[:, 3]
        chart2_spikesd = chart2_data[:, 6]
        fig, (ax0, ax1) = plt.subplots(nrows=2)
        ax0.errorbar(chart2_data[:, 1], chart2_data[:, 4], yerr=chart2_9mmsd, linestyle='None', fmt='-o')
        ax0.set_title('9mM Ca factor Mean and Standard deviation')
        ax0.axhline(0)
        ax0.axvline(0)
        ax0.set_xlim(-4, 6)
        ax0.set_ylim(-4, 6)
        ax0.set_aspect('equal', adjustable='box')
        ax0.grid(True)
        ax0.set_ylabel(r'Spike factor', fontsize=16)
        ax1.errorbar(chart2_data[:, 1], chart2_data[:, 4], yerr=chart2_spikesd, linestyle='None', fmt='o', ecolor='r')
        ax1.set_title('Spike factor Mean and Standard deviation')
        ax1.axhline(0)
        ax1.axvline(0)
        ax1.set_xlim(-4, 6)
        ax1.set_ylim(-4, 6)
        ax1.set_aspect('equal', adjustable='box')
        plt.grid(True)
        plt.xlabel('9mM Ca factor (mean)', fontsize=16)
        plt.ylabel(r'Spike factor (mean)', fontsize=16)

        # PLOT Network vs Species mean and Standard Deviation
        if len(self.network_ehux_gene_expression) > 0 and len(self.network_geph_gene_expression) > 0:
            chart3_ehux9mmsd = self.network_ehux_gene_expression[:, 3]
            chart3_ehuxspikesd = self.network_ehux_gene_expression[:, 6]
            chart3_geph9mmsd = self.network_geph_gene_expression[:, 3]
            chart3_gephspikesd = self.network_geph_gene_expression[:, 6]
            if len(self.network_iso_gene_expression) > 0:
                fig, ((ax0, ax2), (ax1, ax3), (ax4, ax5)) = plt.subplots(nrows=3, ncols=2, sharex=True, sharey=True)
            else:
                fig, ((ax0, ax2), (ax1, ax3)) = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True)
            ax0.errorbar(self.network_ehux_gene_expression[:, 0], self.network_ehux_gene_expression[:, 1],
                         yerr=chart3_ehux9mmsd,
                         linestyle='None', fmt='-o', label="Ehux")
            self.apply_network_plot_style(ax0, title='9mM Ca factor Mean and Standard deviation',
                                          ylabel='9mM Ca factor (mean)',
                                          xlabel='Networks')
            ax1.errorbar(self.network_geph_gene_expression[:, 0], self.network_geph_gene_expression[:, 1],
                         yerr=chart3_geph9mmsd,
                         linestyle='None', fmt='s', label="Geph")
            self.apply_network_plot_style(ax1, ylabel='9mM Ca factor (mean)', xlabel='Networks')
            ax2.errorbar(self.network_ehux_gene_expression[:, 0], self.network_ehux_gene_expression[:, 4],
                         yerr=chart3_ehuxspikesd,
                         linestyle='None', fmt='o', label="Ehux", ecolor='red')
            self.apply_network_plot_style(ax2, title='Spike factor Mean and Standard deviation', ylabel='Spike factor (mean)',
                                          xlabel='Networks')
            ax3.errorbar(self.network_geph_gene_expression[:, 0], self.network_geph_gene_expression[:, 4],
                         yerr=chart3_gephspikesd,
                         linestyle='None', fmt='s', label="Geph", ecolor='red')
            self.apply_network_plot_style(ax3, ylabel='Spike factor (mean)', xlabel='Networks')
            if len(self.network_iso_gene_expression) > 0:
                ax4.errorbar(self.network_iso_gene_expression[:, 0], self.network_iso_gene_expression[:, 4],
                         yerr=self.network_iso_gene_expression[:, 3],
                         linestyle='None', fmt='o', label="Iso")
                self.apply_network_plot_style(ax4, ylabel='9mM ca factor (mean)',
                                              xlabel='Networks')
                ax5.errorbar(self.network_iso_gene_expression[:, 0], self.network_iso_gene_expression[:, 4],
                             yerr=self.network_iso_gene_expression[:, 6],
                             linestyle='None', fmt='s', label="Iso", ecolor='red')
                self.apply_network_plot_style(ax5, ylabel='Spike factor (mean)', xlabel='Networks')
            #fig.subplots_adjust(hspace=0)
        # show all plots
        plt.show()

    def apply_network_plot_style(self, ax, title="", xlabel="", ylabel=""):
        ax.set_title(title)
        ax.set_xlabel(xlabel, fontsize=16)
        ax.set_ylabel(ylabel, fontsize=16)
        ax.legend()
        ax.axhline(0)
        ax.set_xlim(-1, len(self.network_data))
        minor_ticks = np.arange(0, len(self.network_data), 2)
        ax.set_xticks(minor_ticks, minor=True)
        ax.grid(which='both')

    def generate_network_charts_from_ttest(self, ttest_for='ca', minvalue=0.01, show_chart=False):
        nw_gene_expr_stats = []
        # CHECK Network Gene Expression Statistics does not exist
        if os.path.exists(self.output_silix_nw_exp_data_folder_path + 'gene_expression_statistics.csv'):
            nw_gene_expr_stats = np.genfromtxt(
                self.output_silix_nw_exp_data_folder_path + 'gene_expression_statistics.csv', delimiter=',',
                skip_header=1)
        if len(nw_gene_expr_stats) > 0 and not self.multiple_networks:
            nw_gene_expr_stats = nw_gene_expr_stats.reshape((1, 13))

        if len(nw_gene_expr_stats) > 0:
            # Remove network rows with all zeros statistics
            nw_gene_expr_stats = nw_gene_expr_stats[~np.all(nw_gene_expr_stats[:, 1:13] == 0, axis=1)]
            print nw_gene_expr_stats.shape
            # Read ttest pvalue for specific ttest
            print ttest_for
            nw_gene_expr_stats = nw_gene_expr_stats[np.where(nw_gene_expr_stats[:, 8] != 0)]
            if ttest_for == 'ca':
                nw_gene_expr_stats = nw_gene_expr_stats[np.where(nw_gene_expr_stats[:, 8] <= minvalue)]
                # nw_gene_expr_stats = nw_gene_expr_stats[np.where(nw_gene_expr_stats[:, 8] != 'nan')]
                print nw_gene_expr_stats.shape, nw_gene_expr_stats

            elif ttest_for == 'spike':
                nw_gene_expr_stats = nw_gene_expr_stats[np.where(nw_gene_expr_stats[:, 10] <= minvalue)]
            elif ttest_for == 'ca-spike':
                nw_gene_expr_stats = nw_gene_expr_stats[np.where(nw_gene_expr_stats[:, 12] <= minvalue)]
            print nw_gene_expr_stats.shape
            # print nw_gene_expr_stats

        # PlOT network scatter charts
        chartitle = 'Ehux-Geph NWs for ttest <= %s for %s factor' % (minvalue, ttest_for)
        if show_chart and self.multiple_networks:
            self.plot_nw_specific_gene_expr_scatter_charts(nw_gene_expr_stats, chartitle)
        elif show_chart and not self.multiple_networks:
            fig = plt.figure()
            chart_data = np.empty((0, 4))
            for node in self.network_data:
                gene_expr_row = self.gene_expression_manager.map_sequence_to_gene_expression(node)
                chart_data = np.append(chart_data, gene_expr_row, axis=0)
            print chart_data
            x = chart_data[:, 2]
            y = chart_data[:, 3]
            plt.scatter(x, y)
            plt.axhline(0)
            plt.axvline(0)
            plt.ylim([-5, 5])
            plt.xlim([-5, 5])
            fig.suptitle(chartitle, fontsize=20)
            fig.text(0.5, 0.04, '9mM Ca factor (mean)', ha='center')
            fig.text(0.04, 0.5, 'Spike factor (mean)', va='center', rotation='vertical')
            plt.show()

    def plot_nw_specific_gene_expr_scatter_charts(self, nw_gene_expr_stats, chartitle=''):
        # print len(nw_gene_expr_stats)
        number_of_subplots = len(nw_gene_expr_stats)
        cols = 2
        rows = int(math.ceil(number_of_subplots / float(cols)))
        gs = gridspec.GridSpec(rows, cols)
        fig = plt.figure()
        fig.subplots_adjust(wspace=0, hspace=0)
        for i, v in enumerate(xrange(number_of_subplots)):
            ax1 = plt.subplot(gs[i], aspect='equal')
            ttest_value = 'ttest(ca): %.4f\nttest(spike): %.4f' % (nw_gene_expr_stats[i, 8], nw_gene_expr_stats[i, 10])
            self.plot_network_gene_expr_scatter_chart(int(nw_gene_expr_stats[i, 0]), ttest_value, ax1)
        fig.suptitle(chartitle, fontsize=24)
        fig.text(0.5, 0.04, '9mM Ca factor (mean)', ha='center', fontsize=16)
        fig.text(0.04, 0.5, 'Spike factor (mean)', va='center', rotation='vertical', fontsize=16)
        plt.tight_layout()
        plt.show()

    def plot_network_gene_expr_scatter_chart(self, nw_number, ttest_for, ax1):
        # Read GeneExpression data from network file
        print nw_number
        chart_data = np.empty((0, 4))
        for node in self.network_data[nw_number].split(","):
            gene_expr_row = self.gene_expression_manager.map_sequence_to_gene_expression(node)
            chart_data = np.append(chart_data, gene_expr_row, axis=0)
        x = chart_data[:, 2]
        y = chart_data[:, 3]
        ax1.scatter(x, y)
        ax1.axhline(0)
        ax1.axvline(0)
        ax1.set_ylim([-5, 5])
        ax1.set_xlim([-5, 5])
        ax1.text(1, 3, 'NW# %s (%s nodes)\n%s' % (nw_number, len(chart_data), ttest_for), style='italic',
                 bbox={'facecolor': 'green', 'alpha': 0.1, 'pad': 1})

    def plot_single_network_gene_expr_stats_chart(self):
        # PLOT Network Vs 9mM Ca factor Mean and Standard deviation  and Spike factor Mean and Standard Deviation
        chart1_9mmsd = self.network_gene_expressions[:, 3]
        chart1_spikesd = self.network_gene_expressions[:, 6]
        fig, (ax0) = plt.subplots(nrows=1, sharex=True, sharey=True)
        ax0.errorbar(self.network_gene_expressions[:, 1], self.network_gene_expressions[:, 4], xerr=chart1_9mmsd,
                     yerr=chart1_spikesd,
                     linestyle='None', fmt='-o', label='Standard Deviation')
        ax0.axhline(0)
        ax0.set_xlim(-4, 6)
        ax0.set_ylim(-4, 6)
        ax0.set_aspect('equal', adjustable='box')
        ax0.legend()
        plt.grid(True)
        plt.xlabel('9mM Ca factor (mean)', fontsize=16)
        plt.ylabel('Spike factor (mean)', fontsize=16)

        print self.network_ehux_gene_expression, self.network_geph_gene_expression
        # PLOT Network vs Species mean and Standard Deviation
        if len(self.network_ehux_gene_expression) > 0 and len(self.network_geph_gene_expression) > 0:
            chart3_ehux9mmsd = self.network_ehux_gene_expression[:, 3]
            chart3_ehuxspikesd = self.network_ehux_gene_expression[:, 6]
            chart3_geph9mmsd = self.network_geph_gene_expression[:, 3]
            chart3_gephspikesd = self.network_geph_gene_expression[:, 6]
            if len(self.network_geph_gene_expression) > 0:
                fig, ((ax0, ax1, ax3)) = plt.subplots(nrows=3, sharex=True, sharey=True)
            else:
                fig, ((ax0, ax1)) = plt.subplots(nrows=2, sharex=True, sharey=True)
            ax0.errorbar(self.network_ehux_gene_expression[:, 1], self.network_ehux_gene_expression[:, 4],
                         xerr=chart3_ehux9mmsd, yerr=chart3_ehuxspikesd,
                         linestyle='None', fmt='-o', label="Ehux")
            ax0.legend()
            ax1.errorbar(self.network_geph_gene_expression[:, 1], self.network_geph_gene_expression[:, 4],
                         xerr=chart3_geph9mmsd, yerr=chart3_gephspikesd,
                         linestyle='None', fmt='s', label="Geph")
            ax1.legend()
            if len(self.network_geph_gene_expression) > 0:
                ax3.errorbar(self.network_iso_gene_expression[:, 1], self.network_iso_gene_expression[:, 4],
                         xerr=self.network_iso_gene_expression[:, 3], yerr=self.network_iso_gene_expression[:, 6],
                         linestyle='None', fmt='o', label="Iso")
                ax3.legend()
            fig.subplots_adjust(hspace=0)
        # show all plots
        plt.show()


if __name__ == "__main__":
    # Instantiate Gene Expression Manager object
    gene_expression_manager = GeneExpressionManager()
    # LOAD Gene Expression data
    gene_expression_manager.load_gene_expression_data()

    # Create Object of Class
    objNwGE = SilixNetworkGeneExpression(gene_expression_manager)
    # Provide Network File path and output folder path
    objNwGE.silix_nw_file_path = "../data/additiona_analysis/hapto_networks_components.txt"
    objNwGE.output_silix_nw_exp_data_folder_path = "../output/hapto_networks/gene_expr/"
    #objNwGE.silix_nw_file_path = "../output/silix/more_than_six_nodes/node_data_geph_ehux.txt"
    #objNwGE.output_silix_nw_exp_data_folder_path = "../output/silix_more_than_six_nodes/gene_expr/"

    # LOAD Network Cluster data
    objNwGE.load_network_data()

    # MAP Network sequences to gene expression sequences
    objNwGE.map_network_to_gene_expression()

    # GENERATE scatter plot and statistics
    objNwGE.generate_sequence_gene_expression_statistics(show_chart=True)

    # SELECT specific networks based on pvalue statistics and generate charts
    objNwGE.generate_network_charts_from_ttest(ttest_for='ca', minvalue=1, show_chart=True)
    #objNwGE.generate_network_charts_from_ttest(ttest_for='spike', minvalue=1, show_chart=True)
    #objNwGE.generate_network_charts_from_ttest(ttest_for='ca-spike', minvalue=0.001, show_chart=True)

