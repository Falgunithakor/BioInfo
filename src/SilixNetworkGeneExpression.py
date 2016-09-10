import glob
import os
from Bio import SearchIO

import numpy as np
import re
import matplotlib.pyplot as plt

__author__ = 'FalguniT'


class SilixNetworkGeneExpression(object):
    def __init__(self):
        self.network_data = []
        self.geph_gene_expression_data = []
        self.ehux_gene_expression_data = []
        self.iso_gene_expression_file_data = []

        self.geph_gene_expression_file_path = "../data/Geph_voom_coef.txt"
        self.ehux_gene_expression_file_path = "../data/ehux_voom_coef.txt"
        self.iso_gene_expression_file_path = "../data/iso_voom_coef_fix.txt"

        self.silix_nw_file_path = "../output/silix/more_than_six_nodes/node_data_geph_iso_ehux.txt"
        self.silix_nw_exp_data_file_path = "../output/silix/more_than_six_nodes_GephEhuxISO_Gene_Expression.txt"

    def load_gene_expression_data(self):
        self.geph_gene_expression_data = np.genfromtxt(self.geph_gene_expression_file_path, delimiter='\t', dtype=str)
        # print self.gene_expression_data
        print type(self.geph_gene_expression_data), self.geph_gene_expression_data.shape
        self.ehux_gene_expression_data = np.genfromtxt(self.ehux_gene_expression_file_path, delimiter='\t', dtype=str)
        print type(self.ehux_gene_expression_data), self.ehux_gene_expression_data.shape
        self.iso_gene_expression_file_data = np.genfromtxt(self.iso_gene_expression_file_path, delimiter='\t',
                                                           dtype=str)
        print type(self.iso_gene_expression_file_data), self.iso_gene_expression_file_data.shape

    def load_network_data(self):
        self.network_data = np.genfromtxt(self.silix_nw_file_path, delimiter='\n', dtype=str)
        # print self.network_data
        print type(self.network_data), self.network_data.shape

    def map_network_to_gene_expression(self):
        with open(self.silix_nw_exp_data_file_path, 'w') as f_handle:
            for rowdata in self.network_data[:]:
                for node in rowdata.split(","):
                    print node
                    if node in self.ehux_gene_expression_data[:]:
                        for r in (self.ehux_gene_expression_data[
                                      np.where(self.ehux_gene_expression_data[:, 0] == node)])[
                                 :]:
                            f_handle.write(str(r[0] + ',' + r[1] + ',' + r[2] + ',' + r[3]) + '\n')
                    elif node in self.geph_gene_expression_data[:]:
                        for r in (self.geph_gene_expression_data[
                                      np.where(self.geph_gene_expression_data[:, 0] == node)])[
                                 :]:
                            f_handle.write(str(r[0] + ',' + r[1] + ',' + r[2] + ',' + r[3]) + '\n')
                    elif node in self.iso_gene_expression_file_data[:]:
                        for r in (self.iso_gene_expression_file_data[
                                      np.where(self.iso_gene_expression_file_data[:, 0] == node)])[
                                 :]:
                            f_handle.write(str(r[0] + ',' + r[1] + ',' + r[2] + ',' + r[3]) + '\n')
                    else:
                        f_handle.write(node + ',0,0,0' + '\n')
                        print node

    def load_nw_gene_expression_file(self):
        data = np.genfromtxt(self.silix_nw_exp_data_file_path, delimiter=',', dtype=str)
        print type(data), data.shape
        colors = ['b', 'c']
        x = data[:, 1]
        y = data[:, 2]
        plt.scatter(x, y)
        plt.title("More than six node network clusters for Geph Ehux ISO", fontsize=20)
        plt.xlabel('Intercept', fontsize=16)
        plt.ylabel('Ca mm factor', fontsize=16)
        plt.show()

if __name__ == "__main__":
    objSilixNw = SilixNetworkGeneExpression()
    '''
    objSilixNw.load_network_data()
    objSilixNw.load_gene_expression_data()
    objSilixNw.map_network_to_gene_expression()
    '''
    objSilixNw.load_nw_gene_expression_file()

