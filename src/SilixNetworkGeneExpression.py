import glob
import os
from Bio import SearchIO

import numpy as np
from itertools import groupby
import pandas as pd
import re

__author__ = 'FalguniT'


class SilixNetworkGeneExpression(object):
    def __init__(self):
        self.network_data = []
        self.geph_gene_expression_data = []
        self.ehux_gene_expression_data = []
        self.geph_gene_expression_file_path = "../data/Geph_voom_coef.txt"
        self.ehux_gene_expression_file_path = "../data/ehux_voom_coef.txt"
        self.silix_nw_file_path = "../output/silix/node_data_232.txt"
        self.silix_nw_exp_data_file_path = "../output/silix/node_232_geph_gene_exp.txt"
        self.blast_data_path = "../data/testblastdata/*.o"

    def load_gene_expression_data(self):
        self.geph_gene_expression_data = np.genfromtxt(self.geph_gene_expression_file_path, delimiter='\t', dtype=str)
        # print self.gene_expression_data
        print type(self.geph_gene_expression_data), self.geph_gene_expression_data.shape
        self.ehux_gene_expression_data = np.genfromtxt(self.ehux_gene_expression_file_path, delimiter='\t', dtype=str)
        print type(self.ehux_gene_expression_data), self.ehux_gene_expression_data.shape

    def load_network_data(self):
        self.network_data = np.genfromtxt(self.silix_nw_file_path, delimiter=',', dtype=str)
        # print self.network_data
        print type(self.network_data), self.network_data.shape

    def map_network_to_gene_expression(self):
        with open(self.silix_nw_exp_data_file_path, 'w') as f_handle:
            for node in self.network_data:
                if node in self.ehux_gene_expression_data[:]:
                    for r in (self.ehux_gene_expression_data[np.where(self.ehux_gene_expression_data[:, 0] == node)])[
                             :]:
                        f_handle.write(str(r[0] + ',' + r[1] + ',' + r[2] + ',' + r[3]) + '\n')
                elif node in self.geph_gene_expression_data[:]:
                    for r in (self.geph_gene_expression_data[np.where(self.geph_gene_expression_data[:, 0] == node)])[
                             :]:
                        f_handle.write(str(r[0] + ',' + r[1] + ',' + r[2] + ',' + r[3]) + '\n')
                else:
                    f_handle.write(node + ',0,0,0' + '\n')


    def load_nw_gene_expression_file(self):
        data = np.genfromtxt(self.silix_nw_exp_data_file_path, delimiter=' ', dtype=str)
        print type(data), data.shape


if __name__ == "__main__":
    objSilixNw = SilixNetworkGeneExpression()

    objSilixNw.load_network_data()
    objSilixNw.load_gene_expression_data()
    objSilixNw.map_network_to_gene_expression()
    '''
    objSilixNw.load_nw_gene_expression_file()
    '''
