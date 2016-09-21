import numpy as np
import re

__author__ = 'FalguniT'


class GeneExpressionManager(object):
    def __init__(self):
        self.geph_gene_expression_file_path = "../data/gene_expressions/Geph_voom_coef.txt"
        self.ehux_gene_expression_file_path = "../data/gene_expressions/ehux_voom_coef.txt"
        self.iso_gene_expression_file_path = "../data/gene_expressions/iso_voom_coef_fix.txt"
        self.chry_gene_expression_file_path = ""
        self.geph_gene_expression_data = []
        self.ehux_gene_expression_data = []
        self.iso_gene_expression_data = []

    def load_gene_expression_data(self):
        self.geph_gene_expression_data = np.genfromtxt(self.geph_gene_expression_file_path, delimiter='\t', dtype=str,
                                                       skip_header=1)
        self.ehux_gene_expression_data = np.genfromtxt(self.ehux_gene_expression_file_path, delimiter='\t', dtype=str,
                                                       skip_header=1)
        self.iso_gene_expression_data = np.genfromtxt(self.iso_gene_expression_file_path, delimiter='\t',
                                                      dtype=str, skip_header=1)

    def map_sequence_to_gene_expression(self, sequence_no):
        gene_expression_data = []
        if re.match(r'^\d', sequence_no):
            gene_expression_data = self.ehux_gene_expression_data
        elif re.match(r'^evm.model.Contig', sequence_no):
            gene_expression_data = self.geph_gene_expression_data
        elif re.match(r'^evm.model.scaffold', sequence_no):
            gene_expression_data = self.iso_gene_expression_data
        # read gene expression line assuming first column as sequence name
        gene_expression_row = gene_expression_data[np.where(gene_expression_data[:, 0] == sequence_no)]
        return gene_expression_row
