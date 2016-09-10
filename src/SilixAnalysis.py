from collections import Counter
import glob
import os

import numpy as np
from itertools import groupby
import pandas as pd
import re

from sklearn.feature_selection.univariate_selection import f_classif

__author__ = 'FalguniT'


class SilixAnalysis(object):
    def __init__(self):
        self.silix_output_file_path = "../output/silix/output.txt"
        self.input_file_path = "../data/Ehux_DE_IDs.txt"
        self.input_mapping_file_path = "../output/silix/more_than_one_nodes/*.txt"
        self.silix_input_output_mapping_data_file_path = "../output/silix/Ehux_DE_IDs_nw_category_mapping.txt"
        self.input_data = []
        self.silix_output_data = []

    def write_output_file(self, node, data, category):
        category = str(category).replace("../output/silix/more_than_one_nodes\\node_data_", "").replace(".txt", "")
        print node, data, category
        with open(self.silix_input_output_mapping_data_file_path, 'a') as f_handle:
            f_handle.write(str(node) + ', ' + str(data).strip('[]') + ', ' + category + '\n')

    def load_data(self):
        self.input_data = np.genfromtxt(self.input_file_path, delimiter='\t', dtype=str)
        print self.input_data
        print type(self.input_data), self.input_data.shape
        self.silix_output_data = np.genfromtxt(self.silix_output_file_path, delimiter='\t', dtype=None)


    def map_input_to_target_data(self):
        for index, row_data in enumerate(self.input_data[:, 1]):
            print  row_data
            row_data_found = False
            for nw_category_file in glob.glob(self.input_mapping_file_path):
                mapping_data = np.genfromtxt(nw_category_file, delimiter='\n', dtype=str)
                if len([i for i, item in enumerate(mapping_data) if row_data in item]) > 0:
                    self.write_output_file(self.input_data[index, 0] + ", " + row_data,
                                           mapping_data[[i for i, item in enumerate(mapping_data) if row_data in item]],
                                           nw_category_file)
                    row_data_found = True
                    break
            if row_data_found is False:
                # check if exists in silix output file
                if len([i for i, item in enumerate(self.silix_output_data) if row_data in item]):
                    self.write_output_file(self.input_data[index, 0] + ", " +row_data, row_data, "Singled_out")
                else:
                    self.write_output_file(self.input_data[index, 0] + ", " +row_data, row_data, "Not_found_in_Silix_output")


if __name__ == "__main__":
    objSilix = SilixAnalysis()
    objSilix.load_data()
    objSilix.map_input_to_target_data()

