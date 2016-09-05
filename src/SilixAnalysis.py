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
        self.data = []

    @staticmethod
    def write_network_file(file_path, networkdata):
        with open(file_path, 'w') as f_handle:
            for key, value in sorted(networkdata.items()):
                f_handle.write(str(value) + '\n')

    def write_all_network_file(self, node_count, value):
        file_path = '../output/silix/more_than_six_nodes/node_data_%s.txt' % node_count
        #print node_count, value
        with open(file_path, 'a') as f_handle:
            f_handle.write(str(value).replace("[", "").replace("]", "") + '\n')

    def load_actual_data(self):
        self.data = np.genfromtxt(self.silix_output_file_path, delimiter='\t', dtype=None)
        # print self.data

    def get_results_statistics(self):
        total = 0
        count_array = np.zeros(233, dtype=int)
        species_specific = ["all_geph", "all_iso", "all_chrys", "all_ehux",
                            "geph_iso", "geph_chry", "geph_ehux", "iso_chry", "iso_ehux", "chry_ehux",
                            "geph_iso_ehux", "geph_iso_chry", "iso_ehux_chry", "geph_ehux_chry", "species"]
        for file in species_specific:
            file_path = "../output/silix/more_than_six_nodes/node_data_%s.txt" % file
            if os.path.exists(file_path):
                data = np.genfromtxt(file_path, delimiter='\n', dtype=str)
                if len(data.shape) > 0:
                    print "count of", file, ":", data.shape[0], data.shape
                    total = total + data.shape[0]
                    for i in range(0, data.shape[0]):
                        nodes = data[i].count(',') + 1
                        count_array[nodes] += 1
                else:
                    data = np.genfromtxt(file_path, delimiter=',', dtype=None)
                    print "count of", file, ":",  data.shape
                    total += 1
                    count_array[data.shape] += 1

        print "Total All Families : ", total
        for i in range(0, len(count_array)):
            print i, ":", count_array[i]

    def generate_all_network_data_files(self):
        print ('grouping data')
        grouped_data = {}
        single_node_data = {}
        two_node_data = {}
        three_node_data = {}
        four_node_data = {}
        five_node_data = {}
        ten_node_data = {}
        twenty_node_data = {}
        for key, group in groupby(self.data, key=lambda x: x[0]):
            grouped_data[key] = [v[1] for v in group]
        for k, v in grouped_data.items():
            if len(list(filter(None, v))) == 1:
                single_node_data[k] = v
            if len(list(filter(None, v))) == 2:
                two_node_data[k] = v
            if len(list(filter(None, v))) == 3:
                three_node_data[k] = v
            if len(list(filter(None, v))) == 4:
                four_node_data[k] = v
            if len(list(filter(None, v))) == 5:
                five_node_data[k] = v
            if len(list(filter(None, v))) == 10:
                ten_node_data[k] = v
            if len(list(filter(None, v))) == 20:
                twenty_node_data[k] = v
            # Generate individual files for each network size
            self.write_all_network_file(len(list(filter(None, v))), v)
        print 'Total Single Nodes:', len(single_node_data)
        print 'Total Two Node networks:', len(two_node_data)
        print 'Total Three Node networks:', len(three_node_data)
        print 'Total Four Nodes:', len(four_node_data)
        print 'Total Five Node networks:', len(five_node_data)
        print 'Total Ten Node networks:', len(ten_node_data)
        print 'Total Twenty Node networks', len(twenty_node_data)

    def identify_network_group(self, data_array, column_count):
        if len([x for x in data_array if re.match(r'^evm.model.Contig', x)]) == column_count:
            self.write_all_network_file("all_geph", data_array.tolist())
        elif len([x for x in data_array if re.match(r'^evm.model.scaffold', x)]) == column_count:
            self.write_all_network_file("all_iso", data_array.tolist())
        elif len([x for x in data_array if re.match(r'^KOO', x)]) == column_count:
            self.write_all_network_file("all_chrys", data_array.tolist())
        elif len([x for x in data_array if re.match(r'^\d', x)]) == column_count:
            self.write_all_network_file("all_ehux", data_array.tolist())
        elif len([x for x in data_array if
                  re.match(r'^evm.model.Contig|evm.model.scaffold', x)]) == column_count:
            self.write_all_network_file("geph_iso", data_array.tolist())
        elif len([x for x in data_array if re.match(r'^evm.model.Contig|KOO', x)]) == column_count:
            self.write_all_network_file("geph_chry", data_array.tolist())
        elif len([x for x in data_array if re.match(r'^evm.model.Contig|\d', x)]) == column_count:
            self.write_all_network_file("geph_ehux", data_array.tolist())
        elif len([x for x in data_array if re.match(r'^evm.model.scaffold|KOO', x)]) == column_count:
            self.write_all_network_file("iso_chry", data_array.tolist())
        elif len([x for x in data_array if re.match(r'^evm.model.scaffold|\d', x)]) == column_count:
            self.write_all_network_file("iso_ehux", data_array.tolist())
        elif len([x for x in data_array if re.match(r'^KOO|\d', x)]) == column_count:
            self.write_all_network_file("chry_ehux", data_array.tolist())
        elif len([x for x in data_array if
                  re.match(r'^evm.model.Contig|evm.model.scaffold|\d', x)]) == column_count:
            self.write_all_network_file("geph_iso_ehux", data_array.tolist())
        elif len([x for x in data_array if
                  re.match(r'^evm.model.Contig|evm.model.scaffold|KOO', x)]) == column_count:
            self.write_all_network_file("geph_iso_chry", data_array.tolist())
        elif len([x for x in data_array if re.match(r'^\d|evm.model.scaffold|KOO', x)]) == column_count:
            self.write_all_network_file("iso_ehux_chry", data_array.tolist())
        elif len([x for x in data_array if re.match(r'^\d|evm.model.Contig|KOO', x)]) == column_count:
            self.write_all_network_file("geph_ehux_chry", data_array.tolist())
        else:
            self.write_all_network_file("species", data_array.tolist())

    def analyze_node_file(self, node_count):
        file_path = "../output/silix/node_data_%s.txt" % node_count
        if os.path.exists(file_path):
            self.data = np.genfromtxt(file_path, delimiter=',', dtype=str)
            print "data shape: ", self.data.shape, len(self.data.shape)
            column_count = int(node_count)
            if len(self.data.shape) > 1:
                rowcount = len(self.data)
            else:
                rowcount = 1
            print "Network of", column_count, " nodes "
            print "Total Networsk:", rowcount
            # print 'data', len(self.data), self.data.shape
            if rowcount > 1:
                for i in range(0, rowcount):
                    self.identify_network_group(self.data[i, :], column_count)
            else:
                #data_string = np.genfromtxt(file_path, delimiter=',', dtype=str)
                self.identify_network_group(self.data, column_count)
        print "--------------------------------------"

    def count_networks_by_size(self):
        for node_data in glob.glob("../output/silix/node_data*.txt"):
            print node_data
            print 'Total Node networks for ', node_data, len(np.genfromtxt(node_data, dtype=None))


if __name__ == "__main__":
    objSilix = SilixAnalysis()
    '''
    # Read All DATA FROM SILIX output and generate individual files
    objSilix.load_actual_data()
    objSilix.generate_all_network_data_files()
    #Loop thru all files and display size
    objSilix.count_networks_by_size()
    '''
    # read file with more than nodes
    for i in range(7, 233):
        objSilix.analyze_node_file(i)

    #objSilix.analyze_node_file("6")
    objSilix.get_results_statistics()

