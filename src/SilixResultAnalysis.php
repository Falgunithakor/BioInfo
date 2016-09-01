from collections import Counter
import glob

import numpy as np
from itertools import groupby

from sklearn.feature_selection.univariate_selection import f_classif

__author__ = 'FalguniT'


class SilixResultAnalysis(object):
    def __init__(self):
        self.silix_output_file_path = "../output/silix/output.txt"
        self.data = []

    @staticmethod
    def write_network_file(file_path, networkdata):
        with open(file_path, 'w') as f_handle:
            for key, value in sorted(networkdata.items()):
                f_handle.write(str(value) + '\n')

    def write_all_network_file(self, node_count, value):
        file_path = '../output/silix/node_data_%s.txt' % node_count
        with open(file_path, 'a') as f_handle:
            f_handle.write(str(value) + '\n')

    def load_actual_data(self):
        self.data = np.genfromtxt(self.silix_output_file_path, delimiter='\t', dtype=None)
        # print self.data

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

    def analyze_node_file(self, node_count):
        file_path="../output/silix/node_data_%s.txt" % node_count
        self.data = np.genfromtxt(file_path, delimiter=',', dtype=None)

        print "Total", len(self.data)
        for i in range(0, 10):
            print 'data', str(self.data[i][1]) , type(self.data[i][1])
            #print 'lamda', len(list(filter(lambda col: col.startswith('4'), self.data[i])))
            #print 'lem', len([x for x in self.data[i] if x.startswith('evm')])
            if 'Contig' in self.data[i][1]:
                print 'Contig'
            elif self.data[i][1].isdigit():
                print 'ehux'


    def count_networks_by_size(self):
        for node_data in glob.glob("../output/silix/node_data*.txt"):
            print 'Total Node networks for ',node_data , len(np.genfromtxt(node_data, delimiter=',', dtype=None))

if __name__ == "__main__":
    objSilix = SilixResultAnalysis()

    #Read All DATA FROM SILIX output and generate individual files
    '''
    objSilix.load_actual_data()
    objSilix.generate_all_network_data_files()
    objSilix.count_networks_by_size()
    '''
    objSilix.count_networks_by_size()
    #read file with two nodes
    objSilix.analyze_node_file("2")
