import glob
import os
from Bio import SearchIO

import numpy as np
from itertools import groupby
import pandas as pd
import re

__author__ = 'FalguniT'


class SilixNetworkChart(object):
    def __init__(self):
        self.network_data = []
        self.silix_nw_file_path = "../output/silix/node_data_232.txt"
        self.blast_data_path = "../data/testblastdata/*.o"


    def load_network_data(self):
        self.network_data = np.genfromtxt(self.silix_nw_file_path, dtype=str).tolist()
        print self.network_data

    def retrieve_blast_data(self):
        for blast_file in glob.glob(self.blast_data_path):
            print(blast_file)
            print self.network_data
            qresults = SearchIO.parse(blast_file, 'blast-tab', comments=True)
            for qresult in qresults:
                if(qresult.id in self.network_data):
                    print qresult.id
                    #print qresult.hit_filter(lambda hit.description.startswith('Homo sapiens'))

if __name__ == "__main__":
    objSilixNetwork = SilixNetworkChart()
    objSilixNetwork.load_network_data()
    objSilixNetwork.retrieve_blast_data()
