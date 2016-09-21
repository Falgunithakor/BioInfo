import unittest.test
import numpy as np
from src.GeneExpressionManager import GeneExpressionManager
from src.SimilarityNetworks import SimilarityNetworks

__author__ = 'FalguniT'

class TestSimilarityNetworks(unittest.TestCase):
    def test_load_blast_data(self):
        similarity_networks = SimilarityNetworks()
        similarity_networks.generate_blast_data()
        similarity_networks.write_blast_graph_file()
