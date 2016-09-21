import  unittest.test
import numpy as np
from src.GeneExpressionManager import GeneExpressionManager

__author__ = 'FalguniT'

class TestGeneExpressionManager(unittest.TestCase):
    def test_load_gene_expression_data_ehux(self):
        gene_expression_manager = GeneExpressionManager()
        gene_expression_manager.load_gene_expression_data()
        ehux_result_count = len(gene_expression_manager.ehux_gene_expression_data)
        expected_count = 24203
        self.assertEqual(ehux_result_count, expected_count)

    def test_load_gene_expression_data_for_geph(self):
        gene_expression_manager = GeneExpressionManager()
        gene_expression_manager.load_gene_expression_data()
        geph_result_count = len(gene_expression_manager.geph_gene_expression_data)
        expected_count = 24552
        self.assertEqual(geph_result_count, expected_count)

    def test_load_gene_expression_data_for_iso(self):
        gene_expression_manager = GeneExpressionManager()
        gene_expression_manager.load_gene_expression_data()
        iso_result_count = len(gene_expression_manager.iso_gene_expression_data)
        expected_count = 14703
        self.assertEqual(iso_result_count, expected_count)

    def test_map_sequence_to_gene_expression(self):
        gene_expression_manager = GeneExpressionManager()
        gene_expression_manager.load_gene_expression_data()
        # retrieve first row as expected row
        expected_row = gene_expression_manager.ehux_gene_expression_data[1][None, :]
        # pass expected row sequence no to function
        sequence_no = expected_row[0][0]
        result_row = gene_expression_manager.map_sequence_to_gene_expression(sequence_no)
        self.assertEqual(expected_row.tolist(), result_row.tolist())

    def test_map_sequence_to_gene_expression_geph(self):
        gene_expression_manager = GeneExpressionManager()
        gene_expression_manager.load_gene_expression_data()
        # retrieve first row as expected row
        expected_row = gene_expression_manager.geph_gene_expression_data[1][None, :]
        # pass expected row sequence no to function
        sequence_no = expected_row[0][0]
        result_row = gene_expression_manager.map_sequence_to_gene_expression(sequence_no)
        self.assertEqual(expected_row.tolist(), result_row.tolist())

    def test_map_sequence_to_gene_expression_iso(self):
        gene_expression_manager = GeneExpressionManager()
        gene_expression_manager.load_gene_expression_data()
        # retrieve first row as expected row
        expected_row = gene_expression_manager.iso_gene_expression_data[1][None, :]
        # pass expected row sequence no to function
        sequence_no = expected_row[0][0]
        result_row = gene_expression_manager.map_sequence_to_gene_expression(sequence_no)
        self.assertEqual(expected_row.tolist(), result_row.tolist())

    def test_map_sequence_to_gene_expression_chry(self):
        gene_expression_manager = GeneExpressionManager()
        gene_expression_manager.load_gene_expression_data()
        # retrieve first row as expected row
        expected_row = gene_expression_manager.iso_gene_expression_data[1][None, :]
        # pass expected row sequence no to function
        sequence_no = expected_row[0][0]
        result_row = gene_expression_manager.map_sequence_to_gene_expression(sequence_no)
        self.assertEqual(expected_row.tolist(), result_row.tolist())
