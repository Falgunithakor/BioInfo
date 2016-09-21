import unittest.test
import numpy as np
import re
from src.GeneAnnotationManager import GeneAnnotationManager

__author__ = 'FalguniT'


class TestGeneAnnotationManager(unittest.TestCase):
    def test_load_gene_annotation_data_ehux(self):
        gene_annotation_manager = GeneAnnotationManager()
        gene_annotation_manager.load_gene_annotation_data()
        ehux_result_count = len(gene_annotation_manager.ehux_gene_annotation_data)
        expected_count = 30448
        self.assertEqual(ehux_result_count, expected_count)

    def test_load_gene_annotation_data_geph(self):
        gene_annotation_manager = GeneAnnotationManager()
        gene_annotation_manager.load_gene_annotation_data()
        geph_result_count = len(gene_annotation_manager.geph_gene_annotation_data)
        expected_count = 52679
        self.assertEqual(geph_result_count, expected_count)

    def test_load_gene_annotation_data_iso(self):
        gene_annotation_manager = GeneAnnotationManager()
        gene_annotation_manager.load_gene_annotation_data()
        iso_result_count = len(gene_annotation_manager.iso_gene_annotation_data)
        expected_count = 18712
        self.assertEqual(iso_result_count, expected_count)

    def test_map_sequence_to_gene_annotation_ehux(self):
        gene_annotation_manager = GeneAnnotationManager()
        gene_annotation_manager.load_gene_annotation_data()
        expected_row = gene_annotation_manager.ehux_gene_annotation_data[1]
        expected_row = re.split(r'\t+', expected_row.rstrip('\t'))
        sequence_no = expected_row[0]
        expected_annotation_desc = expected_row[1]
        result_row = gene_annotation_manager.map_sequence_to_gene_annotation(sequence_no)
        self.assertEqual(expected_annotation_desc, result_row[1])

    def test_map_sequence_to_gene_annotation_geph(self):
        gene_annotation_manager = GeneAnnotationManager()
        gene_annotation_manager.load_gene_annotation_data()
        expected_row = gene_annotation_manager.geph_gene_annotation_data[1]
        expected_row = re.split(r'\t+', expected_row.rstrip('\t'))
        sequence_no = expected_row[0]
        expected_annotation_desc = expected_row[1]
        result_row = gene_annotation_manager.map_sequence_to_gene_annotation(sequence_no)
        self.assertEqual(expected_annotation_desc, result_row[1])

    def test_map_sequence_to_gene_annotation_iso(self):
        gene_annotation_manager = GeneAnnotationManager()
        gene_annotation_manager.load_gene_annotation_data()
        expected_row = gene_annotation_manager.iso_gene_annotation_data[1]
        expected_row = re.split(r'\t+', expected_row.rstrip('\t'))
        sequence_no = expected_row[0]
        expected_annotation_desc = expected_row[1]
        result_row = gene_annotation_manager.map_sequence_to_gene_annotation(sequence_no)
        self.assertEqual(expected_annotation_desc, result_row[1])

