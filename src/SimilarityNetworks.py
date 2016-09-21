import glob
import os
import numpy as np
import networkx as nx

import time
from Bio import SearchIO
import matplotlib.pyplot as plt

from src.GeneAnnotationManager import GeneAnnotationManager
from src.GeneExpressionManager import GeneExpressionManager

__author__ = 'FalguniT'


class SimilarityNetworks(object):
    def __init__(self, gene_expression_manager, gene_annotation_manager, evalue=1e-05,
                 blast_data_path="../data/blastdata/*.o"):
        self.blast_data_path = blast_data_path
        self.target_sequences = []
        self.blast_output_path = ""
        self.blast_experiment = ""
        self.generate_gml_files = True
        self.e_value = evalue  # 1e-30
        self.blast_graph = nx.Graph()
        self.gene_expression_manager = gene_expression_manager
        self.gene_annotation_manager = gene_annotation_manager

    def initialize_variables(self):
        experiment_timestamp = time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime())
        filters = str(self.e_value)
        self.blast_output_path = "../output/blast_similaritynetworks_{}_{}".format(filters, experiment_timestamp)
        os.makedirs(self.blast_output_path)

    def generate_blast_data(self, target_sequences=[]):
        self.initialize_variables()
        for blast_file in glob.glob(self.blast_data_path):
            # Parse each Blast file
            query_results = SearchIO.parse(blast_file, 'blast-tab', comments=True)
            for query_result in query_results:
                if query_result.id in target_sequences:
                    self.generate_blast_graph(query_result)

    def generate_blast_graph(self, query_result):
        e_value_filter = lambda hsp: hsp.evalue < self.e_value
        # Add node to graph if not present
        if not self.blast_graph.has_node(query_result.id):
            self.add_graph_node(query_result.id)
        # Go to the Hit section of query
        for hit in query_result[:]:
            # Check if Hit has min e-value
            filtered_hit = hit.filter(e_value_filter)
            if filtered_hit is not None:
                if not self.blast_graph.has_node(filtered_hit.id):
                    self.add_graph_node(filtered_hit.id)
                # Add Edge between graph nodes
                self.blast_graph.add_edge(query_result.id, filtered_hit.id,
                                          evalue=str(filtered_hit.hsps[0].evalue),
                                          identpct=filtered_hit.hsps[0].ident_pct,
                                          mismatchno=filtered_hit.hsps[0].mismatch_num,
                                          aln=filtered_hit.hsps[0].aln,
                                          alnspn=filtered_hit.hsps[0].aln_span,
                                          gapopenno=filtered_hit.hsps[0].gapopen_num,
                                          bitscore=filtered_hit.hsps[0].bitscore
                                          )

    def add_graph_node(self, sequence_no):
        # retrieve gene expression and gene annotation data for query result id (query sequence)
        gene_expression_row = self.gene_expression_manager.map_sequence_to_gene_expression(sequence_no)
        gene_annotataion_row = self.gene_annotation_manager.map_sequence_to_gene_annotation(sequence_no)
        self.blast_graph.add_node(sequence_no,
                                  geintercept=gene_expression_row[0][1],
                                  gecacond=gene_expression_row[0][2],
                                  gespikecond=gene_expression_row[0][3],
                                  gadesc=gene_annotataion_row[1],
                                  galength=gene_annotataion_row[2],
                                  gahitscnt=gene_annotataion_row[3],
                                  gaevalue=gene_annotataion_row[4],
                                  gasimmean=gene_annotataion_row[5],
                                  gagono=gene_annotataion_row[6],
                                  gagoname=gene_annotataion_row[7],
                                  gaenzymecodes=gene_annotataion_row[8],
                                  gainterproids=gene_annotataion_row[9]
                                  )

    def generate_sequence_gene_annotation_mapping(self):
        pass

    def generate_sequence_gene_expression_mappings(self):
        pass

    def write_blast_graph_file(self):
        file_name = "{}/blast_graph.txt".format(self.blast_output_path)
        with open(file_name, "wb") as f_handle:
            for u, v, edata in self.blast_graph.edges(data=True):
                print u, v,
                f_handle.write(
                    str(u) + ',' + str([value for value in self.blast_graph.node[u].itervalues()]).strip('[]') + ',')
                f_handle.write(
                    str(v) + ',' + str([value for value in self.blast_graph.node[v].itervalues()]).strip('[]') + ',')
                f_handle.write(
                    str(edata['evalue']) + ',' + str(edata['identpct']) + ',' + str(edata['mismatchno']) + ',' +
                    str(edata['aln']) + ',' + str(edata['alnspn']) + ',' +
                    str(edata['gapopenno']) + ',' + str(edata['bitscore']) + '\n')
                # if self.generate_gml_files:
                # file_name = "{}/blast_graph.gml".format(self.blast_output_path)
                # with open(file_name, "a") as f_handle:
                # nx.write_gml(self.blast_graph, f_handle)

    def blast_result_statistics(self):
        filters = str(self.e_value)
        with open("{}/details_{}.txt".format(self.blast_output_path, filters), "a") as f_handle:
            chrysochromulina_count = 0
            geph_count = 0
            iso_count = 0
            ehux_count = 0
            for node in nx.nodes(self.blast_graph):
                if node.startswith("KOO"):
                    chrysochromulina_count += 1
                elif node.startswith("evm.model.Contig"):
                    geph_count += 1
                elif node.startswith("evm.model.scaffold"):
                    iso_count += 1
                else:
                    ehux_count += 1
            f_handle.write("Experiment for e-value %s" % self.e_value + '\n')
            f_handle.write("Total %d nodes with %d edges"
                           % (nx.number_of_nodes(self.blast_graph),
                              nx.number_of_edges(self.blast_graph)) + '\n')
            f_handle.write("Chry count %d"
                           % (chrysochromulina_count) + '\n')
            f_handle.write("Geph count %d"
                           % (geph_count) + '\n')
            f_handle.write("ISO count %d"
                           % (iso_count) + '\n')
            f_handle.write("Ehux count %d"
                           % (ehux_count) + '\n')


if __name__ == "__main__":
    # Instantiate Gene Expression Manager object
    gene_expression_manager = GeneExpressionManager()
    gene_expression_manager.load_gene_expression_data()
    # Instantiate Gene Annotaiton Manager object
    gene_annotation_manager = GeneAnnotationManager()
    gene_annotation_manager.load_gene_annotation_data()

    similarity_networks = SimilarityNetworks(gene_expression_manager, gene_annotation_manager)
    target_sequences = ['432589', 'evm.model.Contig1207.8', 'evm.model.Contig1623.9', 'evm.model.scaffold_26.48']
    similarity_networks.generate_blast_data(target_sequences)

    similarity_networks.write_blast_graph_file()
