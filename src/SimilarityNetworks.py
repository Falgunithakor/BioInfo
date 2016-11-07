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
                 blast_data_path="../data/blastdata/*.o", by_sequence=[]):
        self.blast_data_path = blast_data_path
        self.blast_output_path = ""
        self.blast_experiment = ""
        self.generate_gml_files = True
        self.e_value = evalue  # 1e-30
        self.blast_graph = nx.Graph()
        self.gene_expression_manager = gene_expression_manager
        self.gene_annotation_manager = gene_annotation_manager
        # Filtering parameters for selection of sequences
        self.target_sequences = by_sequence
        self.output_file_identifier = "";

    def initialize_variables(self):
        experiment_timestamp = time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime())
        filters = str(self.e_value) + '_' + self.output_file_identifier
        self.blast_output_path = "../output/blast_snw_{}_{}".format(filters, experiment_timestamp)
        os.makedirs(self.blast_output_path)

    def generate_blast_data(self):
        self.initialize_variables()
        for blast_file in glob.glob(self.blast_data_path):
            # Parse each Blast file
            query_results = SearchIO.parse(blast_file, 'blast-tab', comments=True)
            filtered_query_results = self.apply_filtering(query_results)
            # Parse each blast record
            for query_result in filtered_query_results:
                print query_result.id
                self.generate_blast_graph(query_result)

    def apply_filtering(self, query_results):
        # Filtering selects specific target sequences from results
        return (x for x in query_results if x.id in self.target_sequences)

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
                if not self.blast_graph.has_edge(query_result.id, filtered_hit.id):
                    self.blast_graph.add_edge(query_result.id, filtered_hit.id,
                                              evalue=str(filtered_hit.hsps[0].evalue),
                                              identpct=filtered_hit.hsps[0].ident_pct,
                                              mismatchno=filtered_hit.hsps[0].mismatch_num,
                                              aln=filtered_hit.hsps[0].aln,
                                              alnspn=filtered_hit.hsps[0].aln_span,
                                              gapopenno=filtered_hit.hsps[0].gapopen_num,
                                              bitscore=filtered_hit.hsps[0].bitscore)
                else:
                    if filtered_hit.hsps[0].evalue < self.blast_graph[query_result.id][filtered_hit.id]['evalue']:
                        self.blast_graph[query_result.id][filtered_hit.id]['evalue'] = str(filtered_hit.hsps[0].evalue)
                        self.blast_graph[query_result.id][filtered_hit.id]['identpct'] = str(
                            filtered_hit.hsps[0].ident_pct)
                        self.blast_graph[query_result.id][filtered_hit.id]['mismatchno'] = str(
                            filtered_hit.hsps[0].mismatch_num)
                        self.blast_graph[query_result.id][filtered_hit.id]['aln'] = str(filtered_hit.hsps[0].aln)
                        self.blast_graph[query_result.id][filtered_hit.id]['alnspn'] = str(
                            filtered_hit.hsps[0].aln_span)
                        self.blast_graph[query_result.id][filtered_hit.id]['gapopenno'] = str(
                            filtered_hit.hsps[0].gapopen_num)
                        self.blast_graph[query_result.id][filtered_hit.id]['bitscore'] = str(
                            filtered_hit.hsps[0].bitscore)
                print query_result.id, filtered_hit.hsps[0]

    def add_graph_node(self, sequence_no):
        # retrieve Gene expression and gene annotation row for sequence
        gene_expression_row, gene_annotation_row = self.generate_sequence_gene_annotation_expression_mapping(
            sequence_no)
        self.blast_graph.add_node(sequence_no,
                                  geintercept=gene_expression_row[0][1],
                                  gecacond=gene_expression_row[0][2],
                                  gespikecond=gene_expression_row[0][3],
                                  gadesc=gene_annotation_row[1],
                                  galength=gene_annotation_row[2],
                                  gahitscnt=gene_annotation_row[3],
                                  gaevalue=gene_annotation_row[4],
                                  gasimmean=gene_annotation_row[5],
                                  gagono=gene_annotation_row[6],
                                  gagoname=gene_annotation_row[7],
                                  gaenzymecodes=gene_annotation_row[8],
                                  gainterproids=gene_annotation_row[9])

    def generate_sequence_gene_annotation_expression_mapping(self, sequence_no):
        # retrieve gene expression and gene annotation data for query result id (query sequence)
        gene_expression_row = self.gene_expression_manager.map_sequence_to_gene_expression(sequence_no)
        gene_annotation_row = self.gene_annotation_manager.map_sequence_to_gene_annotation(sequence_no)
        return gene_expression_row, gene_annotation_row

    @staticmethod
    def get_species_name(name):
        if name.startswith("KOO"):
            return "Chry,"
        elif name.startswith("evm.model.Contig"):
            return "Geph,"
        elif name.startswith("evm.model.scaffold"):
            return "Iso,"
        else:
            return "Ehux,"

    def write_blast_graph_file(self):
        file_name = "{}/blast_graph.txt".format(self.blast_output_path)
        if len(nx.nodes(self.blast_graph)) > 0:
            with open(file_name, "wb") as f_handle:
                f_handle.write(
                    'Source' + ', Species,' + str(
                        [value for value in self.blast_graph.node[nx.nodes(self.blast_graph)[0]].iterkeys()]).strip(
                        '[]') + ',')
                f_handle.write(
                    'Target' + ', Species,' + str(
                        [value for value in self.blast_graph.node[nx.nodes(self.blast_graph)[0]].iterkeys()]).strip(
                        '[]') + ',')
                f_handle.write(
                    'evalue' + ',' + 'identpct' + ',' + 'mismatchno' + ',' + 'aln' + ',' + 'alnspn' + ',' + 'gapopenno' + ',' + 'bitscore' + '\n')

                for u, v, edata in self.blast_graph.edges(data=True):
                    f_handle.write(
                        str(u) + ',' + self.get_species_name(str(u)) + str(
                            [value for value in self.blast_graph.node[u].itervalues()]).strip(
                            '[]') + ',')
                    f_handle.write(
                        str(v) + ',' + self.get_species_name(str(v)) + str(
                            [value for value in self.blast_graph.node[v].itervalues()]).strip(
                            '[]') + ',')
                    f_handle.write(
                        str(edata['evalue']) + ',' + str(edata['identpct']) + ',' + str(edata['mismatchno']) + ',' +
                        str(edata['aln']) + ',' + str(edata['alnspn']) + ',' +
                        str(edata['gapopenno']) + ',' + str(edata['bitscore']) + '\n')
                    # if self.generate_gml_files:
                    # file_name = "{}/blast_graph.gml".format(self.blast_output_path)
                    # with open(file_name, "a") as f_handle:
                    # nx.write_gml(self.blast_graph, f_handle)


if __name__ == "__main__":
    # Instantiate Gene Expression Manager object
    gene_expression_manager = GeneExpressionManager()
    gene_expression_manager.load_gene_expression_data()
    # Instantiate Gene Annotation Manager object
    gene_annotation_manager = GeneAnnotationManager()
    gene_annotation_manager.load_gene_annotation_data()
    # Filter By Target sequences
    target_sequences = np.genfromtxt("../data/protein_sequence_files/elong.txt", delimiter='\n', dtype=str)
    print target_sequences
    print 'Total Target Sequences:', target_sequences.shape
    similarity_networks = SimilarityNetworks(gene_expression_manager, gene_annotation_manager, evalue=1e-05,
                                             by_sequence=target_sequences)
    similarity_networks.output_file_identifier = "_elong_for_all_blast_data"
    similarity_networks.generate_blast_data()
    similarity_networks.write_blast_graph_file()
