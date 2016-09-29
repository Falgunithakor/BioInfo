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
                 blast_data_path="../data/blastdata/*.o", by_sequence=[], by_desc="", by_go_no="",
                 by_go_name="", by_enzyme_code="", by_inter_pro_id=""):
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
        self.target_desc = by_desc
        self.target_GO_no = by_go_no
        self.target_GO_name = by_go_name
        self.target_enzyme_code = by_enzyme_code
        self.target_inter_pro_id = by_inter_pro_id
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
            # Parse each blast record
            for query_result in query_results:
                # retrieve Gene expression and gene annotation row for sequence
                gene_expression_row, gene_annotation_row = self.generate_sequence_gene_annotation_expression_mapping(
                    query_result.id)
                is_valid_sequence = self.apply_filtering(query_result, annotation_row=gene_annotation_row)

                if is_valid_sequence:
                    print query_result.id, is_valid_sequence, gene_annotation_row
                    self.generate_blast_graph(query_result, gene_expression_row, gene_annotation_row )

    def apply_filtering(self, query_result, annotation_row):
        # Filtering applies AND conditions for all filters set
        is_valid_sequence = True;
        if len(self.target_sequences) > 0 and query_result.id not in self.target_sequences:
            is_valid_sequence = False
        if len(self.target_desc) > 0 and self.target_desc not in annotation_row[1]:
            is_valid_sequence = False
        if len(self.target_GO_no) > 0 and self.target_GO_no != int(annotation_row[6]):
            is_valid_sequence = False
        if len(self.target_GO_name) > 0 and self.target_GO_name not in annotation_row[7]:
            is_valid_sequence = False
        if len(self.target_enzyme_code) > 0 and self.target_enzyme_code not in annotation_row[8]:
            is_valid_sequence = False
        if len(self.target_inter_pro_id) > 0 and self.target_inter_pro_id not in annotation_row[9]:
            is_valid_sequence = False
        # if any(ext in url_string for ext in extensionsToCheck):
        return is_valid_sequence

    def generate_blast_graph(self, query_result, gene_expression_row, gene_annotation_row):
        e_value_filter = lambda hsp: hsp.evalue < self.e_value
        # Add node to graph if not present
        if not self.blast_graph.has_node(query_result.id):
            self.add_graph_node(query_result.id, gene_expression_row, gene_annotation_row)
        # Go to the Hit section of query
        for hit in query_result[:]:
            # Check if Hit has min e-value
            filtered_hit = hit.filter(e_value_filter)
            if filtered_hit is not None:
                if not self.blast_graph.has_node(filtered_hit.id):
                    self.add_graph_node(filtered_hit.id, gene_expression_row, gene_annotation_row)
                # Add Edge between graph nodes
                self.blast_graph.add_edge(query_result.id, filtered_hit.id,
                                          evalue=str(filtered_hit.hsps[0].evalue),
                                          identpct=filtered_hit.hsps[0].ident_pct,
                                          mismatchno=filtered_hit.hsps[0].mismatch_num,
                                          aln=filtered_hit.hsps[0].aln,
                                          alnspn=filtered_hit.hsps[0].aln_span,
                                          gapopenno=filtered_hit.hsps[0].gapopen_num,
                                          bitscore=filtered_hit.hsps[0].bitscore)

    def add_graph_node(self, sequence_no, gene_expression_row, gene_annotation_row):
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

    def write_blast_graph_file(self):
        file_name = "{}/blast_graph.txt".format(self.blast_output_path)
        if len(nx.nodes(self.blast_graph)) > 0:
            with open(file_name, "wb") as f_handle:
                f_handle.write(
                    'Source' + ',' + str(
                        [value for value in self.blast_graph.node[nx.nodes(self.blast_graph)[0]].iterkeys()]).strip(
                        '[]') + ',')
                f_handle.write(
                    'Target' + ',' + str(
                        [value for value in self.blast_graph.node[nx.nodes(self.blast_graph)[0]].iterkeys()]).strip(
                        '[]') + ',')
                f_handle.write(
                    'evalue' + ',' + 'identpct' + ',' + 'mismatchno' + ',' + 'aln' + ',' + 'alnspn' + ',' + 'gapopenno' + ',' + 'bitscore' + '\n')

                for u, v, edata in self.blast_graph.edges(data=True):
                    f_handle.write(
                        str(u) + ',' + str([value for value in self.blast_graph.node[u].itervalues()]).strip(
                            '[]') + ',')
                    f_handle.write(
                        str(v) + ',' + str([value for value in self.blast_graph.node[v].itervalues()]).strip(
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
    # Instantiate Gene Annotaiton Manager object
    gene_annotation_manager = GeneAnnotationManager()
    gene_annotation_manager.load_gene_annotation_data()
    # initialize Filter criteria
    # By Target sequences
    # target_sequences = np.genfromtxt("../data/additiona_analysis/Ehux Lipid Metabolism Protein IDs.csv", delimiter='\n',dtype=str)
    # print target_sequences
    # By Description
    target_description = 'transferase'
    similarity_networks = SimilarityNetworks(gene_expression_manager, gene_annotation_manager,
                                             by_desc=target_description)
    similarity_networks.output_file_identifier = "_transferase_for_all_blast_data"
    similarity_networks.generate_blast_data()
    similarity_networks.write_blast_graph_file()
