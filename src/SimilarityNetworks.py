import glob
import os
import networkx as nx
import time
from Bio import SearchIO
import matplotlib.pyplot as plt

__author__ = 'FalguniT'


class SimilarityNetworks(object):
    def __init__(self, evalue = 1e-05, blast_data_path = "../data/testblastdata/*.o"):
        self.blast_data_path = blast_data_path
        self.blast_output_path = ""
        self.blast_experiment = ""
        self.generate_gml_files = True
        self.e_value = evalue  # 1e-30
        self.blast_graph = nx.Graph()


    def initialize_variables(self):
        experiment_timestamp = time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime())
        filters = str(self.e_value)
        self.blast_output_path = "../output/blast_similaritynetworks_{}_{}".format(filters, experiment_timestamp)
        os.makedirs(self.blast_output_path)

    def generate_blast_data(self):
        self.initialize_variables()
        for blast_file in glob.glob(self.blast_data_path):
            # Parse each Blast file
            query_results = SearchIO.parse(blast_file, 'blast-tab', comments=True)
            for query_results in query_results:
                self.generate_blast_graph(query_results)


    def generate_blast_graph(self, query_result):
        e_value_filter = lambda hsp: hsp.evalue < self.e_value
        # Add node to graph if not present
        if not self.blast_graph.has_node(query_result.id):
            self.blast_graph.add_node(query_result.id, )
        # Go to the Hit section of query
        for hit in query_result[:]:
            # Check if Hit has min e-value
            filtered_hit = hit.filter(e_value_filter)
            if filtered_hit is not None:
                if not self.blast_graph.has_node(filtered_hit.id):
                    self.blast_graph.add_node(filtered_hit.id)
                # Add Edge between graph nodes
                self.blast_graph.add_edge(query_result.id, filtered_hit.id, evalue=filtered_hit.hsps[0].evalue)

    def generate_sequence_gene_annotation_mapping(self):
        pass

    def generate_sequence_gene_expression_mappings(self):
        pass

    def write_blast_graph_file(self):
        file_name = "{}/blast_graph.txt".format(self.blast_output_path)
        with open(file_name, "wb") as f_handle:
            nx.write_edgelist(self.blast_graph, f_handle, delimiter=',')

        if self.generate_gml_files:
            file_name = "{}/blast_graph.gml".format(self.blast_output_path)
            with open(file_name, "a") as f_handle:
                nx.write_gml(self.blast_graph, f_handle)

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
