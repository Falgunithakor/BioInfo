import glob
import os
import networkx as nx
import time
from Bio import SearchIO
import matplotlib.pyplot as plt

__author__ = 'FalguniT'


class BlastGraph(object):
    def __init__(self):
        self.blast_data_path = "../data/blastdata/*.o"
        self.blast_connected_graph_statistics_file = ""
        self.blast_output_path = ""
        self.blast_experiment = ""
        self.generate_gml_files = False
        self.evalue = 0.000000000000000000000000000001  # 1e-30
        self.blast_graph = nx.Graph()
        self.initialize_variables()

    def initialize_variables(self):
        experiment_timestamp = time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime())
        self.blast_output_path = "../output/blast_connected_graphs_{}_{}".format(self.evalue, experiment_timestamp)
        os.makedirs(self.blast_output_path)
        self.blast_connected_graph_statistics_file = "{}/blast_connected_graph_statistics.csv".format(
            self.blast_output_path)
        # Write header to generic statics file
        with open(self.blast_connected_graph_statistics_file, "a") as f_handle:
            f_handle.write("Graph no, Total count, Chry count, Geph count, Iso count, Ehux count, \n")

    def generate_blast_graph(self):
        evalue_filter = lambda hsp: hsp.evalue < self.evalue
        file_name = "{}/blast_graph.txt".format(self.blast_output_path)
        for blast_file in glob.glob(self.blast_data_path):
            print("working on " + blast_file)
            # Parse the Blast file
            qresults = SearchIO.parse(blast_file, 'blast-tab', comments=True)
            for qresult in qresults:
                write_line = ""
                write_line += qresult.id + ":"
                # Go to the Hit section of query
                for hit in qresult[:]:
                    if not self.blast_graph.has_node(qresult.id):
                        self.blast_graph.add_node(qresult.id)
                    # Check if Hit has min value
                    filtered_hit = hit.filter(evalue_filter)
                    if filtered_hit is not None:
                        if not self.blast_graph.has_node(filtered_hit.id):
                            self.blast_graph.add_node(filtered_hit.id)
                        # Add Edge between graph nodes
                        self.blast_graph.add_edge(qresult.id, filtered_hit.id)
                        write_line += filtered_hit.id + ","
                if write_line != "":
                    with open(file_name, "a") as f_handle:
                        f_handle.write(write_line + '\n')

        # Write GML files
        if self.generate_gml_files:
            file_name = "{}/blast_graph.gml".format(self.blast_output_path)
            with open(file_name, "a") as f_handle:
                nx.write_gml(self.blast_graph, f_handle)

    def generate_connected_component_graphs(self):
        # Get Connected Component sub graphs
        connected_component_sub_graphs = nx.connected_component_subgraphs(self.blast_graph)
        index = 0
        chrysochromulina_count = 0
        geph_count = 0
        iso_count = 0
        ehux_count = 0
        plt.figure(1, figsize=(25, 25))
        for sub_graph in connected_component_sub_graphs:
            chrysochromulina_count = 0
            geph_count = 0
            iso_count = 0
            ehux_count = 0
            index += 1
            file_name = "{}/graph_{}.txt".format(self.blast_output_path, index)
            with open(file_name, "a") as f_handle:
                f_handle.write(str(nx.nodes(sub_graph)) + ':' + str(nx.edges(sub_graph)) + '\n')

            for node in nx.nodes(sub_graph):
                if node.startswith("KOO"):
                    chrysochromulina_count += 1
                elif node.startswith("evm.model.Contig"):
                    geph_count += 1
                elif node.startswith("evm.model.scaffold"):
                    iso_count += 1
                else:
                    ehux_count += 1

            # generate excel file for histogram graph
            with open(self.blast_connected_graph_statistics_file, "a") as f_stat:
                f_stat.write(str(index) + ',' + str(str(nx.number_of_nodes(sub_graph))) + ',' + str(
                    chrysochromulina_count) + ',' + str(geph_count) + ',' + str(iso_count) + ',' + str(
                    ehux_count) + '\n')

                # Write GML files
            if self.generate_gml_files:
                file_name = "{}/graph_{}.gml".format(self.blast_output_path, index)
                with open(file_name, "a") as f_handle_gml:
                    nx.write_gml(sub_graph, f_handle_gml)


if __name__ == "__main__":
    objblastgraph = BlastGraph()
    objblastgraph.generate_blast_graph()
    print("graph has %d nodes with %d edges"
          % (nx.number_of_nodes(objblastgraph.blast_graph), nx.number_of_edges(objblastgraph.blast_graph)))
    objblastgraph.generate_connected_component_graphs()
    print("connected components by function %s " % nx.number_connected_components(objblastgraph.blast_graph))
    with open("{}/details.txt".format(objblastgraph.blast_output_path), "a") as f_handle:
        chrysochromulina_count = 0
        geph_count = 0
        iso_count = 0
        ehux_count = 0
        for node in nx.nodes(objblastgraph.blast_graph):
            if node.startswith("KOO"):
                chrysochromulina_count += 1
            elif node.startswith("evm.model.Contig"):
                geph_count += 1
            elif node.startswith("evm.model.scaffold"):
                iso_count += 1
            else:
                ehux_count += 1

        f_handle.write("Experiment for e-value %s" % objblastgraph.evalue + '\n')
        f_handle.write("Total %d nodes with %d edges"
                       % (nx.number_of_nodes(objblastgraph.blast_graph),
                          nx.number_of_edges(objblastgraph.blast_graph)) + '\n')
        f_handle.write("Chry count %d"
                       % (chrysochromulina_count) + '\n')
        f_handle.write("Geph count %d"
                       % (geph_count) + '\n')
        f_handle.write("ISO count %d"
                       % (iso_count) + '\n')
        f_handle.write("Ehux count %d"
                       % (ehux_count) + '\n')

        f_handle.write(
            "Total %s connected components " % nx.number_connected_components(objblastgraph.blast_graph))
