import glob
import os

import networkx as nx
import time
from Bio import SearchIO
import matplotlib.pyplot as plt


class BlastGraph(object):
    def __init__(self):
        self.blast_data_path = "../data/blastdata/*.o"
        self.blast_output_path = ""
        self.blast_experiment = ""
        self.generate_gml_files = False;
        self.evalue = 0.000000000001 #1e-12
        self.blast_graph = nx.Graph()
        self.initialize_variables()

    def initialize_variables(self):
        experiment_timestamp = time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime())
        self.blast_output_path = "../output/blast_connected_graphs_{}".format(experiment_timestamp)
        os.makedirs(self.blast_output_path)

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
        plt.figure(1, figsize=(25, 25))
        for sub_graph in connected_component_sub_graphs:
            index += 1
            file_name = "{}/graph_{}.txt".format(self.blast_output_path, index)
            with open(file_name, "a") as f_handle:
                f_handle.write(str(nx.nodes(sub_graph)) + ':' + str(nx.edges(sub_graph)) + '\n')
            # Write GML files
            if self.generate_gml_files:
                file_name = "{}/graph_{}.gml".format(self.blast_output_path,index)
                with open(file_name, "a") as f_handle:
                    nx.write_gml(sub_graph, f_handle)


if __name__ == "__main__":
    objblastgraph = BlastGraph()
    objblastgraph.generate_blast_graph()
    print("graph has %d nodes with %d edges" \
          % (nx.number_of_nodes(objblastgraph.blast_graph), nx.number_of_edges(objblastgraph.blast_graph)))
    objblastgraph.generate_connected_component_graphs()
    print(nx.number_connected_components(objblastgraph.blast_graph), "connected components by function")
