import glob
import os
import networkx as nx
import time
from Bio import SearchIO
import matplotlib.pyplot as plt

__author__ = 'FalguniT'


class BlastGraph(object):
    def __init__(self, evalue = 1e-06):
        self.blast_data_path = "../data/testblastdata/*.o"
        self.blast_connected_graph_statistics_file = ""
        self.blast_output_path = ""
        self.blast_experiment = ""
        self.generate_gml_files = True
        self.evalue = evalue  # 1e-30
        self.perc_identity = 40
        self.blast_graph = nx.Graph()
        self.initialize_variables()

    def initialize_variables(self):
        experiment_timestamp = time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime())
        filters = str(self.evalue) + str(self.perc_identity)
        self.blast_output_path = "../output/blast_connected_graphs_{}_{}".format(filters, experiment_timestamp)
        os.makedirs(self.blast_output_path)
        self.blast_connected_graph_statistics_file = "{}/blast_connected_graph_statistics.csv".format(
            self.blast_output_path)
        # Write header to generic statics file
        with open(self.blast_connected_graph_statistics_file, "a") as f_handle:
            f_handle.write("Graph no, Total count, Chry count, Geph count, Iso count, Ehux count, \n")

    def generate_blast_graph(self):
        evalue_filter = lambda hsp: hsp.evalue < self.evalue
        if self.perc_identity is not '':
            perc_filter = lambda hsp: hsp.ident_pct > self.perc_identity

        file_name = "{}/blast_graph.txt".format(self.blast_output_path)
        for blast_file in glob.glob(self.blast_data_path):
            #print("working on " + blast_file)
            # Parse the Blast file
            qresults = SearchIO.parse(blast_file, 'blast-tab', comments=True)
            for qresult in qresults:
                write_line = ""
                write_line += qresult.id + ":"
                # Add node to graph if not present
                if not self.blast_graph.has_node(qresult.id):
                    self.blast_graph.add_node(qresult.id)
                # Go to the Hit section of query
                for hit in qresult[:]:

                    # Check if Hit has min e-value
                    filtered_hit = hit.filter(evalue_filter)

                    if filtered_hit is not None and self.perc_identity is not '':
                        filtered_hit = filtered_hit.filter(perc_filter)
                        if filtered_hit is not None:
                            print 'filtered with perc ', filtered_hit.id

                    if filtered_hit is not None:
                        if not self.blast_graph.has_node(filtered_hit.id):
                            self.blast_graph.add_node(filtered_hit.id)
                            print filtered_hit.id
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
    #evalue range
    e_value_filter_range = [1e-06]
    for evalue in e_value_filter_range:
        print "Executing for experiment evalue : " + str(evalue)
        objblastgraph = BlastGraph(evalue)
        objblastgraph.generate_blast_graph()
        print("graph has %d nodes with %d edges"
              % (nx.number_of_nodes(objblastgraph.blast_graph), nx.number_of_edges(objblastgraph.blast_graph)))
        objblastgraph.generate_connected_component_graphs()
        print("connected components by function %s " % nx.number_connected_components(objblastgraph.blast_graph))
        filters = str(objblastgraph.evalue) + str(objblastgraph.perc_identity)
        with open("{}/details_{}.txt".format(objblastgraph.blast_output_path, filters), "a") as f_handle:
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
