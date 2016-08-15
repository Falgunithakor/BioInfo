import random
import time

__author__ = 'FalguniT'

import numpy as np
import scipy.sparse as sps
import networkx as nx
from networkx import nx_pydot
from Bio import SearchIO
import matplotlib.pyplot as plt


def main():
    # read first Blast file
    Ehux_Chry_FILE_PATH = "../data/blastdata/blastp_Ehux_Chry.o";
    Ehux_ISO_FILE_PATH = "../data/blastdata/blastp_Ehux_Iso.o";
    Ehux_Geph_FILE_PATH = "../data/blastdata/blastp_Ehux_Geph.o";
    Geph_Chry_FILE_PATH = "../data/blastdata/blastp_Geph_Chry.o";
    Geph_Ehux_FILE_PATH = "../data/blastdata/blastp_Geph_Ehux.o";
    Geph_Iso_FILE_PATH = "../data/blastdata/blastp_Geph_Iso.o";
    Iso_Chry_FILE_PATH = "../data/blastdata/blastp_Iso_Chry.o";
    Iso_Ehux_FILE_PATH = "../data/blastdata/blastp_Iso_Ehux.o";
    Iso_Geph_FILE_PATH = "../data/blastdata/blastp_Iso_Geph.o";

    blast_list = [Ehux_ISO_FILE_PATH,Ehux_Chry_FILE_PATH, Ehux_Geph_FILE_PATH,
                  Geph_Chry_FILE_PATH, Geph_Ehux_FILE_PATH,Geph_Iso_FILE_PATH,
                  Iso_Chry_FILE_PATH,Iso_Ehux_FILE_PATH, Iso_Geph_FILE_PATH]
    g = nx.Graph()
    data = []
    evalue_filter = lambda hsp: hsp.evalue < 1e-10

    for blast_file in blast_list:
        print(blast_file)
        result_count = 0
        qresults = SearchIO.parse(blast_file, 'blast-tab', comments=True)
        for qresult in qresults:
            result_count += 1
            for hit in qresult[:3]:
                if not g.has_node(qresult.id):
                    print(qresult.id)
                    g.add_node(qresult.id)
                filtered_hit = hit.filter(evalue_filter)
                if filtered_hit is not None:
                    if not g.has_node(filtered_hit.id):
                        g.add_node(filtered_hit.id)
                    g.add_edge(qresult.id, filtered_hit.id)

    print("graph has %d nodes with %d edges" \
          % (nx.number_of_nodes(g), nx.number_of_edges(g)))
    print(nx.number_connected_components(g), "connected components")


    C = nx.connected_component_subgraphs(g)
    timestamp = time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime())
    file_name = "../output/blast_connected_graphs_{}.csv".format(timestamp)
    with open(file_name, "a") as f_handle:
        for x in C:
            print("graph has %d nodes with %d edges" % (nx.number_of_nodes(x), nx.number_of_edges(x)))
            f_handle.write( str(nx.degree(x))+ '\n' )
            f_handle.write(str(nx.edges(x)) + '\n')

    plt.figure(1, figsize=(25, 25))
    for x in C:
        c = [random.random()]*nx.number_of_nodes(x) # random color...
        nx.draw(x,
             pos=None,
             node_size=40,
             node_color=c,
             vmin=0.0,
             vmax=1.0,
             with_labels=True
             )
    plt.show()




if __name__ == "__main__":
    main()
