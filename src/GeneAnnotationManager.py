import re
import numpy as np

__author__ = 'FalguniT'


class GeneAnnotationManager(object):
    def __init__(self):
        self.geph_gene_annotation_data = []
        self.ehux_gene_annotation_data = []
        self.iso_gene_annotation_data = []
        # INPUT Gene Annotation Files
        self.geph_gene_annotation_file_path = "../data/gene_annotations/geph_final_b2g.txt"
        self.ehux_gene_annotation_file_path = "../data/gene_annotations/Ehux_JGI_prot_b2g.txt"
        self.iso_gene_annotation_path = "../data/gene_annotations/iso_final_b2g.txt"

    def load_gene_annotation_data(self):
        self.geph_gene_annotation_data = \
            np.genfromtxt(self.geph_gene_annotation_file_path, delimiter='\n', dtype=None, skip_header=1)
        self.ehux_gene_annotation_data = \
            np.genfromtxt(self.ehux_gene_annotation_file_path, delimiter='\n', dtype=None, skip_header=1)
        self.iso_gene_annotation_data = \
            np.genfromtxt(self.iso_gene_annotation_path, delimiter='\n', dtype=None, skip_header=1)

    def map_sequence_to_gene_annotation(self, sequence_no):
        gene_annotation_data = []
        gene_annotation_row = []
        # Retrieve annotation data based on sequence no
        if re.match(r'^\d', sequence_no):
            gene_annotation_data = self.ehux_gene_annotation_data
        elif re.match(r'^evm.model.Contig', sequence_no):
            gene_annotation_data = self.geph_gene_annotation_data
        elif re.match(r'^evm.model.scaffold', sequence_no):
            gene_annotation_data = self.iso_gene_annotation_data
        # Read gene annotation line
        if len(gene_annotation_data) > 0:
            gene_annotation_row = gene_annotation_data[
                [i for i, item in enumerate(gene_annotation_data) if str(item).startswith(sequence_no + '\t')]]
            if len(gene_annotation_row) > 0:
                gene_annotation_row = re.split(r'\t', gene_annotation_row[0])
                gene_annotation_row = [w.replace(',', ';') for w in gene_annotation_row]
        # gene annotation not found then set default string
        if len(gene_annotation_row) is 0:
            gene_annotation_row = np.full((10), '-', dtype=np.str).tolist()
        return gene_annotation_row
