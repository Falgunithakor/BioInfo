import numpy as np
import Bio

print(Bio.__version__)
exit(0)

from src.FileManager import FileManager

__author__ = 'FalguniT'

#Ehux JGI Fasta file
Ehux_JGI_file_path = "../data/Ehux_JGI.fasta"
Ehux_JGI_data = FileManager.load_file(Ehux_JGI_file_path)


#Geph  Blast output file
Geph_file_path = "../data/Ehux_Geph_Blast_060916.txt.p1"
Geph_data = FileManager.load_file(Geph_file_path)
Geph_count = len(Geph_data)
print("geph count", Geph_count)

#ISO blast output file
ISO_file_path = "../data/Ehux_ISO_Blast_060916.txt.p1"
ISO_data = FileManager.load_file(ISO_file_path)
ISO_count = len(ISO_data)
print("ISO_ count", ISO_count)

#Strains 92A blast
strains_92A_file_path = "../data/Ehux_strains_92A_Blast_060816.txt.p1"
strains_92A_data = FileManager.load_file(strains_92A_file_path)
strains_92A_count = len(strains_92A_data)
print("strains_92A_ count", strains_92A_count)
#Strains EH2 blast
strains_EH2_file_path = "../data/Ehux_strains_EH2_Blast_060816.txt.p1"
strains_EH2_data = FileManager.load_file(strains_EH2_file_path)
strains_EH2_count = len(strains_EH2_data)
print("strains_EH2 count", strains_EH2_count)
#Strains Van 556
strains_van556_file_path = "../data/Ehux_strains_van556_Blast_060816.txt.p1"
strains_van556_data = FileManager.load_file(strains_van556_file_path)
strains_van556_count = len(strains_van556_data)
print("strains_van556 count", strains_van556_count)


#TODO - logic to check the file length of each and take smallest file
source_compare_data = {}

#considered ISO data smallest from data manually
print (len(ISO_data))
source_compare_data = ISO_data
source_data_count = ISO_count
#print(source_compare_data)

#check is same query exists in other blast output files
print(source_compare_data[1][0])
match_count = 0

for i in range(1,source_data_count):
    source_query_value = source_compare_data[i][0]
    #print(Geph_data[1:,0])
    #if source_query_value in Geph_data[1:,0]:
    if np.any(Geph_data[1:,0] == source_compare_data[i][0]) and np.any(strains_92A_data[1:,0] == source_compare_data[i][0]) \
        and np.any(strains_EH2_data[1:,0] == source_compare_data[i][0]) and np.any(strains_van556_data[1:,0] == source_compare_data[i][0]):
        print source_query_value
        match_count = match_count + 1
        #print ("Geph data")

print ("Total count :" , match_count)