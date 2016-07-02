import numpy as np
import time
from itertools import (takewhile,repeat)

__author__ = 'FalguniT'

class FileManager(object):
    @staticmethod
    def load_file(file_path):
        loaded_file = np.genfromtxt(file_path,delimiter='\t', dtype=None)
        return loaded_file

    @staticmethod
    def rawincount(filename):
        f = open(filename, 'rb')
        bufgen = takewhile(lambda x: x, (f.read(1024*1024) for _ in repeat(None)))
        return sum( buf.count(b'\n') for buf in bufgen )

    @staticmethod
    def read_csv_file(file_path):
        loaded_file = np.genfromtxt(file_path, delimiter=',', dtype=None)
        return loaded_file


    @staticmethod
    def create_output_file(feature_selection_algorithm ='DEBPSO', model = 'SVR'):
        timestamp = time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime())
        file_name = "../output/{}-{}-{}.csv".format('', '',timestamp)
        with open(file_name, "a") as f_handle:
            f_handle.write('Normalizer,feature_selection,model, Descriptors, Fitness, R2_Train, R2Pred_Validation, R2Pred_Test\n')
        return file_name

    @staticmethod
    def write_model_in_file(file_name, normalizer,feature_selection_algorithm, model , descriptor_ids,  fitness, r2_train, r2pred_validation, r2pred_test):
        with open(file_name,"a") as f_handle:
            f_handle.write( normalizer + ','+ feature_selection_algorithm + ','+model + ','+ str(descriptor_ids).replace(',','-') + ','+ str(fitness) + ','+ str(r2_train) + ','+ str(r2pred_validation) + ','+ str(r2pred_test)+ '\n')
