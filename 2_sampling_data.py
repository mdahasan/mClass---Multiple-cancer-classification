import numpy as np 
import pandas as pd 
import sys
import os
import glob
import sklearn.metrics.pairwise 
# import matplotlib.pyplot as plt 
from scipy.stats import describe
from scipy import spatial

def read_cancer_data(path):
	X = pd.read_table(path, sep = '\t', header = 'infer')
	X.drop('Sample', axis = 1, inplace = True)
	y = X['Cancer_type']
	X.drop('Cancer_type', axis = 1, inplace = True)

	return X

def write_resampled_data(X, cancer_name):
	X.to_csv(str(cancer_name) + '_resampled.txt', sep = '\t', header = True, index = False)
 
def combine_all_resampled_dataset():
	resampled_file_list = glob.glob('*_resampled.txt')
	
	df_list = list()
	all_data = pd.DataFrame()

	for file in range(len(resampled_file_list)):
		X = pd.read_table(resampled_file_list[file], sep = '\t', header = 'infer')
		df_list.append(X)

	all_data = pd.concat(df_list)
	all_data.to_csv('all_data.txt', sep = '\t', header = True, index = False)

def main():
	# read each of the cancer data file and do subsampling by
	# eliminating similar samples based on cosine similarity
	data_directory = sys.argv[1]		# the directory that contains the individual cancer data
	current_path = os.getcwd()

	os.chdir(data_directory)
	cancer_data_list = glob.glob("*.txt")

	# iterate through each cancer dataset for similarity measurement
	for cancer in range(0, len(cancer_data_list)):
		cancer_data_path = os.path.join(data_directory, cancer_data_list[cancer])
		cancer_name = cancer_data_list[cancer].split('.txt')[0]
		X = read_cancer_data(cancer_data_path)
		np_X = X.as_matrix()
		gene_snps_sum = np_X.sum(axis = 1)	# count number of snps in all genes per sample

		sample_drop_list = list()
		for i in range(len(gene_snps_sum)):
			if(gene_snps_sum[i] < 5):
				sample_drop_list.append(i)
		X.drop(sample_drop_list, axis = 0, inplace = True)

		row, col = X.shape
		class_list = [cancer_name] * row
		X['Cancer_type'] = class_list

		print(str(cancer_name) + str(X.shape))
		write_resampled_data(X, cancer_name)

	combine_all_resampled_dataset()


if __name__ == "__main__":
	main()