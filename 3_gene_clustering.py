#
# When computing U, the number of comparisons equals the product of the number of values in group A times the number of values in group B. 
# If the null hypothesis is true, then the value of U should be about half that value. If the value of U is much smaller than that, the P 
# value will be small. The smallest possible value of U is zero. The largest possible value is half the product of the number of values in 
# group A times the number of values in group B.
#

import numpy as np 
import pandas as pd 
import sys
import operator
from sklearn.metrics import pairwise 
# import matplotlib.pyplot as plt 

from scipy import stats

check_dict = dict()

def find_similar_genes(sim_list):
	threshold = 0.55
	similar_genes = list()
	for i in range(len(sim_list)):
		if(sim_list[i] > threshold):
			similar_genes.append(i)
			check_dict[i] = True

	return similar_genes

def similarity_measure_form_cluster(X, gene_name_list):
	row, col = X.shape
	gene_cluster_dict = dict()
	cluster_index = 0

	sim_matrix = pairwise.cosine_similarity(X.T)

	x, y = sim_matrix.shape

	gene_similairty_cluster = dict()

	for i in range(x):
		similar_genes = list()
		similar_genes.append(i)
		check_dict[i] = True

		sim_list = sim_matrix[:, i]
		similar_genes.extend(find_similar_genes(sim_list))
		gene_similairty_cluster[i] = similar_genes

	cluster_index = 1
	for key in gene_similairty_cluster:
		fout = open('Cluster_' + str(cluster_index) + '.clstr', 'w')
		gene_list = gene_similairty_cluster[key]
		for i in range(len(gene_list)):
			fout.write(str(gene_name_list[gene_list[i]]) + '\n')

		print('Building cluster: ' + str(cluster_index))
		cluster_index += 1
		fout.close()

def pairwise_ttest(X):
	row, col = X.shape
	print(row, col)
	gene_cluster_dict = dict()
	cluster_index = 0

	while(row > 1):
		gene_name_list = X.columns.values
		i = 0
		gene1 = X[gene_name_list[i]]
		gene_list = list()
		gene_list.append(gene_name_list[i])

		for j in range(i + 1, row):
			gene2 = X[gene_name_list[j]]

			(t, p) = stats.ttest_ind(gene1, gene2, equal_var = False)
			if(t < 0 and p < 0.05):
				gene_list.append(gene_name_list[j])

		gene_cluster_dict[cluster_index] = gene_list
		cluster_index += 1

		print(gene_list)

		X.drop(gene_list, axis = 1, inplace = True)
		row, col = X.shape
		print(row, col)

	for key in gene_cluster_dict:
		print(key, gene_cluster_dict[key])

	for key in gene_cluster_dict:
		print(key, len(gene_cluster_dict[key]))

def gene_filter_modifying_dataset(X):
	samples, genes = X.shape

	gene_name_list = X.columns.values
	np_X = X.as_matrix()
	gene_snps_sum = np_X.sum(axis = 0)

	gene_snp_dict = dict()
	for i in range(len(gene_name_list)):
		gene_snp_dict[gene_name_list[i]] = gene_snps_sum[i]

	sorted_gene_snp_list_pair = sorted(gene_snp_dict.items(), key = operator.itemgetter(1), reverse = True)
	
	# setting threshold to cutoff genes without minimum number of snps across samples (1% of the samples has snps)
	snp_count_threshold = int((samples * 1)/ 100)

	# write the gene snp frequency
	fout = open('gene_snp_frequency_down.txt', 'w')

	qualified_gene_list = list()
	for i in range(len(sorted_gene_snp_list_pair)):
		if(sorted_gene_snp_list_pair[i][1] >= snp_count_threshold):
			qualified_gene_list.append(sorted_gene_snp_list_pair[i][0])
			fout.write(str(sorted_gene_snp_list_pair[i][0]) + '\t' + str(sorted_gene_snp_list_pair[i][1]) + '\n')
	fout.close()

	# rebuild the dataframe with only the qualified genes with enough snps across samples
	qualified_gene_df = pd.DataFrame()
	for i in range(len(qualified_gene_list)):
		qualified_gene_df[qualified_gene_list[i]] = X[qualified_gene_list[i]].values

	return qualified_gene_df
	

def main():
	data_file = sys.argv[1]		# read the all_data.txt
	X = pd.read_table(data_file, sep = '\t', header = 'infer')
	y = X['Cancer_type']
	X.drop('Cancer_type', axis = 1, inplace = True)

	# filter genes with enough number of snps across samples
	print('Filtering genes... ')
	filtered_df = gene_filter_modifying_dataset(X)

	# run cosine similarity test to cluster genes
	gene_name_list = filtered_df.columns.values
	# for i in range(len(gene_name_list)):
	# 	print(gene_name_list[i])
	print('Pairwise similarity measurement between genes... ')
	similarity_measure_form_cluster(filtered_df, gene_name_list)

if __name__ == "__main__":
	main()

