import numpy as np 
import pandas as pd 
import sys 
import glob
import operator
# import matplotlib.pyplot as plt 
from sklearn import preprocessing as pp  

gene_MI_dict = dict()

def shannon_entropy(c):

	c_normalized = c / float(np.sum(c))
	c_normalized = c_normalized[np.nonzero(c_normalized)]
	H = -sum(c_normalized * np.log2(c_normalized))
	return H

def calc_MI(X, Y, bins):

	c_XY = np.histogram2d(X, Y, bins)[0]
	c_X = np.histogram(X, bins)[0]
	c_Y = np.histogram(Y, bins)[0]

	H_X = shannon_entropy(c_X)
	H_Y = shannon_entropy(c_Y)
	H_XY = shannon_entropy(c_XY)

	MI = H_X + H_Y - H_XY

	return MI

def mutual_info(clstr_genes, original_D, bins):
	cancer_label = original_D['Cancer_type'].tolist()
	le = pp.LabelEncoder()
	le.fit(cancer_label)
	encoded_label = le.transform(cancer_label)

	important_gene_list = list()
	for i in range(len(clstr_genes)):
		gene = clstr_genes[i]
		gene_snps = original_D[gene].tolist()

		# if only non_zero values are selected to calcluate MI
		
		# non_zero_indices = [j for j, x in enumerate(gene_snps) if x != 0]

		# X = list()
		# Y = list()

		# for j in range(len(non_zero_indices)):
		# 	X.append(gene_snps[non_zero_indices[j]])
		# 	Y.append(encoded_label[non_zero_indices[j]])	

		MI = calc_MI(gene_snps, encoded_label, bins)
		# setting a threshold for cutoff value of MI (0.001)
		if(gene not in gene_MI_dict):
			gene_MI_dict[gene] = MI

		if(MI > 0.001):	
			# print(clstr_genes[i], MI)
			important_gene_list.append(clstr_genes[i])

	return list(set(important_gene_list))

# to see the distribution of MI across gene
# def plot_gene_MI():
# 	MI_list = list()
# 	for key in gene_MI_dict:
# 		MI_list.append(gene_MI_dict[key])

# 	plt.hist(MI_list, bins = 'auto')
# 	plt.xlabel('MI', fontsize = 14)
# 	plt.ylabel('Number of genes', fontsize = 14)
# 	plt.title('Gene MI distribution')
# 	plt.xlim(0, 0.03)
# 	plt.grid()
# 	plt.show()

def filter_genes_from_cluster(clstr_genes, gene_snp_count_dict):

	# select top 20% of the genes based on the snp count
	clst_gene_dict = dict()
	for i in range(len(clstr_genes)):
		if(clstr_genes[i] in gene_snp_count_dict):
			clst_gene_dict[clstr_genes[i]] = gene_snp_count_dict[clstr_genes[i]]

	# sort the genes based on their snp count
	sorted_clstr_genes = sorted(clst_gene_dict.items(), key = operator.itemgetter(1), reverse = True)

	top_gene = int(len(clstr_genes) * 10 / 100) + 1		# 10 is for top 10% of the genes in the cluster

	selected_genes = list()
	for i in range(top_gene):
		selected_genes.append(sorted_clstr_genes[i][0])

	return selected_genes


def main():
	data_file = sys.argv[1]		# all_data.txt (original data file)
	clstr_file_list = glob.glob('*.clstr')	# list of all cluster files
	bins = 5					# fixed number of bins 5?

	original_D = pd.read_table(data_file, sep = '\t', header = 'infer')

	# read gene snp count file
	gene_snp_count_dict = dict()
	with open('gene_snp_frequency_down.txt') as f:
		for line in f:
			line = line.strip()
			cols = line.split()
			gene_snp_count_dict[cols[0]] = int(cols[1])

	# read each cluster file
	cluster_gene_mi_dict = dict()
	for i in range(len(clstr_file_list)):
		file = clstr_file_list[i]
		clst = file.split('.clstr')[0]
		print(clst)

		clstr_genes = pd.read_table(file, header = None)[0].tolist()
		selected_genes_from_cluster = filter_genes_from_cluster(clstr_genes, gene_snp_count_dict)
		# print(selected_genes_from_cluster)
		cluster_gene_mi_dict[clst] = mutual_info(selected_genes_from_cluster, original_D, bins)

	# plot_gene_MI()

	f_gene = open('selected_genes_down.txt', 'w')
	selected_genes = list()
	for key in cluster_gene_mi_dict:
		# print(key, cluster_gene_mi_dict[key])
		if(len(cluster_gene_mi_dict[key]) > 0):
			for i in range(len(cluster_gene_mi_dict[key])):
				selected_genes.append(cluster_gene_mi_dict[key][i])

	unique_genes = list(set(selected_genes))
	for i in range(len(unique_genes)):
		f_gene.write(str(unique_genes[i]) + '\n')
	f_gene.close()

if __name__ == "__main__":
	main()