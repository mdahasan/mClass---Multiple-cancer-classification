import pandas as pd 
import numpy as np 
import sys

from sklearn import preprocessing as pp

def shannon_entropy(c):

	c_normalized = c / float(np.sum(c))
	c_normalized = c_normalized[np.nonzero(c_normalized)]
	H = -sum(c_normalized * np.log2(c_normalized))
	return H

def calc_mi(x, y, bins = 5):

	c_xy = np.histogram2d(x, y, bins)[0]
	c_x = np.histogram(x, bins)[0]
	c_y = np.histogram(y, bins)[0]

	h_x = shannon_entropy(c_x)
	h_y = shannon_entropy(c_y)
	h_xy = shannon_entropy(c_xy)

	mi = h_x + h_y - h_xy

	return mi, h_x, h_y

def class_feature_MI(X, y):
	sample_size, feature_size = X.shape

	# check each gene to find the max I(f_i, y)
	max_mi = 0

	# writing in a file (class - feature MI)
	f_class_feat = open('All_Class_feature_MI.txt', 'w')

	for i in range(feature_size):
		feat_i = X.iloc[:, i].tolist()

		mi, _, _ = calc_mi(feat_i, y)

		f_class_feat.write(str(X.columns[i]) + '\t' + str(mi) + '\n')

	f_class_feat.close()

def main():

	data_file = sys.argv[1]		# original snp count file (all_data.txt)
	selected_genes_file = sys.argv[2]	# filtered genes file (selected_gene_file.txt)

	original_X = pd.read_table(data_file, sep = '\t', header = 'infer')
	original_y = original_X['Cancer_type']

	original_y = original_X['Cancer_type'].values
	le = pp.LabelEncoder()
	le.fit(original_y)
	y = le.transform(original_y)

	original_X.drop('Cancer_type', axis = 1, inplace = True)

	selected_genes = pd.read_table(selected_genes_file, header = None)[0].tolist()

	X = original_X.filter(selected_genes, axis = 1)

	gene_names = X.columns.values

	sample_size, feature_size = X.shape
	print('Filtered Data Dim: ' + str(X.shape))

	# compute the entropy of all genes
	H_dict = dict()
	for i in range(0, feature_size):
		gene = X.iloc[:, i].tolist()
		_, h_x, _ = calc_mi(gene, gene)
		H_dict[gene_names[i]] = h_x

	# write the entropy values for each gene in a file
	f_h = open('All_Gene_entropy_values.txt', 'w')
	for key in H_dict:
		f_h.write(str(key) + '\t' + str(H_dict[key]) + '\n')
	f_h.close()

	# first compute all gene pair mutual information
	f_pair_mi = open('All_Gene_pair_MI.txt', 'w')
	MI_mat = np.zeros((feature_size, feature_size))

	for i in range(0, feature_size):
		gene1 = X.iloc[:, i].tolist()
		gene1_name = gene_names[i]

		for j in range(i + 1, feature_size):
			gene2 = X.iloc[:, j].tolist()
			gene2_name = gene_names[j]

			mi, h_g1, h_g2 = calc_mi(gene1, gene2)
			f_pair_mi.write(str(gene1_name) + '\t' + str(gene2_name) + '\t' + str(mi) + '\n')

			MI_mat[i, j] = mi

		if(i % 100 == 0 and i > 1):
			print(str(i) + ' genes processed...')

	f_pair_mi.close()

	# write the pair-wise MI values in a file
	f_mi = open('All_Pair_wise_gene_MI.txt', 'w')
	for i in range(feature_size):
		for j in range(feature_size):
			f_mi.write(str(MI_mat[i, j]) + '\t')
		f_mi.write('\n')
	f_mi.close()

	
	# get the MI for each feature with class f(f, c)
	class_feature_MI(X, y)

if __name__ == "__main__":
	main()