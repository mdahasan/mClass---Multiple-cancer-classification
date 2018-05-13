import pandas as pd 
import numpy as np 
import sys
import operator

def pair_wise_mi_matrix_to_dict(X, selected_gene_names, feature_size):

	# converting pair_wise_mi matrix to a dict(), redundant 
	Gene_pair_MI_dict = dict()

	for i in range(0, feature_size):
		gene1_name = selected_gene_names[i]
		for j in range(i + 1, feature_size):
			gene2_name = selected_gene_names[j]

			pair = (gene1_name, gene2_name)

			Gene_pair_MI_dict[pair] = X.iloc[i][j]

	f_pair_mi = open('Pair_wise_mi_dict.txt', 'w')
	for key in Gene_pair_MI_dict:
		pair = key
		f_pair_mi.write(str(pair[0]) + '\t' + str(pair[1]) + '\t' + str(Gene_pair_MI_dict[key]) + '\n')

	return Gene_pair_MI_dict

def greedy_feature_selection(gene_names, class_feat_mi_dict, Gene_pair_MI_dict, Gene_entropy_dict, S):

	# formula: G = I(c, f) = (1/|s|)*sum(I(f_i, f_j)/min(H(f_i), H(f_j)))
	max_G = -9999

	for i in range(len(gene_names)):
		# first check if the gene has already been selected or not
		if(gene_names[i] in S):
			continue
		else:
			current_gene = gene_names[i]

			if(current_gene in class_feat_mi_dict and current_gene in Gene_entropy_dict):
				I_cf = class_feat_mi_dict[current_gene]

				I_fi_fs_tot = 0
				for j in range(len(S)):
					selected_gene = S[j]
					pair = (current_gene, selected_gene)
					current_gene_entropy = Gene_entropy_dict[current_gene]
					selected_gene_entropy = Gene_entropy_dict[selected_gene]

					I_fi_fs_tot += Gene_pair_MI_dict[pair] / float(min(current_gene_entropy, selected_gene_entropy))

				G = I_cf - (I_fi_fs_tot / float(len(S)))

				if(G > max_G):
					max_G = G
					gene_index = i

	# print(gene_names[gene_index], G, max_G)

	return gene_names[gene_index]

def write_selected_genes_to_file(S, t = 0.55):

	# writing the selected genes to file for classification
	# f_sel_gene = open('Final_feature_set_' + str(t) + '.txt', 'w')
	f_sel_gene = open('Final_feature_set.txt', 'w')
	for i in range(len(S)):
		f_sel_gene.write(str(S[i]) + '\n')
	f_sel_gene.close()

def main():
	# pair_wise_mi_file = sys.argv[1]		# pair_wise_mi_dict.txt file
	selected_genes_file = sys.argv[1]	# selected_genes.txt file
	# sim_threshold = sys.argv[2]			# this is temporary
	# class_feat_mi_file = sys.argv[3]	# class_feat_mi.txt file

	# list for selected feature
	S = list()

	# first read the pair_wise_MI file and store the values in a dictionary
	# X = pd.read_table(pair_wise_mi_file, sep = '\t', header = None)
	# sample_size, feature_size = X.shape

	selected_gene_names = pd.read_table(selected_genes_file, header = None)[0].tolist()
	
	# read class_feat_mi into a dict
	class_feat_mi_dict = dict()
	with open('All_Class_feature_MI_down.txt') as f:
		for line in f:
			line = line.strip()
			cols = line.split()

			class_feat_mi_dict[cols[0]] = float(cols[1])

	sorted_class_feat_mi = sorted(class_feat_mi_dict.items(), key = operator.itemgetter(1), reverse = True)
	sorted_class_feat_mi_dict = dict()
	for i in range(len(sorted_class_feat_mi)):
		sorted_class_feat_mi_dict[sorted_class_feat_mi[0]] = sorted_class_feat_mi[1]

	# write the mair wise MI values to a file
	# Gene_pair_MI_dict = pair_wise_mi_matrix_to_dict(X, selected_gene_names, feature_size - 1) # number of feature is one less than the dimension of columns (genes)

	# reading the MI value for each gene pair from a file
	# this can be done in program wihthout reading from file

	Gene_pair_MI_dict = dict()
	with open('All_Gene_pair_MI_down.txt', 'r') as f:
		for line in f:
			line = line.strip()
			cols = line.split()

			gene_pair = (cols[0], cols[1])
			gene_rev_pair = (cols[1], cols[0])
			Gene_pair_MI_dict[gene_pair] = float(cols[2])
			Gene_pair_MI_dict[gene_rev_pair] = float(cols[2])

	# reading the gene entropy file for H(gene) 
	Gene_entropy_dict = dict()
	with open('All_Gene_entropy_values_down.txt') as f:
		for line in f:
			line = line.strip()
			cols = line.split()

			Gene_entropy_dict[cols[0]] = float(cols[1])



	features = list(sorted_class_feat_mi_dict.keys())
	first_feature = features[0][0]

	# number of desired feature
	k = len(selected_gene_names)

	print('Total number of genes ' + str(len(selected_gene_names)))

	S.append(first_feature)						# adding the feature that has been selected as important
	# selected_gene_names.remove(first_feature)	# removing the feature that has been selected from overall gene list

	while(len(S) < k):
		next_feature = greedy_feature_selection(selected_gene_names, class_feat_mi_dict, Gene_pair_MI_dict, Gene_entropy_dict, S)
		S.append(next_feature)
		
		if(len(S) % 100 == 0):
			print(str(len(S)) + ' genes processed...')

	write_selected_genes_to_file(S, sim_threshold)

if __name__ == "__main__":
	main()