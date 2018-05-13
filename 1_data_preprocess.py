import numpy as np 
import pandas as pd
import sys
import os
import glob

def read_gene_names(path):
	# if the directory contains the file 'all_gene.txt' then take the gene names from there
	# if not, use the MANIFEST.txt file to get the name of all sample and iterate through 
	# them to extract the gene names
 	os.chdir(path)
	genes = list()
	if(glob.glob('all_gene.txt')):
		with open('all_gene.txt', 'r') as f:
			for line in f:
				line = line.strip()
				genes.append(line)
	else:
		file_name_list = list()
		with open('MANIFEST.txt', 'r') as mf:
			for line in mf:
				line = line.strip()
				cols = line.split()

				file_name_list.append(cols[1])

		# read each sample file to get the gene names
		for file in range(len(file_name_list)):
			current_gene_list = list()
			with open(file_name_list[file]) as f:
				title_line = f.readline()
				for line in f:
					line = line.strip()
					col = line.split()
					current_gene_list.append(col[0])

			genes.extend(set(current_gene_list) - set(genes))	# only including the genes that are not previously found

	return genes

def read_sample_files(path, gene_index_dict):
	os.chdir(path)
	file_name_list = list()
	with open('MANIFEST.txt', 'r') as mf:
		for line in mf:
			line = line.strip()
			cols = line.split()

			file_name_list.append(cols[1])

	# for first sample 
	# other samples will be vstack upon this on
	first_file_name = file_name_list[0]
	all_sample_cancer_snp_data = np.zeros(len(gene_index_dict), dtype = np.int)
	with open(first_file_name, 'r') as f:
		title_line = f.readline()
		title_cols = title_line.split('\t')

		hugo_symbol_index = title_cols.index('Hugo_Symbol')
		variant_type_index = title_cols.index('Variant_Type')

		for line in f:
			line = line.strip()
			cols = line.split()

			gene = cols[hugo_symbol_index]
			variant = cols[variant_type_index]

			if(variant == 'SNP'):
				gene_index = gene_index_dict[gene]
				all_sample_cancer_snp_data[gene_index] += 1

	# from the second sample file, read the file and create snp data
	# and then stack this array upon the first sample snp data (all_sample_cancer_snp_data)
	for file in range(1, len(file_name_list)):
		file_name = file_name_list[file]
		sample_snp_data = np.zeros(len(gene_index_dict), dtype = np.int)
		with open(file_name, 'r') as f:
			title_line = f.readline()
			title_cols = title_line.split('\t')

			hugo_symbol_index = title_cols.index('Hugo_Symbol')
			variant_type_index = title_cols.index('Variant_Type')

			for line in f:
				line = line.strip()
				cols = line.split()

				gene = cols[hugo_symbol_index]
				variant = cols[variant_type_index]

				if(variant == 'SNP'):
					gene_index = gene_index_dict[gene]
					sample_snp_data[gene_index] += 1

		all_sample_cancer_snp_data = np.vstack((all_sample_cancer_snp_data, sample_snp_data))
	
	return (all_sample_cancer_snp_data, file_name_list)

def write_cancer_snp_data_into_file(data, all_cancer_gene, sample_name_list, cancer_name):
	# write file for individual cancer dataset with snp frequency. also added 'Sample' and 'Cancer_type' columns
	fopen = open(cancer_name + '.txt', 'w')
	fopen.write('Sample' + '\t')
	for gene in range(len(all_cancer_gene)):
		fopen.write(all_cancer_gene[gene] + '\t')
	fopen.write('Cancer_type' + '\n')

	row, col = data.shape
	print(str(cancer_name) + ' ' + str(row) + ' ' + str(col))
	for i in xrange(row):
		fopen.write(str(sample_name_list[i].split('.maf')[0]) + '\t')
		for j in xrange(col):
			fopen.write(str(data[i, j]) + '\t')
		fopen.write(str(cancer_name) + '\n')

	fopen.close()

def main():
	data_path = sys.argv[1]			# path for the dataset directory
	current_path = os.getcwd()		# path for current working directory
	os.chdir(data_path)
	cancer_directoreis = [name for name in os.listdir(".") if os.path.isdir(name)]
	os.chdir(current_path)

	# collect the name of all genes from all cancer dataset
	# the file 'all_gene.txt' contains the name of genes associated with each cancer
	all_cancer_gene = list()
	for cancer in range(len(cancer_directoreis)):
		cancer_directory_path = os.path.join(data_path, cancer_directoreis[cancer])
		gene_list = read_gene_names(cancer_directory_path)
		os.chdir(current_path)

		all_cancer_gene.extend(set(gene_list) - set(all_cancer_gene))	# only includes non-repeating genes

	gene_index_dict = dict()
	gene_index = 0
	for gene in range(len(all_cancer_gene)):
		gene_index_dict[all_cancer_gene[gene]] = gene_index
		gene_index += 1

	# get the name of each sample from the MANIFEST.txt file
	# iterate through each file by name and update cancer matrix
	for cancer in range(0, len(cancer_directoreis)):
		cancer_name = cancer_directoreis[cancer]
		cancer_directory_path = os.path.join(data_path, cancer_directoreis[cancer])
		(individual_cancer_snp_data, sample_name_list) = read_sample_files(cancer_directory_path, gene_index_dict)
		# print(cancer_name, len(sample_name_list), np.sum(np.sum(individual_cancer_snp_data, axis = 1)))
		os.chdir(data_path)
		write_cancer_snp_data_into_file(np.array(individual_cancer_snp_data), all_cancer_gene, sample_name_list, cancer_name)
		os.chdir(current_path)		# changing directory to current one

if __name__ == "__main__":
	main()