import pandas as pd
import numpy as np 
import sys
import math
import random
import string
# import matplotlib.pyplot as plt 
import itertools

from sklearn import preprocessing as pp  
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn import svm
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.metrics import confusion_matrix, precision_score, recall_score, f1_score, r2_score, classification_report, accuracy_score

# def plot_confusion_matrix(cm, classes,
#                       normalize=False,
#                       title='Confusion matrix',
#                       cmap=plt.cm.Blues):
#     """
#     This function prints and plots the confusion matrix.
#     Normalization can be applied by setting `normalize=True`.
#     """
#     if normalize:
#         cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
#         print("Normalized confusion matrix")
#     else:
#         print('Confusion matrix, without normalization')

#     print(cm)

#     plt.imshow(cm, interpolation='nearest', cmap=cmap)
#     plt.title(title)
#     plt.colorbar()
#     tick_marks = np.arange(len(classes))
#     plt.xticks(tick_marks, classes, rotation=45)
#     plt.yticks(tick_marks, classes)

#     fmt = '.2f' if normalize else 'd'
#     thresh = cm.max() / 2.
#     for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
#         plt.text(j, i, format(cm[i, j], fmt),
#                  horizontalalignment="center",
#                  color="white" if cm[i, j] > thresh else "black")

#     plt.tight_layout()
#     plt.ylabel('True label', fontsize = 14)
#     plt.xlabel('Predicted label', fontsize = 14)

# def classify_data(X, y, class_names):

# 	X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.1, random_state = 42)

# 	# clf = RandomForestClassifier(n_estimators = 1000, class_weight = "balanced")
# 	clf = LogisticRegression(penalty = 'l2', class_weight = 'balanced')
# 	# clf = svm.SVC(kernel='linear', C = 2e3)
# 	y_pred = clf.fit(X_train, y_train).predict(X_test)

# 	print('Precision Score:')
# 	print(precision_score(y_test, y_pred, average = 'macro'))

# 	print('Recall Score:')
# 	print(recall_score(y_test, y_pred, average = 'micro'))

# 	print('F1 Score: ')
# 	print(f1_score(y_test, y_pred, average = 'weighted'))

# 	print('R2 Score:')
# 	print(r2_score(y_test, y_pred, multioutput = 'variance_weighted'))

# 	print('Classification report:')
# 	print(classification_report(y_test, y_pred))

# 	print('Accuracy score: ')
# 	print(accuracy_score(y_test, y_pred))
# 	# return accuracy_score(y_test, y_pred)

# 	# Compute confusion matrix
# 	cnf_matrix = confusion_matrix(y_test, y_pred)
# 	np.set_printoptions(precision=2)

# 	# Plot non-normalized confusion matrix
# 	plt.figure()
# 	plot_confusion_matrix(cnf_matrix, classes=class_names,
# 	                      title='Confusion matrix, without normalization')
# 	plt.tight_layout()
# 	# plt.savefig('Confusion_matrix_not_normalized.png', dpi = 120)
# 	plt.show()
	
# 	# Plot normalized confusion matrix
# 	plt.figure()
# 	plot_confusion_matrix(cnf_matrix, classes=class_names, normalize=True,
# 	                      title='Normalized confusion matrix')
# 	plt.tight_layout()
# 	# plt.savefig('Confusion_matrix_normalized.png', dpi = 120)
# 	plt.show()

def cross_validation(X, new_y):
	# clf = svm.SVC(kernel = 'rbf', C = 2e4, gamma = 2e-5)
	clf = LogisticRegression(penalty = 'l2', class_weight = 'balanced')
	# clf = KNeighborsClassifier(p = 2)
	# clf = GaussianNB()
	scores = cross_val_score(clf, X, new_y, cv = 10)
	# print(scores)
	return scores


def main():
	data_file = sys.argv[1]		# data file containing all snps info (all_data.txt)
	filtered_gene_file = sys.argv[2]	# filtered gene file (Final_feature.txt)
	number_of_feature = sys.argv[3]		# number of features to be selected for classification
	# dist_threshold = sys.argv[4]		# temporary (55, 60, 65, ..., 75)

	# file to store ffs results
	# f_cv_ffs = open('FFS_cv_acc_' + str(dist_threshold) + '.txt', 'w')
	f_cv_ffs = open('FFS_cv_acc.txt', 'w')

	original_X = pd.read_table(data_file, sep = '\t', header = 'infer')
	original_y = original_X['Cancer_type'].values
	original_y_label = list(set(original_y))
	# integer_label = [12, 5, 23, 19, 10, 2, 34, 21, 8, 54, 7, 13]
	integer_label = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
	original_X.drop('Cancer_type', axis = 1, inplace = True)

	new_y = list()
	cancer_name_label_dict = dict()
	for i in range(len(original_y)):
		label_index = original_y_label.index(original_y[i])
		new_y.append(integer_label[label_index])
		cancer_name_label_dict[integer_label[label_index]] = original_y[i]

	# print(cancer_name_label_dict)
	
	gene_list = pd.read_table(filtered_gene_file, header = None)[0].tolist()

	selected_gene_list = gene_list[0: int(number_of_feature)]

	X = original_X.filter(selected_gene_list, axis = 1)

	cancer_name = list()
	for i in range(len(integer_label)):
		cancer_name.append(cancer_name_label_dict[integer_label[i]])

	# classify_data(X, new_y, cancer_name)

	if(int(number_of_feature) > len(gene_list)):
		number_of_feature = len(gene_list)

	# check for every added 50 features
	for num_feat in range(50, number_of_feature):

		selected_gene_list = gene_list[0: int(num_feat)]
		
		X = original_X.filter(selected_gene_list, axis = 1)
		# print(X.shape)

		# acc = classify_data(X, new_y, integer_label)
		# f_ffs.write(str(num_feat) + '\t' + str(acc) + '\n')
		# print(num_feat, acc)
		cross_acc = cross_validation(X, new_y)
		print(num_feat, np.mean(cross_acc), np.max(cross_acc))
		# print(num_feat, acc)
		f_cv_ffs.write(str(num_feat) + '\t' + str(np.mean(cross_acc)) + '\t' + str(np.max(cross_acc)) + '\n')

if __name__ == "__main__":
	main()