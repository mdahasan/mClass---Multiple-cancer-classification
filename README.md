mClass - Multiple cancer type classifier with somatic point mutation data

mClass is a feature (gene) selection method for multiple cancer type classification using only somatic point mutation data. For this experiment, we have collected the dataset from http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/

First we process the dataset for each cancer type. The following code process each cancer mutation dataset and produce file with the caner type name and mutation count for each gene across all samples for that cancer.

```
python 1_data_preprocess.py <path for the cancer data directory>
``` 
Each caner data is then processed to eliminate samples that contains below threshold mutation count for all genes. This step produces files with reduced number of samples. All resampled data are then combined and `all_data.txt` file is generated.

```
python 2_sampling_data.py <path for the cancer data directory>
```

The next step cluster the genes based on cosine similarity and produces cluster files.

```
python 3_gene_clustering.py all_data.txt
```

The gene filtering step then filter out the genes based on a threshold (v) and generate `selected_genes.txt` file for further experiment.

```
python 4_gene_filtering.py all_data.txt
```

This step calculates normalized mutual information for all pairs of genes and produces a file named `All_Pair_wise_gene_MI.txt`

```
python 5_pair_wise_mi.py all_data.txt selected_genes.txt
```

Genes are then ranked using the selected gene list and their normalized mutual information values cacluated from the previous step. This step produces the final set of features (genes) that are ranked and fed for evaluating the classification accuracy in feed-forward feature selection way.

```
python 6_feature_selection_with_mi.py selected_genes.txt
```

The last step is evaluation of the mClass feature selection. We have used logistic regression in a OvR fashion for the multiple class classification. The code contains methods for cross validation and methods to find the testing accuracy. The feed forward feature selection stars with 50 genes and measures the accuracy by adding one feature (gene) at a time.

```
python 7_classification_result.py all_data.txt Final_feature_set.txt <number of genes to use>
```

