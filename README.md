mClass - Multiple cancer type classifier with somatic point mutation data

mClass is a feature (gene) selection method for multiple cancer type classification using only somatic point mutation data. For this experiment, we have collected the dataset from http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/

First we process the dataset for each cancer type. The following code process each cancer mutation dataset and produce file with the caner type name and mutation count for each gene across all samples for that cancer.

```
python 1_data_preprocess.py <path for data directory>
``` 