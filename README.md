# Watershed


Watershed is an unsupervised probabilistic framework that integrates genomic annotations and RNA-seq outlier calls into one model to identify the probability a rare variant has a functional effect on a particular RNA-seq outlier phenotype (examples of phenotypes can be, but are not limited to total expression, splicing, or ASE). Watershed extends our previous model RIVER (which can also be run from this package) by incoorperating information from multiple outlier phenoytpes into one model, where predictions for functional effects in one outlier phenotype are informed by observed outlier calls in another phenotype. Please see [ref bioRxiv preprint] for more details.


## Input file
Watershed models instances of gene-individual pairs. Therefor each line in the Watershed input file must be a gene-individual pair. Each line must include a minimum of one genomic annotations describing the rare variants nearby the gene for that individual. Each line must also include a minimum of one outlier score (p-value) describing the gene-individual pair. To evaluate Watershed performance we rely on N2 pairs, or pairs of individuals that share the same rare genetic variant nearby a particular gene (see bioRxiv preprint for more details). Therefor, each line in the input file must contain information concerning whether the gene-individual is an N2 pair. If the gene-individual is not an N2 pair, it can be represented as "NA". If the gene-individual is in an N2 pair, it can be represented as a positive integer that represents the unique identifier for that N2 pair. Note: if there are more than two individauls with the same rare genetic variant, we randomly select two individuals and ignore all other individuals. 

More formally, assuming G genomic annotations and E outlier p-values, the columns of the Watershed input file are as follows (column ordering is essential):
1. "SubjectID": Identifier for the individual in the gene-individual pair
2. "GeneName": Identifier for the gene in the gene-individual pair
3. A column for each of the G genomic annotations. Header names of these columns can be entirely user specified. Values are real valued.
4. A column for each of the E outlier p-values. Header names of these columns can be entirely user specified. Values are p-values (ie [0,1]). 
5. "N2pair": Identifier for whether gene-individual is in an N2 pair. If it is not, use "NA". If it is an N2 pair, use a unique positive integer to identify the pair.

An example input file with 18 genomic annotations and 3 outlier p-values can be found in "example_data/watershed_example_data.txt".
Another example input file with 18 genomic annotations and 1 outlier p-value can be found in "example_data/river_example_data_pheno_1.txt"

Note on outlier p-values: We can explicitely model outlier calls with sign (ie. total expression outliers are generated from a zscore and can therefor be over-expression or under-expression). To model this, if a particular outlier signal has sign/direction, simply put a negative sign in front of the p-value for instances that are under-expression. An exmample of this can be seen in the "pheno_2_pvalue" column of "example_data/watershed_example_data.txt".


## Running Watershed
There are two main scripts useful to user looking to apply Watershed to their data:
1. "evaluate_watershed.R": Briefly, this script will train a Watershed model on non-N2 pairs, and evaluate (via precision recall curves) on held-out N2 pairs. This will allow the user to get an idea of the accuracy of Watershed applied to their data.
2. "predict_watershed.R": Briefly, this script allows a user to train a Watershed model on training data, and then predict Watershed posterior probabilities (using Watershed parameters optimized in training) on all gene-individual in a much larger prediction data set.


## Running evaluate_watershed.R
hello











## Authors

* **Ben Strober** -- [BennyStrobes](https://github.com/BennyStrobes) -- bstrober3@gmail.com