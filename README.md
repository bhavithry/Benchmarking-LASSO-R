# Benchmarking penalized regression algorithmsin machine learning for single cell RNAsequencing dataset

# Table of content
* [Abstract](#abstract)
* [Technologies](#technologies)
* [Data](#data)
* [Step by step implementation](#step-by-step-implementation)
* [References](#references)


# Abstract

With the advent of Single-Cell RNA Sequencing (scRNA-seq) technology, gene expression can now be examined at a single-cell resolution. Processing scRNA-seq data has its own challenges, such as the high dimensionality of the data. Machine learning can be applied to do gene (feature) selection from high-dimensional scRNA-seq data. Even though there are numerous machine learning methods such as penalized regression for feature selection, there is no methodical comparison of their performances across different data sets. Therefore, we analyze penalized regression methods that are suitable for scRNA-seq data. With the 4 scRNA-seq data sets that we analyzed, results indicate that Sparse Group Lasso (SGL) outperforms the other 6 methods in terms of area under the receiver operating curve (AUC) and average computation time. Based on these findings, we propose a new method that is an ensemble of penalized regression methods. The proposed method works by creating a small pre-selected subset of genes and then applying SGL to it to select the differentially expressed genes in scRNA-seq data. By using hierarchical clustering for grouping genes, the proposed method eliminates the need for prior knowledge of gene grouping information. Furthermore, the proposed method has consistently improved AUC across data sets.

Key Words: Single Cell RNA Sequencing; Machine Learning; LASSO; Feature Selection; High Dimensional Data; R Programming Language.

# Technologies
Software: R Version 4.1.2

Operating System: GNU/Linux 5.4.0-96-generic x86_64

Cloud Server: Compute Canada

# Data
GSE60749 data set is downloaded from Conquer website (http://imlspenticton.uzh.ch:3838/conquer/) is of
species Mus musculus. This dataset was generated in the study of gene expression
variability in pluripotent stem cells (PSCs) by single-cell expression profiling of
PSCs under different chemical and genetic perturbations. Gene
expression levels are quantified as transcripts per million reads (TPM) in this
dataset. For our research we selected 183 individual v6.5 mouse embryonic stem
cells (mESCs) and 84 Dgcr8 -/- mESCs that lack mature miRNAs (knockout of
a miRNA processing factor) from this dataset. The 183 individual mESCs are
assigned to class 1 and 84 Dgcr8 -/- mESCs are assigned to class 0. The dataset
included 22443 genes initially which reduced to 15508 after data preprocessing
in which all the genes with no variance in expression across all the cells were
removed. Some of the genes in this dataset are noncoding piRNAs with numbers
as names. Such numbers are converted to text by prefixing them with ’RNA ’,
before loading data to R to make sure that the names are not converted as dates
in R.

# Step by step implementation
The new pipeline include below steps,

1. Load dataset to R and assign classes 1 and 0 to the two selected group of cells
to form a binary classification problem
2. Shuffle cells within each class
3. Remove genes with no variability in expression across all cells from the dataset
4. Split the dataset into training (90%) and testing data (10%) for 10-fold cross
validation
Repeat steps 5 to 9 for 10 folds

        5. Generate ridge, lasso, elastic net and drop lasso models with the training
        data
        6. Extract coefficients for each model, sort the coefficients and select
        those which are above the mean of coefficients for each model
        7. Form a gene pool by taking union of only the top genes from all 4
        models. Top genes are genes for which magnitude of coefficient is greater
        than mean of coefficients of the model. For example, Fig. 1 represents
        gene pool of data set GSE123818.
        8. Run SGL with new gene pool grouped by hierarchical clustering.
        9. Save the coefficients of SGL model
        
10. Average and sort the coefficients, and then visualize the gene Vs coefficients
plot
11. Select final number of genes by finding the elbow in plot
12. Visualize K-means clustering of the data set with only the final selection of
genes


![Venn Diagram](./bhavithry/Benchmarking-LASSO-R/venn_123818.png?raw=true "Fig 1")

# Acknowledgement
The authors would like to acknowledge the funding for this research from,
1. TRU Internal Research Fund (IRF) awarded to Dr. Jabed Tomal, Department
of Mathematics and Statistics, Thompson Rivers University, and
Dr Yan Yan, Department of Mathematics and Statistics, Thompson Rivers
University.
2. Natural Sciences and Engineering Research Council of Canada (NSERC)
awarded to Dr. Jabed Tomal, Department of Mathematics and Statistics, Thompson
Rivers University and Dr. Yan Yan, Department of Computing Science,
Thompson Rivers University.

The authors also acknowledges Compute Canada for hosting the 32GB Linux
remote server which is used for computation in this research.

# References
Puliparambil, B. S., Tomal, J., & Yan, Y. (2022). Benchmarking Penalized Regression Methods in Machine Learning for Single Cell RNA Sequencing Data. In RECOMB International Workshop on Comparative Genomics (pp. 295-310). Springer, Cham.
