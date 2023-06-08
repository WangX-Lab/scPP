# ScPP

A simple and effective algorithm for recognizing cell subpopulations with specific phenotypes based on the expression profiles of phenotype-associated marker genes in bulks and single cells.

# Installation

Users can install the released version of ScPP with:

```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("WangX-Lab/ScPP")
```

# Description

To infer phenotypes of single cells from scRNA-seq data, ScPP requires three types of data as input, including single-cell transcriptomes, bulk transcriptomes, and phenotypic features in bulk data.The phenotypic features can be categorical variables, continuous variables, or clinical survival.

# Details

+ The function `singcell_Preprocessing()` is for single cell expression data preprocessing. Its input is a count matrix with cell names in columns and gene names in rows.

+ The function `marker_Binary()` is used for binary variables to generating marker genes or signatures correlated with two groups, containing 4 parameters: bulk_data, features, ref_group and Log2FC_cutoff.
  + "bulk_data" is a Log2-normalized bulk expression data with genes in row and samples in column.
  + "features" is the feature data of bulk samples, column1 are sample names (colname is "Sample") and column2 are feature (colname is "Feature") labels of each sample.
  + "ref_group" is a character to indicate which feature is the control group.
  + "Log2FC_cutoff" is the absolute cutoff value of fold change, with default 0.585.

+ The function `marker_Continuous()` is used for continuous variables to generating marker genes or signatures correlated with this feature, containing 4 parameters: bulk_data, features, method and estimate_cutoff.
  + "bulk_data" is a Log2-normalized bulk expression data with genes in row and samples in column.
  + "features" is the feature data of bulk samples, such as TMB or CNA values of each sample.
  + "method" is the method uses for cor.test(), with default "spearman", and another choice is "pearson".
  + "estimate_cutoff" is the absolute cutoff value of correlation coefficient, with default 0.2.

+ The function `marker_Survival()` is used for survival data to generating marker genes or signatures correlated with patients'prognosis, containing 2 parameters: bulk_data and survival_data.
  + "bulk_data" is a Log2-normalized bulk expression data with genes in row and samples in column.
  + "survival_data" is the survival data with time in column1 and status in column2. Row names of survival_data are sample names.

+ The function `ScPP()` is used for  for Single Cells’Phenotype Prediction,  containing 3 parameters: sc_dataset, geneList and probs.
  + "sc_dataset" is a seurat object of single cell RNA sequencing data, it can the output of function `singcell_Preprocessing()`.
  + "geneList" is a gene list correlated with interested features, it can be the output of functions `marker_Binary()`, `marker_Continuous()` and `marker_Survival()` for binary variables, continuous varaibles and survival data, respectively.
  + "probs" is the α value of ScPP, with default 0.2.

# Examples

## **Binary variables**

```R
library(ScPP)
load(system.file("data/binary.RData",package = "ScPP"))
sc = sc_Preprocess(sc_count)
geneList = marker_Binary(bulk, binary, ref_group = "Normal")
metadata = ScPP(sc, geneList)
head(metadata)
sc$ScPP = metadata$ScPP
Idents(sc) = "ScPP"

#Visualization of ScPP-identified cells
DimPlot(sc, group = "ScPP", cols = c("grey","red","blue"))
```

<img width="637" alt="image" src="https://github.com/WangX-Lab/ScPP/assets/54932820/b8784ad9-c911-47ff-b921-667de946ecdc">



## **Continuous variables**

```R
library(ScPP)
load(system.file("data/continuous.RData",package = "ScPP"))
sc = sc_Preprocess(sc_count)
geneList = marker_Continuous(bulk, continuous$TMB_non_silent)
metadata = ScPP(sc, geneList)
sc$ScPP = metadata$ScPP
Idents(sc) = "ScPP"

#Visualization of ScPP-identified cells
DimPlot(sc, group = "ScPP", cols = c("grey","red","blue"))
```

<img width="637" alt="image" src="https://github.com/WangX-Lab/ScPP/assets/54932820/2796dd06-a015-4852-9574-2a7b6d65772c">


## **Survival data**

```R
library(ScPP)
load(system.file("data/survival.RData",package = "ScPP"))
sc = sc_Preprocess(sc_count)
geneList = marker_Survival(bulk, survival)
metadata = ScPP(sc, geneList)
sc$ScPP = metadata$ScPP
Idents(sc) = "ScPP"

#Visualization of ScPP-identified cells
DimPlot(sc, group = "ScPP", cols = c("grey","red","blue"))
```

<img width="637" alt="image" src="https://github.com/WangX-Lab/ScPP/assets/54932820/d3ca4d5f-fbe4-4867-926c-84c84482917a">

# Contact

E-mail any questions to Xiaosheng Wang (xiaosheng.wang@cpu.edu.cn)
