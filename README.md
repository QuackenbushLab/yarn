[![Travis-CI Build Status](https://travis-ci.org/jnpaulson/yarn.svg?branch=master)](https://travis-ci.org/QuackenbushLab/yarn)

# YARN: Robust Multi-Tissue RNA-Seq Preprocessing and Normalization

The goal of yarn is to expedite large RNA-seq analyses using a combination of previously developed tools. Yarn is meant to make it easier for the user to perform accurate comparison of conditions by leveraging many Bioconductor tools and various statistical and normalization techniques while accounting for the large heterogeneity and sparsity found in very large RNA-seq experiments. 

## Installation

You can install yarn from github with:

```R
# install.packages("devtools")
devtools::install_github("quackenbushlab/yarn")
```

## Example


This is a basic workflow in terms of code: 


0. First always remember to have the library loaded.
```R
library(yarn)
```

1.  Download the GTEx gene count data as an ExpressionSet object or load the sample skin dataset.
```R
library(yarn)
data(skin)
```

2. Check mis-annotation of gender or other phenotypes using group-specific genes
```R
checkMisAnnotation(skin,"GENDER",controlGenes="Y",legendPosition="topleft")
```

3. Decide what sub-groups should be merged
```R
checkTissuesToMerge(skin,"SMTS","SMTSD")
```

4. Filter lowly expressed genes
```R
skin_filtered = filterLowGenes(skin,"SMTSD")
dim(skin)
dim(skin_filtered)
# Or group specific genes
tmp = filterGenes(skin,labels=c("X","Y","MT"),featureName = "chromosome_name")
# Keep only the sex names
tmp = filterGenes(skin,labels=c("X","Y","MT"),featureName = "chromosome_name",keepOnly=TRUE)
```

5. Normalize in a tissue or group-aware manner
```R
plotDensity(skin_filtered,"SMTSD",main="log2 raw counts")
skin_filtered = normalizeTissueAware(skin_filtered,"SMTSD")
plotDensity(skin_filtered,"SMTSD",normalized=TRUE,main="Normalized")
```
