# yarn - Yet Another RNa-seq package

The goal of yarn is to expedite large RNA-seq analyses using a combination of previously developed tools. Yarn is meant to make it easier for the user to perform accurate comparison of conditions by leveraging many Bioconductor tools and various statistical and normalization techniques while accounting for the large heterogeneity and sparsity found in very large RNA-seq experiments. 

## Installation

You can install yarn from github with:

```R
# install.packages("devtools")
devtools::install_github("jnpaulson/yarn")
```

## Example

This is a basic workflow: 

1.  Download the GTEx gene count data as an ExpressionSet object or load the sample skin dataset.
```R
data(skin)
obj = downloadGTEx()
```
2. Check mis-annotation of gender
```R
checkMisAnnotation(obj,"GENDER",controlGenes="Y",legendPosition="topleft")
```
3. Decide what sub-groups should be merged
```R
checkTissuesToMerge(skin,"SMTS","SMTSD")
```
4. Filter lowly expressed genes or condition specific genes
```R
obj = filterLowGenes(skin,"SMTS")
obj = filterGenes(skin,labels=c("X","Y","MT"),featureName = "chromosome_name")
# Keep only the sex names
obj = filterGenes(skin,labels=c("X","Y","MT"),featureName = "chromosome_name",keepOnly=TRUE)
```
5. Normalize in a group-aware manner
```R
skin = normalizeTissueAware(skin,"SMTSD")
```
