#' Filter specific genes
#'
#' The main use case for this function is the removal of sex-chromosome genes.
#' Alternatively, filter genes that are not protein-coding.
#'
#' @param obj ExpressionSet object.
#' @param labels Labels of genes to filter or keep, eg. X, Y, and MT
#' @param featureName FeatureData column name, eg. chr
#' @param keepOnly Filter or keep only the genes with those labels
#'
#' @return Filtered ExpressionSet object
#' @export
#'
#' @importFrom Biobase exprs
#' @importFrom Biobase fData
#'
#' @examples
#' data(skin)
#' filterGenes(skin,labels = c('X','Y','MT'),featureName='chromosome_name')
#' filterGenes(skin,labels = 'protein_coding',featureName='gene_biotype',keepOnly=TRUE)
#'
filterGenes <- function(obj, labels = c("X", "Y", "MT"), featureName = "chromosome_name",
                        keepOnly = FALSE) {
  features <- fData(obj)[, featureName]
  if (keepOnly == FALSE) {
    throwAwayGenes <- which(features %in% labels)
  } else {
    throwAwayGenes <- which(!features %in% labels)
  }
  obj <- obj[-throwAwayGenes, ]
  obj
}
