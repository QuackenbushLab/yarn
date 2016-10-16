#' Filter genes not expressed in any sample
#'
#' The main use case for this function is the removal of missing genes.
#'
#' @param obj ExpressionSet object.
#' @param threshold Minimum sum of gene counts across samples -- defaults to zero.
#'
#' @return Filtered ExpressionSet object
#' @export
#'
#' @importFrom Biobase exprs
#' @importFrom Biobase fData
#'
#' @examples
#' data(skin)
#' filterMissingGenes(skin)
#'
filterMissingGenes <- function(obj, threshold = 0) {
  sumGenes <- rowSums(exprs(obj))
  throwAwayGenes <- which(sumGenes <= threshold)
  if (length(which(sumGenes <= 0)) > 0) {
    obj <- obj[-throwAwayGenes, ]
  }
  obj
}
