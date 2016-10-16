#' Normalize in a tissue aware context
#'
#' This function provides a wrapper to various normalization methods developed.
#' Currently it only wraps qsmooth and quantile normalization returning a log-transformed
#' normalized matrix. qsmooth is a normalization approach that normalizes samples in
#' a condition aware manner.
#'
#' @param obj ExpressionSet object
#' @param groups Vector of labels for each sample or a column name of the phenoData slot
#' for the ids to filter. Default is the column names
#' @param normalizationMethod Choice of 'qsmooth' or 'quantile'
#' @param ... Options for \code{\link{qsmooth}} function or \code{\link[limma]{normalizeQuantiles}}
#'
#' @return ExpressionSet object with an assayData called normalizedMatrix
#' @export
#'
#' @source The function qsmooth comes from the qsmooth packages
#' currently available on github under user 'kokrah'.
#'
#' @importFrom limma normalizeQuantiles
#' @importFrom Biobase storageMode
#' @importFrom Biobase storageMode<-
#' @importFrom Biobase assayData
#' @importFrom Biobase assayData<-
#' @importFrom preprocessCore normalize.quantiles
#' @importClassesFrom Biobase eSet
#' @importClassesFrom Biobase ExpressionSet
#'
#' @examples
#' data(skin)
#' normalizeTissueAware(skin,"SMTSD")
#'
normalizeTissueAware <- function(obj, groups, normalizationMethod = c("qsmooth",
                                                                      "quantile"), ...) {
  normalizationMethod <- match.arg(normalizationMethod)
  if (length(groups) == 1) {
    groups <- factor(pData(obj)[, groups])
  }
  storageMode(obj) <- "environment"
  if (normalizationMethod == "qsmooth") {
    normalizedMatrix <- qsmooth(obj, groups = groups,
                                ...)
  } else if (normalizationMethod == "quantile") {
    if(length(unique(groups))>1){
      normalizedMatrix <- sapply(unique(groups), function(i) {
        cnts <- exprs(obj[, which(pData(obj)$our %in% i)])
        nmat <- normalize.quantiles(cnts)
        colnames(nmat) <- colnames(cnts)
        nmat
      })
      normalizedMatrix <- Reduce("cbind", normalizedMatrix)
      normalizedMatrix <- normalizedMatrix[, match(colnames(obj),
                                                   colnames(normalizedMatrix))]
    } else {
      normalizedMatrix <- normalize.quantiles(exprs(obj))
      colnames(normalizedMatrix) <- colnames(obj)
    }
  }
  assayData(obj)[["normalizedMatrix"]] <- normalizedMatrix
  storageMode(obj) <- "lockedEnvironment"
  obj
}
