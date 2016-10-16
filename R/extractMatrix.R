#' Extract the appropriate matrix
#'
#' This returns the raw counts, log2-transformed raw counts, or normalized expression.
#' If normalized = TRUE then the log paramater is ignored.
#'
#' @param obj ExpressionSet object or objrix.
#' @param normalized TRUE / FALSE, use the normalized matrix or raw counts
#' @param log TRUE/FALSE log2-transform.
#'
#' @importFrom Biobase assayData
#' @return matrix
#' @examples
#'
#' data(skin)
#' head(yarn:::extractMatrix(skin,normalized=FALSE,log=TRUE))
#' head(yarn:::extractMatrix(skin,normalized=FALSE,log=FALSE))
#'
extractMatrix <- function(obj, normalized = FALSE, log = TRUE) {
  if (class(obj) == "ExpressionSet") {
    if (!normalized) {
      obj <- exprs(obj)
    } else {
      if (!"normalizedMatrix" %in% names(assayData(obj)))
        stop("normalizedMatrix assayData missing")
      obj <- assayData(obj)[["normalizedMatrix"]]
      if (log & normalized)
        message("normalizedMatrix is assumed to already be log-transformed")
      log <- FALSE
    }
  }
  if (log == TRUE) {
    obj <- log2(obj + 1)
  }
  obj
}
