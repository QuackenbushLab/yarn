#' Check tissues to merge based on gene expression profile
#'
#' @param obj ExpressionSet object.
#' @param majorGroups Column name in the phenoData slot that describes the general body region or site of the sample.
#' @param minorGroups Column name in the phenoData slot that describes the specific body region or site of the sample.
#' @param filterFun Filter group specific genes that might disrupt PCoA analysis.
#' @param plotFlag TRUE/FALSE whether to plot or not
#' @param ... Parameters that can go to \code{\link[yarn]{checkMisAnnotation}}
#'
#' @return CMDS Plots of the majorGroupss colored by the minorGroupss. Optional matrix of CMDS loadings for each comparison.
#' @export
#'
#' @seealso checkTissuesToMerge
#'
#' @examples
#' data(skin)
#' checkTissuesToMerge(skin,'SMTS','SMTSD')
#'
checkTissuesToMerge <- function(obj, majorGroups, minorGroups,
                                filterFun = NULL, plotFlag = TRUE, ...) {
  if (length(majorGroups) == 1) {
    region <- factor(pData(obj)[, majorGroups])
  } else {
    region <- factor(majorGroups)
  }
  if (!is.null(filterFun)) {
    obj <- filterFun(obj)
  }
  result <- lapply(levels(region), function(i) {
    keepSamples <- which(region == i)
    objSubset <- obj[, keepSamples]
    objSubset <- objSubset[which(rowSums(exprs(objSubset)) > 0), ]
    res <- checkMisAnnotation(objSubset, phenotype = minorGroups,
                              controlGenes = "all", plotFlag = plotFlag, main = i,
                              ...)
    res
  })
  invisible(result)
}
