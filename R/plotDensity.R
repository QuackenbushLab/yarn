#' Density plots of columns in a matrix
#'
#' Plots the density of the columns of a matrix. Wrapper for \code{\link[quantro]{matdensity}}.
#'
#' @param obj ExpressionSet object
#' @param groups Vector of labels for each sample or a column name of the phenoData slot
#' for the ids to filter. Default is the column names.
#' @param normalized TRUE / FALSE, use the normalized matrix or log2-transformed raw counts
#' @param legendPos Legend title position. If null, does not create legend by default.
#' @param ... Extra parameters for \link[quantro]{matdensity}.
#'
#' @return A density plot for each column in the ExpressionSet object colored by groups
#' @export
#'
#' @importFrom quantro matdensity
#' @importFrom Biobase assayData
#' @importFrom Biobase storageMode
#' @importFrom graphics legend
#'
#' @examples
#' data(skin)
#' filtData <- filterLowGenes(skin,"SMTSD")
#' plotDensity(filtData,groups="SMTSD",legendPos="topleft")
#' # to remove the legend
#' plotDensity(filtData,groups="SMTSD")
#'
plotDensity <- function(obj, groups = NULL, normalized = FALSE,
                        legendPos = NULL, ...) {
  if (length(groups) == 1) {
    groups <- factor(pData(obj)[, groups])
  }
  mat <- extractMatrix(obj, normalized, log = TRUE)
  matdensity(mat, groupFactor = groups, ...)
  if (!is.null(legendPos)) {
    legend(legendPos, legend = levels(groups), fill = 1:length(levels(groups)),
           box.col = NA)
  }
}
