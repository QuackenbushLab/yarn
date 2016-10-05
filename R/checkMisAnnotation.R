#' Check for wrong annotation of a sample using classical MDS and control genes.
#'
#' @param obj ExpressionSet object.
#' @param phenotype phenotype column name in the phenoData slot to check.
#' @param controlGenes Name of controlGenes, ie. 'Y' chromosome. Can specify 'all'.
#' @param columnID Column name where controlGenes is defined in the featureData slot if other than 'all'.
#' @param plotFlag TRUE/FALSE Whether to plot or not
#' @param legendPosition Location for the legend.
#' @param ... Extra parameters for \code{\link{plotCMDS}} function.
#'
#' @importFrom graphics legend
#'
#' @return Plots a classical multi-dimensional scaling of the 'controlGenes'. Optionally returns co-ordinates.
#' @export
#'
#' @examples
#' data(bladder)
#' checkMisAnnotation(bladder,'GENDER',controlGenes='Y',legendPosition='topleft')
#'
checkMisAnnotation <- function(obj, phenotype, controlGenes = "all",
                               columnID = "chromosome_name", plotFlag = TRUE, legendPosition = NULL,
                               ...) {
  if (tolower(controlGenes) != "all") {
    obj <- filterGenes(obj, labels = controlGenes, featureName = columnID,
                       keepOnly = TRUE)
  }
  if (length(phenotype) == 1) {
    phenotype <- factor(pData(obj)[, phenotype])
  }
  res <- plotCMDS(obj, pch = 21, bg = phenotype, plotFlag = plotFlag,
                  ...)
  if (!is.null(legendPosition))
    legend(legendPosition, legend = levels(phenotype), fill = 1:length(levels(phenotype)))
  invisible(res)
}
