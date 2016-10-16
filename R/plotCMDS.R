#' Plot classical MDS of dataset
#'
#' This function plots the MDS coordinates for the "n" features of interest. Potentially uncovering batch
#' effects or feature relationships.
#'
#' @param obj ExpressionSet object or objrix.
#' @param comp Which components to display.
#' @param normalized TRUE / FALSE, use the normalized matrix or raw counts.
#' @param distFun Distance function, default is dist.
#' @param distMethod The distance measure to be used. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". Any unambiguous substring can be given.
#' @param n Number of features to make use of in calculating your distances.
#' @param samples Perform on samples or genes.
#' @param log TRUE/FALSE log2-transform raw counts.
#' @param plotFlag TRUE/FALSE whether to plot or not.
#' @param ... Additional plot arguments.
#' @return coordinates
#'
#' @importFrom matrixStats rowSds
#' @importFrom stats dist
#' @importFrom stats cmdscale
#' @importFrom graphics plot
#'
#' @export
#' @examples
#' data(skin)
#' res <- plotCMDS(skin,pch=21,bg=factor(pData(skin)$SMTSD))
#' \donttest{
#' # library(calibrate)
#' # textxy(X=res[,1],Y=res[,2],labs=rownames(res))
#' }
plotCMDS <- function(obj, comp = 1:2, normalized = FALSE, distFun = dist,
                     distMethod = "euclidian", n = NULL, samples = TRUE, log = TRUE,
                     plotFlag = TRUE, ...) {
  if (is.null(n))
    n <- min(nrow(obj), 1000)
  obj <- extractMatrix(obj, normalized, log)
  genesToKeep <- which(rowSums(obj) > 0)
  geneVars <- rowSds(obj[genesToKeep, ])
  geneIndices <- genesToKeep[order(geneVars, decreasing = TRUE)[seq_len(n)]]
  obj <- obj[geneIndices, ]

  if (samples == TRUE) {
    obj <- t(obj)
  }
  d <- distFun(obj, method = distMethod)
  ord <- cmdscale(d, k = max(comp))
  xl <- paste("MDS component:", comp[1])
  yl <- paste("MDS component:", comp[2])

  if (plotFlag == TRUE)
    plot(ord[, comp], ylab = yl, xlab = xl, ...)
  invisible(ord[, comp])
}
