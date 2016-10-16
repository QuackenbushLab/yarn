#' Plot heatmap of most variable genes
#'
#' This function plots a heatmap of the gene expressions forthe "n" features of interest.
#'
#' @param obj ExpressionSet object or objrix.
#' @param n Number of features to make use of in plotting heatmap.
#' @param fun Function to sort genes by, default \code{\link[stats]{sd}}.
#' @param normalized TRUE / FALSE, use the normalized matrix or raw counts.
#' @param log TRUE/FALSE log2-transform raw counts.
#' @param ... Additional plot arguments for \code{\link[gplots]{heatmap.2}}.
#' @return coordinates
#'
#' @importFrom gplots heatmap.2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats sd
#'
#' @export
#' @examples
#' data(skin)
#' tissues <- pData(skin)$SMTSD
#' plotHeatmap(skin,normalized=FALSE,log=TRUE,trace="none",n=10)
#' # Even prettier
#' \donttest{
#' # library(RColorBrewer)
#' data(skin)
#' tissues <- pData(skin)$SMTSD
#' heatmapColColors <- brewer.pal(12,"Set3")[as.integer(factor(tissues))]
#' heatmapCols <- colorRampPalette(brewer.pal(9, "RdBu"))(50)
#' plotHeatmap(skin,normalized=FALSE,log=TRUE,trace="none",n=10,
#'  col = heatmapCols,ColSideColors = heatmapColColors,cexRow = 0.6,cexCol = 0.6)
#'}
plotHeatmap <- function(obj, n = NULL, fun = stats::sd, normalized = TRUE,
                        log = TRUE, ...) {
  if (is.null(n))
    n <- min(nrow(obj), 100)
  mat <- extractMatrix(obj, normalized, log)
  genesToKeep <- which(rowSums(mat) > 0)
  geneStats <- apply(mat[genesToKeep, ], 1, fun)
  geneIndices <- genesToKeep[order(geneStats, decreasing = TRUE)[seq_len(n)]]
  mat <- mat[geneIndices, ]
  heatmap.2(mat, ...)
  invisible(mat)
}
