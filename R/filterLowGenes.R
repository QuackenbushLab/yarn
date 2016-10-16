#' Filter genes that have less than a minimum threshold CPM for a given group/tissue
#'
#' @param obj ExpressionSet object.
#' @param groups Vector of labels for each sample or a column name of the phenoData slot.
#' for the ids to filter. Default is the column names.
#' @param threshold The minimum threshold for calling presence of a gene in a sample.
#' @param minSamples Minimum number of samples - defaults to half the minimum group size.
#' @param ... Options for \link[edgeR]{cpm}.
#' @seealso \link[edgeR]{cpm} function defined in the edgeR package.
#'
#' @return Filtered ExpressionSet object
#' @export
#'
#' @importFrom edgeR cpm
#' @importFrom Biobase exprs
#' @importFrom Biobase pData
#'
#' @examples
#' data(skin)
#' filterLowGenes(skin,'SMTSD')
#'
filterLowGenes <- function(obj, groups, threshold = 1, minSamples = NULL,
                           ...) {
  if (is.null(minSamples)) {
    if (length(groups) == 1) {
      minSamples <- min(table(pData(obj)[, groups]))/2
    } else {
      minSamples <- min(table(groups))/2
    }
  }
  counts <- cpm(exprs(obj), ...)
  keep <- rowSums(counts > threshold) >= minSamples
  obj <- obj[keep, ]
  obj
}
