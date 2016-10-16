#' Filter samples
#'
#' @param obj ExpressionSet object.
#' @param ids Names found within the groups labels corresponding to samples to be removed
#' @param groups Vector of labels for each sample or a column name of the phenoData slot
#' for the ids to filter. Default is the column names.
#' @param keepOnly Filter or keep only the samples with those labels.
#'
#' @return Filtered ExpressionSet object
#' @export
#'
#' @importFrom Biobase pData
#'
#' @examples
#' data(skin)
#' filterSamples(skin,ids = "Skin - Not Sun Exposed (Suprapubic)",groups="SMTSD")
#' filterSamples(skin,ids=c("GTEX-OHPL-0008-SM-4E3I9","GTEX-145MN-1526-SM-5SI9T"))
#'
filterSamples <- function(obj, ids, groups = colnames(obj), keepOnly = FALSE) {
  if (length(groups) == 1) {
    groups <- pData(obj)[, groups]
  }
  throwAway <- which(groups %in% ids)
  if (keepOnly) {
    obj <- obj[, throwAway]
  } else {
    obj <- obj[, -throwAway]
  }
  obj
}
