#' Skin RNA-seq data from the GTEx consortium
#'
#' Skin RNA-seq data from the GTEx consortium. V6 release. Random selection of 20 skin samples.
#' 13 of the samples are fibroblast cells, 5 Skin sun exposed, 2 sun unexposed.
#'
#' @docType data
#'
#' @usage data(skin)
#'
#' @format An object of class \code{"ExpressionSet"}; see \code{\link[Biobase]{ExpressionSet}}.
#'
#' @keywords datasets
#'
#' @return ExpressionSet object
#'
#' @references GTEx Consortium, 2015. The Genotype-Tissue Expression (GTEx) pilot analysis: Multitissue gene regulation in humans. Science, 348(6235), pp.648-660.
#' (\href{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4547484/}{PubMed})
#'
#' @source GTEx Portal
#'
#' @examples
#' \donttest{data(skin);
#' checkMissAnnotation(skin,"GENDER");}
"skin"
