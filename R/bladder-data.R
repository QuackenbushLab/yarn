#' Bladder RNA-seq data from the GTEx consortium
#'
#' Bladder RNA-seq data from the GTEx consortium. V6 release.
#'
#' @docType data
#'
#' @usage data(bladder)
#'
#' @format An object of class \code{"ExpressionSet"}; see \code{\link[Biobase]{ExpressionSet}}.
#'
#' @keywords datasets
#'
#' @references GTEx Consortium, 2015. The Genotype-Tissue Expression (GTEx) pilot analysis: Multitissue gene regulation in humans. Science, 348(6235), pp.648-660.
#' (\href{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4547484/}{PubMed})
#'
#' @source GTEx Portal
#'
#' @return ExpressionSet object
#'
#' @examples
#' \donttest{data(bladder);
#' checkMissAnnotation(bladder);}
"bladder"
