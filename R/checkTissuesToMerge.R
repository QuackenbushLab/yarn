#' Check tissues to merge based on gene expression profile
#'
#' @param obj ExpressionSet object.
#' @param majorGroups Column name in the phenoData slot that describes the general body region or site of the sample.
#' @param minorGroups Column name in the phenoData slot that describes the specific body region or site of the sample.
#' @param ... Parameters that can go to \link{checkMissAnnotation}
#'
#' @return CMDS Plots of the majorGroupss colored by the minorGroupss. Optional matrix of CMDS loadings for each comparison.
#' @export
#'
#' @examples
#' data(skin)
#' checkTissuesToMerge(skin,"SMTS","SMTSD")
checkTissuesToMerge <- function(obj,majorGroups,minorGroups,...){
  region = pData(obj)[,majorGroups]
  result = sapply(unique(region),function(i){
    keepSamples = which(region == i)
    objSubset = obj[,keepSamples]
    objSubset = objSubset[which(rowSums(exprs(objSubset))>0),]
    res = checkMisAnnotation(objSubset,phenotype=minorGroups,controlGenes = "all",...)
    res
  })
  invisible(result)
}
