#' Check for wrong annotation of a sample using classical MDS and control genes.
#'
#' @param obj ExpressionSet object.
#' @param phenotype phenotype column name in the phenoData slot to check.
#' @param controlGenes Name of controlGenes, ie. 'Y' chromosome. Can specify 'all'.
#' @param columnID Column name where controlGenes is defined in the featureData slot if other than 'all'.
#' @param legendPosition Location for the legend.
#' @param ... Extra parameters for plotMDS function.
#'
#' @return Plots a classical multi-dimensional scaling of the 'controlGenes'. Optionally returns co-ordinates.
#' @export
#'
#' @examples
#' \donttest{data(bladder);
#' checkMisAnnotation(bladder,"GENDER",controlGenes="Y",legendPosition="topleft");}
checkMisAnnotation <- function(obj,phenotype,controlGenes="all",
                                columnID="chromosome_name",legendPosition=NULL,...){
  if(is.null(controlGenes) | is.na(controlGenes) | controlGenes=="all" | controlGenes == "ALL"){
    keepGenes = 1:nrow(obj)
  } else {
    keepGenes = which(fData(obj)[,columnID]==controlGenes)
  }
  obj = obj[keepGenes,]
  if(length(phenotype)==1){
    pd = factor(pData(obj)[,phenotype])
  } else {
    pd = phenotype
  }
  res = plotCMDS(obj,pch=21,bg=pd,...)
  if(!is.null(legendPosition)) legend(legendPosition,legend=levels(pd),fill=1:length(levels(pd)))
  invisible(res)
}
