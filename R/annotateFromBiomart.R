#' Annotate your Expression Set with biomaRt
#'
#' @param obj ExpressionSet object.
#' @param genes Genes or rownames of the ExpressionSet.
#' @param filters getBM filter value, see getBM help file.
#' @param attributes getBM attributes value, see getBM help file.
#' @param biomart BioMart database name you want to connect to. Possible database names can be retrieved with teh function listMarts.
#' @param dataset Dataset you want to use. To see the different datasets available within a biomaRt you can e.g. do: mart = useMart('ensembl'), followed by listDatasets(mart).
#' @param ... Values for useMart, see useMart help file.
#'
#' @return ExpressionSet object with a fuller featureData.
#' @export
#'
#' @importFrom biomaRt useMart
#' @importFrom biomaRt getBM
#' @importFrom Biobase featureNames
#' @importFrom Biobase fData
#' @importFrom Biobase fData<-
#' @importFrom Biobase ExpressionSet
#' @importClassesFrom Biobase ExpressionSet
#'
#' @examples
#'
#' data(skin)
#' # subsetting and changing column name just for a silly example
#' skin <- skin[1:10,]
#' colnames(fData(skin)) = paste("names",1:6)
#' biomart<-"ENSEMBL_MART_ENSEMBL";
#' genes <- sapply(strsplit(rownames(skin),split="\\."),function(i)i[1])
#' newskin <-annotateFromBiomart(skin,genes=genes,biomar=biomart,host=host)
#' head(fData(newskin)[,7:11])
#'
annotateFromBiomart <- function(obj,genes=featureNames(obj),filters="ensembl_gene_id",
                                attributes=c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position"),
                                biomart="ensembl",dataset="hsapiens_gene_ensembl",...){
  mart <- useMart(biomart=biomart,dataset=dataset,...)
  anno <- getBM(attributes = attributes, filters = filters,
                values = genes, mart = mart)
  if(nrow(anno)<length(genes)){
    warning("getBM returned fewer rows than genes queried.")
  }
  if(nrow(anno)>length(genes)){
    warning(sprintf("getBM returned more rows than genes queried. Using first call of %s.",colnames(anno)[1]))
    throw = which(duplicated(anno[,1]))
    anno  = anno[-throw,]
  }
  anno = anno[match(genes,anno[,"ensembl_gene_id"]),]
  if(!is.null(fData(obj))) anno = cbind(fData(obj),anno)
  fData(obj) = anno
  obj
}
