#' Download GTEx files and turn them into ExpressionSet object
#'
#' Downloads the V6 GTEx release and turns it into an ExpressionSet object.
#'
#' @param type Type of counts to download - default genes.
#' @param file File path and name to automatically save the downloaded GTEx expression set. Saves as a RDS file.
#'
#' @return Organized ExpressionSet set.
#' @export
#'
#' @importFrom downloader download
#' @importFrom readr read_tsv
#' @importFrom Biobase AnnotatedDataFrame
#'
downloadGTEx<-function(type="genes",file=NULL,...){
  phenoFile = "http://www.gtexportal.org/static/datasets/gtex_analysis_v6/annotations/GTEx_Data_V6_Annotations_SampleAttributesDS.txt"
  pheno2File= "http://www.gtexportal.org/static/datasets/gtex_analysis_v6/annotations/GTEx_Data_V6_Annotations_SubjectPhenotypesDS.txt"
  geneFile  = "http://www.gtexportal.org/static/datasets/gtex_analysis_v6/rna_seq_data/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct.gz"

  message("Downloading and reading files")
  pdFile = tempfile("phenodat",fileext=".txt")
  download(phenoFile,destfile=pdFile)
  pd = read_tsv(pdFile)
  rownames(pd) = pd[,"SAMPID"]
  ids = sapply(strsplit(pd[,"SAMPID"],"-"),function(i)paste(i[1:2],collapse="-"))

  pd2File = tempfile("phenodat2",fileext=".txt")
  download(pheno2File,destfile=pd2File)
  pd2 = read_tsv(pd2File)
  rownames(pd2) = pd2[,"SUBJID"]
  pd2 = pd2[ids,]

  pdfinal = AnnotatedDataFrame(data.frame(cbind(pd,pd2)))

  if(type=="genes"){
    countsFile = tempfile("counts",fileext=".gz")
    download(geneFile,destfile=countsFile)
    cnts = suppressWarnings(read_tsv(geneFile,skip=2))
    genes = cnts[,1]
    geneNames = cnts[,2]
    counts = cnts[,-c(1:2)]
    rownames(counts) = genes
    for(i in 1:nrow(problems(cnts))){
      counts[problems(cnts)$row[i],problems(cnts)$col[i]] = 100000
    }
    throwAway = which(rowSums(counts)==0)
    counts = counts[-throwAway,]
    genes <- sub("\\..*","",rownames(counts))

    host="dec2013.archive.ensembl.org"; biomart="ENSEMBL_MART_ENSEMBL"; dataset="hsapiens_gene_ensembl"
    attributes=c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position","gene_biotype")
  }

  message("Creating ExpressionSet")
  pdfinal = pdfinal[match(colnames(counts),rownames(pdfinal)),]
  es = ExpressionSet(as.matrix(counts))
  phenoData(es) = pdfinal
  pData(es)["GTEX-YF7O-2326-101833-SM-5CVN9","SMTS"] = "Skin"
  pData(es)["GTEX-YEC3-1426-101806-SM-5PNXX","SMTS"] = "Stomach"

  message("Annotating from biomaRt")
  es = annotateFromBiomart(obj=es,genes=genes,host=host,biomart=biomart,dataset=dataset,attributes=attributes)

  message("Cleaning up files")
  unlink(pdFile)
  unlink(pd2File)
  unlink(countsFile)

  if(!is.null(file)) saveRDS(es,file=file)
  return(es)
}
