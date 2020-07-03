#' Read in DaPars outputs and filter out low expressed genes.
#'
#' \code{readin} returns filtered isoform count matrix for APA identification.
#'
#' @param path Path to the folder containing (only) outputs from DaPars for all samples.
#' @param NAcutoff The number of NAs could be tolerated for each gene. \code{readin} will only keep genes with number of NAs <= NAcutoff.
#' @param CPMcutoff The cutoff for gene-wise sum of CPM. \code{readin} will calculate CPM for each long/short Isoform from each sample separately by CPM=(raw count)*10^6/lib_size. Then, only genes with rowSums(CPM)>CPMcutoff will be kept.
#' @return The returns will be an isoform specific count matrix based on the DaPars estimation.
#'
#' @examples
#' ISOMatrix  <- readin(path=getwd(),NAcutoff=3)
#' ISOMatrix  <- readin(path="Directory/storing/result/tables/from/previous/step",NAcutoff=3,CPMcutoff=10)
#'
#'
readin <- function(path, NAcutoff, CPMcutoff=10){
  ## check input parameters
  if (!dir.exists(path)) stop(paste('directory',path,'does not exist.'))
  if (!is.numeric(NAcutoff)) stop("'NAcutoff' is not numeric.")
  if (!is.numeric(CPMcutoff)) stop("'CPMcutoff' is not numeric.")
  if (CPMcutoff<=0) stop("'CPMcutoff' should be larger than 0, the default value is 10.")
  if (NAcutoff<=0) stop("'NAcutoff' should be larger than 0.")

  readin_list <- list.files(path = path)
  DaPars_list <- list()
  length(DaPars_list) <- length(readin_list)
  for (i in 1:length(readin_list)){
    DaPars_list[[i]] <- read.delim(readin_list[i], stringsAsFactors=FALSE)
    DaPars_list[[i]] <- DaPars_list[[i]][,c(1,5,6)]
  }
  names(DaPars_list) <- readin_list
  PDUI <- suppressWarnings(Reduce(function(x,y) merge(x,y,1,all=T),DaPars_list))
  PDUI_cpm <- data.frame(C1_exp=rowSums(PDUI[,c(2,3)], na.rm = T))
  for (i in 2:((ncol(PDUI)-1)/2)){
    PDUI_cpm <- cbind(PDUI_cpm, rowSums(PDUI[,c((2*i),((2*i)+1))], na.rm = T))
    colnames(PDUI_cpm)[i] <- paste0("C",i,"_exp")
  }
  lib_size <- colSums(PDUI_cpm)
  PDUI_cpm <- t((t(PDUI_cpm) * 10^6)/lib_size)
  PDUI_keep <- rowSums(PDUI_cpm)> CPMcutoff
  ## CPM cutoff
  PDUI <- PDUI[PDUI_keep,]
  PDUI <- PDUI[rowSums(is.na(PDUI[,-1])) <= NAcutoff,]
  ## NA cutoff
  PDUI_rd <- round(PDUI[,-1])
  rownames(PDUI_rd) <- PDUI$Gene
  colnames(PDUI_rd) <- paste0(rep(do.call(rbind,str_split(readin_list,"_"))[,2],each=2), c("_long_exp","_short_exp"))
  return(PDUI_rd)
}
