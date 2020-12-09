#' Generate filtered isofrom-specific count matrix
#'
#' @description Read in the output from DaPars2, select clusters of interests, filter out genes with low quality, and generate a isoform-specific count matrix for APA detection.
#'
#' @param path Path to the folder containing (only) outputs from DaPars for all clusters.
#' @param NAcutoff The lower limit for number of clusters in which the expression of a gene is detected (not NA).
#' @param CPMcutoff_L The gene-wise cutoff for average of CPM values of long isoforms across clusters. Default value is 10.
#' @param CPMcutoff_S The gene-wise cutoff for average of CPM values of short isoforms across clusters. Default value is 10.
#' @param clusterOfInterests The clusters should be included in the downstream analysis. It should be "all" to include all clusters or a vector of integers indicating the index of clusters of interest. At least 2 clusters should be provided.
#' @return The return will be an isoform specific count matrix based on the DaPars2 estimation.
#'
#' @details \code{readinPAsites} will read in the output from BAMprocess programs and do selection and filtering based on parameters described above.
#' @details To filter out genes with low quality, only genes that expressed in number of clusters >= NAcutoff will be kept. \code{readinPAsites} will calculate CPM for long and short Isoforms separately by CPM=(raw count)*10^6/lib_size. Then, only genes with rowMeans(CPM of Long isoform)>CPMcutoff_L & rowMeans(CPM of Short isoform)>CPMcutoff_S will be kept. The flexibility of setting different cutoffs for long and short isoforms allows user to customize the QC by the feature of their samples.
#'
#' @examples
#' ISOMatrix  <- readinPAsites(path=getwd(), NAcutoff=3, CPMcutoff_L = 10, CPMcutoff_S = 10, clusterOfInterests = "all")
#' ISOMatrix  <- readin(path="Directory/storing/result/tables/from/previous/step")
#'
#'
readinPAsites <- function(path, NAcutoff = 3, CPMcutoff_L = 10, CPMcutoff_S = 10, clusterOfInterests = "all"){
  if (!dir.exists(path))
    stop(paste("directory", path, "does not exist."))
  if (!is.numeric(NAcutoff))
    stop("'NAcutoff' is not numeric.")
  if (!(is.numeric(CPMcutoff_L) & is.numeric(CPMcutoff_S)))
    stop("'CPMcutoff' is not numeric.")
  if (CPMcutoff_L <= 0 | CPMcutoff_S <= 0)
    stop("'CPMcutoff' should be larger than 0, the default value is 10.")
  readin <- list.files(path = path)
  DaPars_list <- list()
  length(DaPars_list) <- length(readin)
  for (i in 1:length(readin)){
    DaPars_list[[i]] <- read.delim(readin[i], stringsAsFactors=FALSE)
    numClusters = (ncol(DaPars_list[[i]])-4)/3
    DaPars_list[[i]] <- DaPars_list[[i]][,sort(c(1,seq(5,5+(numClusters-1)*3,by=3), seq(6,6+(numClusters-1)*3,by=3)))]
  }
  names(DaPars_list) <- readin
  PDUI <- do.call(rbind,unname(DaPars_list))
  if (clusterOfInterests[1] != "all"){
    if (length(clusterOfInterests)<2)
      stop("At least 2 clusters should be provided.")
    ## if only two clusters are going to be compared, there should be no NAs in each row.
    if (length(clusterOfInterests)==2){
      NAcutoff = 2
    }
    if (!(class(clusterOfInterests[1]) %in% c("numeric","integer")))
      stop("'clusterOfInterests' should be a vector of integers indicating the index of clusters of interest.")
    select_clus <- integer()
    for (i in 1:length(clusterOfInterests)){
      select_clus = c(select_clus,clusterOfInterests[i]*2)
      select_clus = c(select_clus,clusterOfInterests[i]*2+1)
    }
    PDUI <- PDUI[,c(1, select_clus)]
  }
  ## filter by rowwise mean of CPM
  PDUI_cpm_L <- PDUI[,seq(2, ncol(PDUI), by=2)]
  PDUI_cpm_L[is.na(PDUI_cpm_L)] <- 0
  PDUI_cpm_S <- PDUI[,seq(3, ncol(PDUI), by=2)]
  PDUI_cpm_S[is.na(PDUI_cpm_S)] <- 0
  #PDUI_cpm <- data.frame(C1_exp=rowSums(PDUI[,c(2,3)], na.rm = T))
  #for (i in 2:((ncol(PDUI)-1)/2)){
  #  PDUI_cpm <- cbind(PDUI_cpm, rowSums(PDUI[,c((2*i),((2*i)+1))], na.rm = T))
  #  colnames(PDUI_cpm)[i] <- paste0("C",i,"_exp")
  #}
  lib_size_L <- colSums(PDUI_cpm_L)
  lib_size_S <- colSums(PDUI_cpm_S)
  #lib_size <- colSums(PDUI_cpm)
  PDUI_cpm_L <- t((t(PDUI_cpm_L) * 10^6)/lib_size_L)
  PDUI_cpm_S <- t((t(PDUI_cpm_S) * 10^6)/lib_size_S)
  #PDUI_cpm <- t((t(PDUI_cpm) * 10^6)/lib_size)

  #sum(rowMeans(PDUI_cpm_L) > 10)
  #sum(rowMeans(PDUI_cpm_S) > 10)
  #sum(rowMeans(PDUI_cpm_L) > 10 & rowMeans(PDUI_cpm_S) > 10)
  #sum(rowMeans(PDUI_cpm) > 10)
  PDUI_keep <- rowMeans(PDUI_cpm_L) > CPMcutoff_L & rowMeans(PDUI_cpm_S) > CPMcutoff_S
  #PDUI_keep <- rowMeans(PDUI_cpm)>CPMcutoff
  PDUI <- PDUI[PDUI_keep,]
  PDUI <- PDUI[rowSums(!(is.na(PDUI[,-1]))) >= (NAcutoff*2),]
  PDUI_rd <- round(PDUI[,-1])
  rownames(PDUI_rd) <- PDUI$Gene
  ##remove genes where only 1 cluster has very large counts so that it passes gene-wise filtering
  rowsum_PDUIrd <- rep(0, nrow(PDUI_rd))
  for (i in 1:(ncol(PDUI_rd)/2)){
    a <- rowSums(PDUI_rd[,c((i*2-1),i*2)],na.rm = T)>0
    rowsum_PDUIrd <- rowsum_PDUIrd + a
  }
  PDUI_rd <- PDUI_rd[!(rowsum_PDUIrd <= 1),]
  return(PDUI_rd)
}
