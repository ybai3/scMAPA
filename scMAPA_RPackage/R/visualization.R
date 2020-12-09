#' Heatmap of cluster-specific APA genes
#'
#' \code{clusterAPAheatmap} draws heatmap of cluster-specific APA genes to show the degree and direction of APA. This function first selects genes that are significant cluster-specific APA in any cluster. Then, it draws heatmap with estimated coefficients from logistic regression.
#'
#' @param ECoeffSig_Mat significant APA table returned from \code{APAtest}.
#' @param FDR_P_cutoff The cutoff for FDR-controlled P values of Wald tests. Default to 0.05.
#' @param CoeffCutoff The cutoff for estimated coefficients of logistic regression.
#' @return \code{clusterAPAheatmap} returns a heatmap of cluster-specific APA genes
#'
#' @examples
#' clusterAPAheatmap(ECoeffSig_Mat = result_from_estimateSig$ECoeffSig_Mat, FDR_P_cutoff=0.05, CoeffCutoff=log(2))
#'

clusterAPAheatmap <- function(ECoeffSig_Mat, FDR_P_cutoff=0.05, CoeffCutoff=log(2)){
  ECoeffSig_Mat.pval <- ECoeffSig_Mat[,str_detect(colnames(ECoeffSig_Mat), ".pval")]
  EP.qval <- matrix(nrow = nrow(ECoeffSig_Mat),ncol = ncol(ECoeffSig_Mat.pval))
  for (i in 1:ncol(ECoeffSig_Mat.pval)){
    pval <- ECoeffSig_Mat.pval[,i]
    EP.qval[,i] <- p.adjust(pval)
  }
  keep.EP.coef<-(abs(ECoeffSig_Mat[,str_detect(colnames(ECoeffSig_Mat), ".coef")])>=CoeffCutoff)>0
  keep.EP.pval<-(EP.qval<=FDR_P_cutoff)>0
  keep.EP <- keep.EP.coef+keep.EP.pval
  keep.EP <- rowSums(keep.EP==2, na.rm = T)>=1

  EP.sig.clus <- ECoeffSig_Mat[keep.EP,]
  y <- as.matrix(EP.sig.clus[,str_detect(colnames(ECoeffSig_Mat), ".coef")])
  y[is.na(y)] <- 0
  colnames(y) <- unique(do.call(rbind, str_split(colnames(ECoeffSig_Mat), "\\."))[,1])[-1]
  colfunc <- colorRampPalette(c("blue", "white", "red"))
  heatmap.2(y, col=colfunc(15),
            dendrogram="both", srtCol=45,
            scale="none", density.info="none", trace="none",labRow = NA, key.xlab = "Coefficient")
}

#' Dot plot of selected APA genes
#'
#' @description Draw dot plot of cluster-specific APA genes selected by user. The size of dots shows the deviation of the proportion of long isoforms from grand mean of all transcripts. The color of dots shows the direction of 3' UTR processing (lengthening or shortening).
#'
#' @param ECoeffSig_Mat significant APA table returned from \code{APAtest}.
#' @param FDR_P_cutoff The cutoff for FDR-controlled P values of Wald tests. Default to 0.05.
#' @param CoeffCutoff The cutoff for estimated coefficients of logistic regression.
#' @param APAgenes Character vector, a list of gene IDs user would like to display in the dot plot. The format should be consistent with geneID part of rownames in the significant APA table, transcript IDs, chrIDs, and all other parts of APA ID should be removed. For better interpretation, gene IDs will be converted to gene symbols in the plot.
#' @return \code{APAdotplot} returns a dot plot showing the status of 3' UTR legnthening or shortening for user selected genes.
#'
#' @examples
#' APAdotplot(ECoeffSig_Mat, FDR_P_cutoff=0.05, CoeffCutoff=log(2), APAgenes=str_extract(rownames(ECoeffSig_Mat), "ENSMUSG[:digit:]*"))
#'
APAdotplot <- function(ECoeffSig_Mat, FDR_P_cutoff=0.05, CoeffCutoff=log(2), APAgenes){
  if (!is.character(APAgenes)) stop("'APAgenes' is not character vector.")
  clustername <- unique(do.call(rbind, str_split(colnames(ECoeffSig_Mat), "\\."))[,1])[-1]
  clustername <- clustername[-c(1,2)]
  ECoeffSig_Mat.pval <- ECoeffSig_Mat[,str_detect(colnames(ECoeffSig_Mat), ".pval")]
  EP.qval <- matrix(nrow = nrow(ECoeffSig_Mat),ncol = ncol(ECoeffSig_Mat.pval))
  for (i in 1:ncol(ECoeffSig_Mat.pval)){
    pval <- ECoeffSig_Mat.pval[,i]
    EP.qval[,i] <- p.adjust(pval)
  }
  EP.sig.list <- list()
  length(EP.sig.list) <- ((ncol(ECoeffSig_Mat)-3)/4)
  ECoeffSig_Mat <- ECoeffSig_Mat[,-c(1,2)]
  for (i in 1:length(EP.sig.list)){
    EP.sig.list[[i]] <- data.frame(Genes=ECoeffSig_Mat$Genes,coef=ECoeffSig_Mat[,(i+1)], qval=EP.qval[,i],cluster=rep(clustername[i], nrow(ECoeffSig_Mat)))
    EP.sig.list[[i]] <- na.omit(EP.sig.list[[i]][(abs(ECoeffSig_Mat[,(i+1)])>CoeffCutoff)&(EP.qval[,i]<FDR_P_cutoff),])
  }
  EP.clus.dotplot <- do.call(rbind, EP.sig.list)
  EP.clus.dotplot$abs_Coef <- abs(EP.clus.dotplot$coef)
  EP.clus.dotplot$isoform <- ifelse(EP.clus.dotplot$coef>0, "long", "short")
  if (str_detect(EP.clus.dotplot$Genes, "ENSMUSG")[1]){
    EP.clus.dotplot$Genes <- str_extract(EP.clus.dotplot$Genes, "ENSMUSG[:digit:]*")
  } else if(str_detect(EP.clus.dotplot$Genes, "ENSG")[1]){
    EP.clus.dotplot$Genes <- str_extract(EP.clus.dotplot$Genes, "ENSG[:digit:]*")
  }else{
    stop("The input gene ID is not ENSEMBL human or mouse.")
  }
  EP.clus.dotplot <- EP.clus.dotplot[tolower(EP.clus.dotplot$Genes) %in% tolower(APAgenes),]
  if (str_detect(EP.clus.dotplot$Genes[1],"ENSMUSG")){
    EP.clus.dotplot$GeneSymbol <- mapIds(org.Mm.eg.db, keys = str_extract(EP.clus.dotplot$Genes, "ENSMUSG[:digit:]*"), keytype = "ENSEMBL", column="SYMBOL")
  } else if (str_detect(EP.clus.dotplot$Genes[1],"ENSG")){
    EP.clus.dotplot$GeneSymbol <- mapIds(org.Hs.eg.db, keys = str_extract(EP.clus.dotplot$Genes, "ENSG[:digit:]*"), keytype = "ENSEMBL", column="SYMBOL")
  } else {
    stop("The species of gene ID is not ENSEMBL human or mouse.")
  }
  ggplot(EP.clus.dotplot,aes(x=GeneSymbol, y=cluster, size=abs_Coef, color=isoform))+
    geom_point()+theme_classic()+theme(axis.text.x = element_text(angle = 45,hjust = 1))
}
