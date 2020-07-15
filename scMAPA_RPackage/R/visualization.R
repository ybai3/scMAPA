#' Heatmap of cluster-specific APA genes
#'
#' \code{clusterAPAheatmap} draws heatmap of cluster-specific APA genes to show the degree and direction of APA. This function first selects genes that are significant cluster-specific APA in any cluster. Then, it draws heatmap with estimated coefficients from logistic regression.
#'
#' @param ECoeffSig_Mat Matrix containing estimated coefficients, unadjusted P values of Wald tests on coefficients, and SE returned by model mode \code{estimateSig}.
#' @param FDR_P_cutoff The cutoff for FDR-controlled P values of Wald tests. Default to 0.05. Only cluster-specific APA events with FDR-controlled P values of Wald test smaller than this number will be considered as 3'UTR lengthening or shortening.
#' @param CoeffCutoff The cutoff for estimated coefficients of logistic regression. Only cluster-specific APA events with absoluste value of estimated coeffcients greater than this cutoff will be considered as 3'UTR lengthening or shortening.
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
#' Draw dot plot of genes selected by user. The size of dots shows the deviation of the proportion of long isoforms from grand mean of all transcripts. The color of dots shows the direction of 3' UTR processing (lengthening or shortening).
#'
#' @param ECoeffSig_Mat Matrix containing estimated coefficients, unadjusted P values of Wald tests on coefficients, and SE returned by model mode \code{estimateSig}.
#' @param FDR_P_cutoff The cutoff for FDR-controlled P values of Wald tests. Default to 0.05. Only cluster-specific APA events with FDR-controlled P values of Wald test smaller than this number will be considered as 3'UTR lengthening or shortening.
#' @param CoeffCutoff The cutoff for estimated coefficients of logistic regression. Only cluster-specific APA events with absoluste value of estimated coeffcients greater than this cutoff will be considered as 3'UTR lengthening or shortening.
#' @param APAgenes Character vector, a list of gene symbols/gene IDs user would like to display in the dot plot.
#' @return \code{APAdotplot} returns a dot plot showing the status of 3' UTR legnthening or shortening for user selected genes.
#'
#' @examples
#' APAdotplot(ECoeffSig_Mat, FDR_P_cutoff=0.05, CoeffCutoff=log(2), APAgenes=c("gene","symbols","gene","IDs"))
#'
APAdotplot <- function(ECoeffSig_Mat, FDR_P_cutoff=0.05, CoeffCutoff=log(2), APAgenes){
  if (!is.character(APAgenes)) stop("'APAgenes' is not character vector.")
  clustername <- unique(do.call(rbind, str_split(colnames(ECoeffSig_Mat), "\\."))[,1])[-1]
  ECoeffSig_Mat.pval <- ECoeffSig_Mat[,str_detect(colnames(ECoeffSig_Mat), ".pval")]
  EP.qval <- matrix(nrow = nrow(ECoeffSig_Mat),ncol = ncol(ECoeffSig_Mat.pval))
  for (i in 1:ncol(ECoeffSig_Mat.pval)){
    pval <- ECoeffSig_Mat.pval[,i]
    EP.qval[,i] <- p.adjust(pval)
  }
  EP.sig.list <- list()
  length(EP.sig.list) <- ((ncol(ECoeffSig_Mat)-1)/3)
    for (i in 1:length(EP.sig.list)){
      EP.sig.list[[i]] <- data.frame(Genes=ECoeffSig_Mat$Genes,coef=ECoeffSig_Mat[,(i+1)], qval=EP.qval[,i],cluster=rep(clustername[i], nrow(ECoeffSig_Mat)))
      EP.sig.list[[i]] <- EP.sig.list[[i]][(abs(ECoeffSig_Mat[,(i+1)])>CoeffCutoff)&(EP.qval[,i]<FDR_P_cutoff),]
    }
  EP.clus.dotplot <- do.call(rbind, EP.sig.list)
  EP.clus.dotplot$abs_Coef <- abs(EP.clus.dotplot$coef)
  EP.clus.dotplot$isoform <- ifelse(EP.clus.dotplot$coef>0, "long", "short")
  EP.clus.dotplot$Genes <- do.call(rbind,str_split(EP.clus.dotplot$Genes, "\\|"))[,2]
  EP.clus.dotplot <- EP.clus.dotplot[tolower(EP.clus.dotplot$Genes) %in% tolower(APAgenes),]
  ggplot(EP.clus.dotplot,aes(x=Genes, y=cluster, size=abs_Coef, color=isoform))+
    geom_point()+theme_classic()+theme(axis.text.x = element_text(angle = 45,hjust = 1))
}
