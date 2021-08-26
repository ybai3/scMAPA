#' Estimate significance of APA event at gene level
#'
#' @description  \code{APAtest} estimates significance of APA event at gene level.
#'
#' @param countMatrix Isoform-specific count matrix returned by \code{readinPAsites}.
#' @param coverageCutoff The cutoff for coverage of each gene (sum of long and short isoforms of all clusters). Default value is 20.
#' @param ORcutoff The cutoff for odds ratio of long isoform from each cluster against the grand mean of all clusters. Default to 0.25.
#' @param adPval The cutoff for FDR-controlled P values. Default to 0.05.
#' @return A data.frame containing significant APA genes with 1) their gene-level FDR controlled P values. 2) estimated coefficients (log(OR)), P values from Wald test, SE, and OR for each gene.
#'
#' @details There are two criteria to be met to call a gene to be APA. 1) FDR-controlled p value < 0.05. 2) At least one cluster whose odds of long isoform divided by odds of long isoform of all clusters > OR cutoff.
#'
#'
#' @examples
#' sigAPA_table  <- APAtest(countMatrix=Matix_returned_from_readin, coverageCutoff = 20, ORcutoff = 0.25, adPval = 0.05)
#'
APAtest <- function(countMatrix, coverageCutoff = 20, ORcutoff = 0.25, adPval = 0.05){
  PDUI_rd <- countMatrix
  pvalue_LRT <- numeric()
  EP <- matrix(ncol = ncol(PDUI_rd) * 4/2)
  for (i in 1:nrow(PDUI_rd)) {
    test <- as.vector(as.matrix(PDUI_rd[i, ]))
    test[which(is.na(test))] <- 0
    for (t in 1:(length(test)/2)){
      if (sum(test[c(t*2-1,t*2)])<coverageCutoff){
        test[c(t*2-1,t*2)] <- c(0,0)
      }
    }
    if (sum(test != 0)>=4){
      g1 <- data.frame(isoform = c(rep(1, test[1]), rep(0, test[2])), cluster = c(rep("C1", sum(test[c(1, 2)]))))
      for (k in 2:(ncol(PDUI_rd)/2)) {
        g1 <- rbind(g1, data.frame(isoform = c(rep(1, test[(2 * k - 1)]), rep(0, test[(2 * k)])), cluster = c(rep(paste0("C", k), sum(test[c((2 * k - 1), (2 * k))])))))
      }
      g1$cluster <- factor(g1$cluster)
      c_ref <- rownames(as.matrix(table(g1$cluster)))[1:2]
      contrasts(g1$cluster) <- contr.wec(g1$cluster, omitted = c_ref[1])
      greg <- glm(isoform ~ cluster, family = "binomial", data = g1)
      greg_mat <- summary(greg)$coefficients
      greg_mat <- cbind(greg_mat, exp(coef(greg)))
      colnames(greg_mat)[ncol(greg_mat)] <- "OR"
      contrasts(g1$cluster) <- contr.wec(g1$cluster, omitted = c_ref[2])
      greg_a <- glm(isoform ~ cluster, family = "binomial", data = g1)
      greg_a_mat <- summary(greg_a)$coefficients
      greg_a_mat <- cbind(greg_a_mat, exp(coef(greg_a)))
      colnames(greg_a_mat)[ncol(greg_a_mat)] <- "OR"
      greg_sum <- rbind(greg_a_mat, greg_mat[!(rownames(greg_mat) %in% rownames(greg_a_mat)), ])
      rownames(greg_sum)[nrow(greg_sum)] <- rownames(greg_mat)[!(rownames(greg_mat) %in% rownames(greg_a_mat))]
      greg_null <- glm(isoform ~ 1, family = "binomial", data = g1)
      pvalue_LRT <- c(pvalue_LRT, lrtest(greg, greg_null)$`Pr(>Chisq)`[2])
      names(pvalue_LRT)[length(pvalue_LRT)] <- rownames(PDUI_rd)[i]
      estimate_p <- numeric(length = ncol(PDUI_rd) * 4/2)
      cluster_name <- unique(str_remove(colnames(PDUI_rd), "_long_exp|_short_exp"))
      names(estimate_p) <- c(paste0(cluster_name, ".coef"), paste0(cluster_name, ".pval"), paste0(cluster_name, ".se"),paste0(cluster_name, ".OR"))
      for (j in 1:(ncol(PDUI_rd)/2)) {
        c <- greg_sum[which(rownames(greg_sum) == paste0("clusterC", j)), 1]
        estimate_p[j] <- ifelse(length(c) == 0, NA,c)
        p <- greg_sum[which(rownames(greg_sum) == paste0("clusterC", j)), 4]
        estimate_p[(j + ncol(PDUI_rd)/2)] <- ifelse(length(p) ==  0, NA, p)
        s <- greg_sum[which(rownames(greg_sum) == paste0("clusterC", j)), 2]
        estimate_p[(j + ncol(PDUI_rd))] <- ifelse(length(s) == 0, NA, s)
        r <- greg_sum[which(rownames(greg_sum) == paste0("clusterC", j)), 5]
        estimate_p[(j + ncol(PDUI_rd)*3/2)] <- ifelse(length(r) == 0, NA, r)
      }
      EP <- rbind(EP, estimate_p)
      rownames(EP)[length(rownames(EP))] <- rownames(PDUI_rd)[i]
    }
  }
  EP <- EP[-1, ]
  EP <- data.frame(Genes = rownames(EP), EP)
  pval_LRT.adjusted <- p.adjust(pvalue_LRT, method = "fdr")

  EP <- cbind(pvalue_LRT,pval_LRT.adjusted, EP)
  EP.OR <- EP[,str_detect(colnames(EP), ".OR")]  ## number need to be corrected
  EP.OR <- abs(1-EP.OR)
  sum(rowSums(EP.OR > ORcutoff, na.rm = T) > 0)
  EP <- EP[rowSums(EP.OR > ORcutoff, na.rm = T) > 0,]
  EP.sig <- EP[EP$pval_LRT.adjusted < adPval,]
  return(EP.sig)
}

#' Identify gene-cluster-specific 3'UTR shortening or lengthening
#'
#' @description \code{IdentifyClusterAPA} Identifies gene-cluster-specific 3'UTR shortening or lengthening based on ECoeffSig_Mat output from \code{APAtest}.
#'
#' @param ECoeffSig_Mat significant APA table returned from \code{APAtest}.
#' @param WaldP_cutoff The cutoff for P values of Wald tests. Default to 0.05.
#' @param CoeffCutoff The cutoff for estimated coefficients of logistic regression. Default to log(2).
#' @return \code{IdentifyClusterAPA} returns a list consists of tables for every cluster. Tables contain the gene IDs, APA event IDs of 3'UTR shortening or lengthening that passed filters.
#'
#' @details After identifying APA genes, scMAPA can identify which cluster contribute to the significant APA dynamic across clusters and estimate what degree, what direction and the significance of the cluster-specific APA event. There are two criteria to be met for a gene to be considered as gene-cluster-specific APA gene. 1) Gene with P values of Wald test on a certain cluster smaller than WaldP_cutoff. 2) Gene with absoluste estimated coeffcient for same cluster greater than CoeffCutoff. Direction will be assigned based on sign of coefficients.
#'
#' @examples
#' result_list  <- IdentifyClusterAPA(ECoeffSig_Mat = sigAPA_table, WaldP_cutoff=0.05, CoeffCutoff=log(2))
#'
IdentifyClusterAPA <- function (ECoeffSig_Mat, WaldP_cutoff = 0.05, CoeffCutoff = log(2))
{
  ECoeffSig_Mat.pval <- ECoeffSig_Mat[, str_detect(colnames(ECoeffSig_Mat), ".pval")]
  if (str_detect(ECoeffSig_Mat$Genes[1],"ENSMUSG")){
    ECoeffSig_Mat$GeneSymbol <- mapIds(org.Mm.eg.db, keys = str_extract(ECoeffSig_Mat$Genes, "ENSMUSG[:digit:]*"), keytype = "ENSEMBL", column="SYMBOL")
  } else if (str_detect(ECoeffSig_Mat$Genes[1],"ENSG")){
    ECoeffSig_Mat$GeneSymbol <- mapIds(org.Hs.eg.db, keys = str_extract(ECoeffSig_Mat$Genes, "ENSG[:digit:]*"), keytype = "ENSEMBL", column="SYMBOL")
  } else {
    stop("The species of gene ID is not ENSEMBL human or mouse.")
  }
  ECoeffSig_Mat <- ECoeffSig_Mat[,-c(1,2,3)]
  list_clusterAPA <- list()
  length(list_clusterAPA) <- ncol(ECoeffSig_Mat.pval)
  for (j in 1:ncol(ECoeffSig_Mat.pval)) {
    list_clusterAPA[[j]] <- data.frame(GeneSymbol = c(ECoeffSig_Mat$GeneSymbol[ECoeffSig_Mat[, j] > CoeffCutoff & ECoeffSig_Mat.pval[, j] <= WaldP_cutoff], ECoeffSig_Mat$GeneSymbol[ECoeffSig_Mat[, j] < (-CoeffCutoff) &  ECoeffSig_Mat.pval[, j] <= WaldP_cutoff]),
                                       UTR = c(rep("lengthening", length(ECoeffSig_Mat$GeneSymbol[ECoeffSig_Mat[, j] > CoeffCutoff & ECoeffSig_Mat.pval[, j] <= WaldP_cutoff])), rep("shortening", length(ECoeffSig_Mat$GeneSymbol[ECoeffSig_Mat[, j] < (-CoeffCutoff) & ECoeffSig_Mat.pval[, j] <= WaldP_cutoff]))))
    list_clusterAPA[[j]] <- list_clusterAPA[[j]][!is.na(list_clusterAPA[[j]]$GeneSymbol), ]
  }
  names(list_clusterAPA) <- paste0(unique(do.call(rbind, str_split(colnames(ECoeffSig_Mat.pval), "\\."))[, 1]))
  return(list_clusterAPA)
}
