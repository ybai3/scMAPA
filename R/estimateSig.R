#' Estimate significance of APA event at gene level
#'
#' \code{estimateSig} estimates significance of APA event at gene level using either model-based approach or test-based approach. Model-based approach will fit a logistic regression model for each gene and estimate significance using Likelihood ratio test. Test-based approach will condunct Fisher's exact test on all possible pairs of clusters (i.e. 1 vs. 2; 1 vs. 3; 2 vs. 3 if 3 clusters) for each gene and use FDR to adjust for multiple-testing issue. Test-based approach defines significant APA genes if they pass significance test in any pair of the clusters.
#'
#' @param IsoMatrix Isoform-specific count matrix returned by \code{readin}.
#' @param mode Character parameter. Either "test" or "model", indicating which approach should be used to estimate the significance of APA events. Default to model-based approach. Default to "model".
#' @param FDR_P_cutoff The cutoff for FDR-controlled P values. Default to 0.05.
#' @return Model-based mode will return a list object containing three elements. 1) siglist_FDRp contains significant APA event IDs and their gene-level FDR controlled P values. 2) ECoeff_Mat contains estimated coefficients, P values, and SE from logsitic regression model for all genes. 3) ECoeffSig_Mat contains estimated coefficients, P values, and SE from logsitic regression model for significant APA events only and will be used for visualization. \cr \cr Test-based mode will return a list of IDs of significant APA events.
#'

#' @examples
#' result_list  <- estimateSig(ISOMatrix=Matix_returned_from_readin)
#' estimateSig(ISOMatrix=Matrix , mode = "model", FDR_P_cutoff=0.05)
#' estimateSig(ISOMatrix=Matrix , mode = "test", FDR_P_cutoff=0.05)
#'
estimateSig <- function(ISOMatrix , mode = "model", FDR_P_cutoff=0.05){
  ## check input parameters
  if (!(mode=='model'||mode=='test')) stop("'mode' must be set to 'model' or 'test'.")
  if (!is.numeric(FDR_P_cutoff)) stop("'FDR_P_cutoff' is not numeric.")

  ## model-based approach
  if (mode == 'model'){
    allshort <- character()
    alllong <- character()
    pvalue_cluster <- numeric()
    pvalue_indep <- numeric()
    EP <- matrix(ncol = ncol(ISOMatrix)*3/2)
    for(i in 1:nrow(ISOMatrix)){
      test <- as.vector(as.matrix(ISOMatrix[i,]))
      ## assign 0 to NA clusters
      test[which(is.na(test))] <- 0
      longsum <- sum(test[seq(from=1,to=length(test),by=2)])
      shortsum <- sum(test[seq(from=2,to=length(test),by=2)])
      if (longsum == 0){
        allshort <- c(allshort, rownames(ISOMatrix)[i])
      } else if(shortsum == 0){
        alllong <- c(alllong, rownames(ISOMatrix)[i])
      } else{
        g1 <- data.frame(isoform = c(rep(1, test[1]),rep(0, test[2])), cluster=c(rep("C1",sum(test[c(1,2)]))))
        for(k in 2:(ncol(ISOMatrix)/2)){
          g1 <- rbind(g1, data.frame(isoform = c(rep(1, test[(2*k-1)]),rep(0, test[(2*k)])), cluster=c(rep(paste0("C",k),sum(test[c((2*k-1),(2*k))])))))
        }
        g1$cluster <- factor(g1$cluster)
        c_ref <- rownames(as.matrix(table(g1$cluster)))[1:2]
        ## weighted average contrast
        contrasts(g1$cluster) <- contr.wec(g1$cluster, omitted = c_ref[1])
        greg <- glm(isoform~cluster, family="binomial", data = g1)
        greg_mat <- summary(greg)$coefficients
        contrasts(g1$cluster) <- contr.wec(g1$cluster, omitted = c_ref[2])
        greg_a <- glm(isoform~cluster, family="binomial", data = g1)
        greg_a_mat <- summary(greg_a)$coefficients

        greg_sum <- rbind(greg_a_mat ,greg_mat[!(rownames(greg_mat) %in% rownames(greg_a_mat)),])
        rownames(greg_sum)[nrow(greg_sum)] <- rownames(greg_mat)[!(rownames(greg_mat) %in% rownames(greg_a_mat))]

        greg_null <- glm(isoform~1, family="binomial", data = g1)
        pvalue_cluster <- c(pvalue_cluster,lrtest(greg,greg_null)$`Pr(>Chisq)`[2])
        names(pvalue_cluster)[length(pvalue_cluster)] <- rownames(ISOMatrix)[i]
        ## output estimated coefficients, pval and standard error.
        estimate_p <- numeric(length = ncol(ISOMatrix)*3/2)
        names(estimate_p) <- c(paste0(unique(do.call(rbind,str_split(colnames(ISOMatrix), "_"))[,1]),".coef"),
                               paste0(unique(do.call(rbind,str_split(colnames(ISOMatrix), "_"))[,1]),".pval"),
                               paste0(unique(do.call(rbind,str_split(colnames(ISOMatrix), "_"))[,1]),".se"))
        for(j in 1:(ncol(ISOMatrix)/2)){
          c <- greg_sum[which(rownames(greg_sum)==paste0("clusterC",j)),1]
          estimate_p[j] <- ifelse(length(c)==0, NA, c)
          ## assign Ci esimate
          p <- greg_sum[which(rownames(greg_sum)==paste0("clusterC",j)),4]
          estimate_p[(j+ncol(ISOMatrix)/2)] <- ifelse(length(p)==0, NA, p)
          ## assign Ci p-value
          s <- greg_sum[which(rownames(greg_sum)==paste0("clusterC",j)),2]
          estimate_p[(j+ncol(ISOMatrix))] <- ifelse(length(s)==0, NA, s)
          ## assign Ci standard error
        }
        EP <- rbind(EP, estimate_p)
        rownames(EP)[length(rownames(EP))] <- rownames(ISOMatrix)[i]
      }
    }
    EP <- EP[-1,]
    EP <- data.frame(Genes=rownames(EP), EP)
    pvalue_cluster <- pvalue_cluster[-1]
    pval.adjusted <- p.adjust(pvalue_cluster,method = "fdr")
    sum(pval.adjusted<=FDR_P_cutoff)
    siglist_scMAPA <- pval.adjusted[pval.adjusted<=FDR_P_cutoff]
    EP.sig <- EP[rownames(EP) %in% names(pval.adjusted)[pval.adjusted<FDR_P_cutoff],]
    return(list(siglist_FDRp=siglist_scMAPA, ECoeff_Mat=EP, ECoeffSig_Mat=EP.sig))
  }

  ## test-based approach
  if(mode == "test"){
    pval_naiveFDR <- matrix(nrow = nrow(ISOMatrix), ncol = choose((ncol(ISOMatrix)/2), 2))
    for (k in 1:nrow(ISOMatrix)){
      comp_slot =1
      for (i in 1:((ncol(ISOMatrix)/2)-1)){
        for (j in (i+1):(ncol(ISOMatrix)/2)){
          test <- as.integer(ISOMatrix[k,c((i*2-1):(i*2),(j*2-1):(j*2))])
          test[is.na(test)] <- 0
          pval_naiveFDR[k,comp_slot] <- fisher.test(matrix(test, nrow=2),simulate.p.value = T)$p.value
          comp_slot =comp_slot +1
        }
      }
    }
    pval_naiveFDR_adjust <- matrix(p.adjust(as.numeric(pval_naiveFDR), method = "fdr"),nrow = nrow(ISOMatrix), ncol = choose((ncol(ISOMatrix)/2), 2))
    siglist_pairwiseDaPars <- rownames(ISOMatrix)[rowSums(pval_naiveFDR_adjust <= FDR_P_cutoff)>=1]
    return(siglist_pairwiseDaPars)
  }
}

#' Identify gene-cluster-specific 3'UTR shortening or lengthening
#'
#' \code{IdentifyClusterAPA} Identifies gene-cluster-specific 3'UTR shortening or lengthening based on ECoeffSig_Mat output from model mode \code{estimateSig}.
#'
#' @param ECoeffSig_Mat Matrix containing estimated coefficients, unadjusted P values of Wald tests on coefficients, and SE returned by model mode \code{estimateSig}.
#' @param FDR_P_cutoff The cutoff for FDR-controlled P values of Wald tests. Default to 0.05. Only cluster-specific APA events with FDR-controlled P values of Wald test smaller than this number will be considered as 3'UTR lengthening or shortening.
#' @param CoeffCutoff The cutoff for estimated coefficients of logistic regression. Default to log(2). Only cluster-specific APA events with absoluste value of estimated coeffcients greater than this cutoff will be considered as 3'UTR lengthening or shortening.
#' @return \code{IdentifyClusterAPA} returns a list consists of tables for every cluster. Tables contain the gene IDs, APA event IDs of 3'UTR shortening or lengthening that passed filters.
#'
#' @examples
#' result_list  <- IdentifyClusterAPA(ECoeffSig_Mat = result_from_estimateSig$ECoeffSig_Mat, FDR_P_cutoff=0.05, CoeffCutoff=log(2))
#'
IdentifyClusterAPA <- function(ECoeffSig_Mat, FDR_P_cutoff=0.05, CoeffCutoff=log(2)){
  ECoeffSig_Mat.pval <- ECoeffSig_Mat[,str_detect(colnames(ECoeffSig_Mat), ".pval")]
  EP.qval <- matrix(nrow = nrow(ECoeffSig_Mat),ncol = ncol(ECoeffSig_Mat.pval))
  for (i in 1:ncol(ECoeffSig_Mat.pval)){
    pval <- ECoeffSig_Mat.pval[,i]
    EP.qval[,i] <- p.adjust(pval)
  }
  ECoeffSig_Mat$GeneSymbol <- do.call(rbind, str_split(ECoeffSig_Mat$Genes, "\\|"))[,2]
  ECoeffSig_Mat$TranscriptID <- do.call(rbind, str_split(ECoeffSig_Mat$Genes, "\\|"))[,1]
  list_clusterAPA <- list()
  length(list_clusterAPA) <- ncol(ECoeffSig_Mat.pval)
  ECoeffSig_Mat <- ECoeffSig_Mat[,-1]
  for (j in 1:ncol(ECoeffSig_Mat.pval)){
    list_clusterAPA[[j]] <- data.frame(GeneSymbol=c(ECoeffSig_Mat$GeneSymbol[ECoeffSig_Mat[,j]>CoeffCutoff & EP.qval[,j]<=FDR_P_cutoff], ECoeffSig_Mat$GeneSymbol[ECoeffSig_Mat[,j]<(-CoeffCutoff) & EP.qval[,j]<=FDR_P_cutoff]), UTR=c(rep("lengthening", length(ECoeffSig_Mat$GeneSymbol[ECoeffSig_Mat[,j]>CoeffCutoff & EP.qval[,j]<=FDR_P_cutoff])), rep("shortening", length(ECoeffSig_Mat$GeneSymbol[ECoeffSig_Mat[,j]<(-CoeffCutoff) & EP.qval[,j]<=FDR_P_cutoff]))), TranscriptID= c(ECoeffSig_Mat$TranscriptID[ECoeffSig_Mat[,j]>CoeffCutoff & EP.qval[,j]<=FDR_P_cutoff], ECoeffSig_Mat$TranscriptID[ECoeffSig_Mat[,j]<(-CoeffCutoff) & EP.qval[,j]<=FDR_P_cutoff]))
    list_clusterAPA[[j]] <- list_clusterAPA[[j]][!is.na(list_clusterAPA[[j]]$GeneSymbol),]
  }
  names(list_clusterAPA) <- paste0(unique(do.call(rbind,str_split(colnames(ECoeffSig_Mat.pval), "\\."))[,1]))
  return(list_clusterAPA)
}
