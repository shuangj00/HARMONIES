##########################################################################################
#  utility functions used for HARMONIES 
#  maintaniner: Shuang Jiang, <shuangj@smu.edu>
##########################################################################################

# Description: 
#' summarize the output from function "run.ZINBDPP"

# Parameters:
#' @param mcmc.output Output from function "run.ZINBDPP" 
#' @param count.matrix: original n-by-p count matrix where n is the sample size and p is the number of taxa
#' @param phenotype: phenotype indicator vector (all 0 or binary (0, 1) if have 2 groups)
#' @param bayes.fdr Bayesian false discovery rate controled in the zero imputation step

# Output:
#' @return nm.alpha0 Normalized abundance matrix for group 0 (the default group)
#' @return taxa.remove A binary vector indicating if a taxon is removed for fitting the ZINB-DPP model
#' @return nm.alpha1 (Optional) Normalized abundance matrix for 1. Only available if we have 2 phenotypes
#' @return gamma.ppi (Optional) posterior probability of inclusion of the differential abundance indicator gamma. 
#'         Only available if we have 2 phenotypes
#' @return fold.change (Optional) Fold change of the normalized abundance of group 1 over group 0 for the taxa included in 
#'         fitting the ZINB-DPP model. Only available if we have 2 phenotypes
#' 
summarize.zinbdpp = function(mcmc.output, 
                             count.matrix, 
                             phenotype,
                             bayes.fdr = 0.01){
  # check input 
  source("check.input.R")
  if(! all(count.matrix == floor(count.matrix)) | any(count.matrix < 0)){
    stop("Elements in the input count matrix must be nonnegative integers")
  }
  phenotype = check.phenotype(phenotype)
  if(bayes.fdr < 0 | bayes.fdr > 0.1){
    stop("bayes.fdr controls the FDR of zero imputation (suggested value: 0-0.1)")
  }
  
  # summarize the mcmc output 
  alpha.mat = mcmc.output$alpha.matrix
  taxa.remove = mcmc.output$remove.idx
  # remove the taxa not included by the ZINB-DPP fitting process
  alpha.mat = alpha.mat[, taxa.remove == 0]
  # get the normalized abundance
  nm.alpha = log(alpha.mat)
  # zero imputation
  H.mat = mcmc.output$H.ppi
  H.mat = H.mat[, taxa.remove == 0]
  impute.cutoff = get.BayFDR(as.vector(H.mat[H.mat!=0]), bayes.fdr)
  nm.alpha[count.matrix[, taxa.remove == 0] == 0 & H.mat < impute.cutoff] = 0
  
  # get diffrerential abundance information estimation if we have 2 groups 
  if(length(unique(phenotype)) == 2){
    gamma.ppi = mcmc.output$gamma.ppi[taxa.remove == 0]
    fold.change = apply(mcmc.output$fold.change, 2, mean)[taxa.remove == 0]
    nm.alpha0 = nm.alpha[phenotype == 0, ]
    nm.alpha1 = nm.alpha[phenotype == 1, ]
    # return a list #
    summary.list = list(nm.alpha0 = nm.alpha0,
                        taxa.remove = taxa.remove,
                        nm.alpha1 = nm.alpha1,
                        gamma.ppi = gamma.ppi,
                        fold.change = fold.change)
  }else{
    summary.list = list(nm.alpha0 = nm.alpha, 
                        taxa.remove = taxa.remove)
  }
  return(summary.list)
}

# Description:
#' Estimate a sparse partial correlation matrix

# Parameters:
#' @param alpha.matrix Normalized relative abundance matrix used as the input for prescison matrix estimation 
#' @param beta.stars Stability parameter used in StARS method
#' @param n.rep Number of subsamples used in StARS method
#' @param seed Random seed

# Output:
#' @return pcorr Estimated partial correlation matrix

est.pcor = function(alpha.matrix, 
                    beta.stars = 0.1, 
                    n.rep = 20, 
                    seed = 123){
  # load libraries 
  if (!require(qgraph)) {install.packages("qgraph", version = "1.6.3")}
  if (!require(huge)) {install.packages("huge", version = "1.3.2")}
  if (!require(pulsar)) {install.packages("pulsar", version = "0.3.5")}
  
  # initialize the input for glasso
  cor.mat = cor(alpha.matrix)
  diag(cor.mat) = 0
  lmd.grid = seq(max(cor.mat), max(cor.mat)/100, by = -0.01)
  hugeargs = list(lambda = lmd.grid, method = "glasso", verbose = F)
  
  # infer the sparse network by glasso
  cat("Estimating the sparse precision matrix... \n")
  pulsar.output = pulsar(alpha.matrix, 
                         fun=huge, 
                         fargs=hugeargs,
                         rep.num = n.rep,
                         thresh = beta.stars,
                         criterion='stars', 
                         lb.stars=TRUE, 
                         ub.stars=TRUE,
                         seed = seed)
  pulsar.fit  = refit(pulsar.output)
  
  # obtain the estimated precision matrix
  icov.refit = pulsar.fit[["est"]][["icov"]]
  sel.num = unlist(base::lapply(icov.refit,function(x){sum(x!=0) - ncol(alpha.matrix)} ))
  lam.ind.sel = grep(sum(pulsar.fit[["refit"]]$stars), sel.num)[1]
  precision.est = icov.refit[[lam.ind.sel]]
  
  # convert the precion matrix to partial correlation matrix
  pcorr = qgraph::wi2net(precision.est)
  pcorr = as.matrix(pcorr)
  return(pcorr)
}

# Description:
#' Extract the information (source, destination and partial correlation) from a partial correlation matrix

# Parameters:
#' @param pcorr A partial correlation matrix
#' @param taxa.name A vector of the taxa names for the correlation matrix. The length of \code{taxa.name} should match with the dimension
#'                  of the square matrix \code{pcorr}

# Output:
#' @return edge.df A data frame of source, destination and partial correlation of the taxa pairs that are partially correlated

summarize.edge = function(pcorr, 
                          taxa.name){
  # check input 
  if(nrow(pcorr) != ncol(pcorr)){
    stop("pcorr must be a square matrix")
  }
  if(nrow(pcorr) != length(taxa.name)){
    stop("Length of taxa.name should match with the dimension of pcorr")
  }
  
  # summarize the data frame of the edge information (source, destination and partial correlation)
  if(sum(abs(pcorr)) == 0){
    edge.df = data.frame(Source = character(), 
                         Destination = character(), 
                         PartialCorrelation= numeric(), stringsAsFactors = F)
    cat('No edge detected under the current parameter setting. Try with a larger stability parameter next time. \n')
  }else{
    # extract the nonzero partial correlation coefficients
    pcor.nonzero = which(pcorr != 0, arr.ind = T)
    
    # remove the duplicates
    rm.idx = apply(pcor.nonzero, 1, function(x){ifelse(x[2]>x[1], F, T)})
    pcor.ls = pcor.nonzero[!rm.idx, ]
    
    edge.df = t(apply(pcor.ls, 1, get.edge, tnames = taxa.name,
                      pcorr = pcorr))
    edge.df = data.frame(edge.df, stringsAsFactors = F)
    colnames(edge.df) = c("Source", "Destination", "PartialCorrelation")
  }
  return(edge.df)
}


# Description:
#' Merge the abundance information for the node from 2 phenotype groups. Additional information such as fold change, posterior
#' probability of inclusion will also be supplied

# Parameters:
#' @param taxa Name of taxa that we want to get the abundance information 
#' @param abundance.df0 A data frame of the taxa's abundance in group 0. Must contain the abundance of the input taxa 
#' @param abundance.df1 A data frame of the taxa's abundance in group 1. Must contain the abundance of the input taxa 
#' @param fold.change Estimated fold change from the ZINB-DPP model
#' @param gamma.ppi Posterior probability of inclusion of the discriminating indicator gamma obtained from the ZINB-DPP model

# Output:
#' @return abundance.df A data frame containing the estimated abundance of the input taxa in both groups, fold change and PPI

merge.abundance = function(taxa, 
                           abundance.df0,
                           abundance.df1, 
                           fold.change,
                           gamma.ppi){
  # check input #
  if(is.null(names(gamma.ppi)) | is.null(names(fold.change))){
    stop("mssing the name information in gamma.ppi and fold.change")
  }
  names.df0 = abundance.df0$Taxon
  names.df1 = abundance.df1$Taxon
  if(any(!(taxa %in% names.df0))){
    stop("missing the abundance information for the input taxa in group 0")
  }
  if(any(!(taxa %in% names.df1))){
    stop("missing the abundance information for the input taxa in group 1")
  }
  
  # get the abundance in group 0 and 1
  grp0.idx = sapply(taxa, function(x){which(abundance.df0$Taxon == x)})
  abundance.grp0 = abundance.df0[grp0.idx, ]
  grp1.idx = sapply(taxa, function(x){which(abundance.df1$Taxon == x)})
  abundance.grp1 = abundance.df1[grp1.idx, ]
  abundance.df = data.frame(Taxon = taxa, 
                            AbundanceGroup0 = abundance.grp0$Abundance, 
                            AbundanceGroup1 = abundance.grp1$Abundance, stringsAsFactors = F)
  
  # add fold change & PPI information #
  fc.idx = sapply(taxa, function(x){which(names(fold.change) == x)})
  fc.tmp = fold.change[fc.idx]
  ppi.idx = sapply(taxa, function(x){which(names(gamma.ppi) == x)})
  ppi.tmp = gamma.ppi[ppi.idx]
  abundance.df$FoldChange = fc.tmp
  abundance.df$PPI = ppi.tmp
  
  return(abundance.df)
}


# Description:
#' Get the cut-off value to control the Bayesian false discovery rate for the input posetrior probability of inclusion 

# Parameter:
#' @param PPI: posterior probability of inclusion obtained from the ZINB-DPP model
#' @param alpha: pre-specified Bayesian false discovery rate to be controlled

# Output:
#' @return cutoff A cut-off value that controls the Bayesian false discovery rate to be lower than alpha
#'
get.BayFDR <- function(PPI, alpha){
  PPI.sorted = sort(PPI,decreasing = TRUE)
  k = 1
  fdr = 0
  while(fdr < alpha){
    fdr = mean(1 - PPI.sorted[1:k])
    k = k+1
    if(k > length(PPI.sorted)){
      k = length(PPI.sorted)
      break
    }
  }
  cutoff = PPI.sorted[k]
  cutoff = ifelse(is.na(cutoff), 0, cutoff)
  return(cutoff)
}

# Description:
#' Obtain the log-scale abundance determined by both  the ZINB-DPP normalized abundance and the observed data

# Parameters:
#' @param taxa.names Name of taxa to obtain the abundance matrix
#' @param alpha.matrix Normalized abundance matrix from the ZINB-DPP model. The columns should match with \code{taxa.names}
#' @param count.matrix Orginal count matrix with the sample and taxa matching with \code{alpha.matrix}

# Output:
#' @return abundance.df A data frame of the input taxa names and their abundances

get.abundance = function(taxa.names, 
                         alpha.matrix, 
                         count.matrix){
  # check input 
  if(length(taxa.names) != ncol(alpha.matrix)){
    stop("Length of the taxa.names should be the same as the number of columns of alpha.matrix")
  }
  if(nrow(alpha.matrix) != nrow(count.matrix) | ncol(alpha.matrix) != ncol(count.matrix)){
    stop("Dimension of the alpha.matrix and count.matrix should be the same")
  }
  
  # calculate the abundance for the input taxa
  alpha.matrix[count.matrix == 0] = 0
  logA.vec = apply(alpha.matrix, 2, mean)
  abundance.df = data.frame(Taxon = taxa.names, Abundance = logA.vec,stringsAsFactors = F)
  return(abundance.df)
}


# Description:
#' Filter the taxa by the proportion of zeros across all the samples

# Parameter:
#' @param count.matrix A sample-by-taxon count matrix
#' @param sparsity.cutoff A threshold between 0-1. Taxa with proportions of zeros larger than the threshold
#'                        will be dropped for the network inference
                      
# Output:
#' @return filter.idx A vector containing the indices of filtered taxa 
#' 
filter.sparsetaxa = function(count.matrix, 
                             sparsity.cutoff = 0.5){
  zero.count.by.taxa = apply(count.matrix, 2, function(x){sum(x == 0)})
  gsize = nrow(count.matrix)
  drop.size = ceiling(sparsity.cutoff * gsize) + 1
  filter.idx = sapply(zero.count.by.taxa, function(x){ifelse(x > drop.size , T, F)})
  which(filter.idx == T)
}

# Description:
#' Create a vector containing one record of the selected edge information in the partial correlation

# Parameters:
#' @param idx A vector of row and column indices of a nonzero partial correlation in the partial correlation matrix
#' @param tnames Taxa name vector
#' @param pcorr Partial correlation matrix 

# Output:
#' @return edge.info A vector of the name of source, destination and partial correlation for the none zero entry in the partial correlation matrix
#' 
get.edge = function(idx, 
                    tnames, 
                    pcorr){
  name1 = tnames[idx[1]]
  name2 = tnames[idx[2]]
  pcor.res = pcorr[idx[1], idx[2]]
  if(pcor.res == 0){
    stop("There is no edge between the selected taxa pair")
  }
  edge.vec = c(name1, name2, pcor.res)
  return(edge.vec)
}

