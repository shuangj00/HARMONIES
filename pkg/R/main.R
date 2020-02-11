##########################################################################################
#  HARMONIES: A Hybrid Approach for Microbiome Networks Inference via Exploiting Sparsity
#  maintaniner: Shuang Jiang, <shuangj@smu.edu>
##########################################################################################

#' HARMONIES: a hybrid model to infer the microbiome network based on a sparse estimation
#' of the precision matrix.
#'
#' The model consists of two major components: (1) a zero-inflated
#' negative binomial model with Dirichlet process prior for microbiome count data normalization;
#' (2) the graphical lasso algorithm with a stability-based parameter tuning process to obtain a sparse network.

# Dependencies:
#' @import huge
#' @import pulsar
#' @import Rcpp
#' @import RcppArmadillo
#' @importFrom qgraph wi2net

# Parameters:
#' @param count.matrix An \eqn{n*p} sequencing count matrix.
#'                     Columns represent the taxonomic features and rows represents the samples.
#'                     Set the column name to be the name of taxa, and all the taxa here should be on the same taxonomic level.
#' @param phenotype A phenotype indicator vector for all the samples (should be a vector of 0s if only have
#'                  one phenotype, or a vector of 0 and 1 if have 2 phenotypes. The current method can handle
#'                  at most 2 phenotypes)
#' @param N.mcmc Number of MCMC iterations
#' @param b Shape hyper-parameter for the variance term
#' @param h Scale hyper-parameter for the variance term
#' @param sparsity.cutoff A threshold between 0-1. Taxa with proportions of zeros larger than the threshold
#'                        will be dropped for the network inference
#' @param beta.stars Stability threshold for selection criteria
#' @param n.rep Number of random subsamples to take for network estimation
#' @param bayes.fdr Bayesian false discovery rate controlled for zero imputation (in the normalization step)
#' @param seed: Random seed

# Outputs:
#' @return A list with the following results:
#' \itemize{
#'   \item partial.corr: A list of the partical correlation estimations for group 0 and group 1 (only for two phenotypes) by HARMONIES
#'   \item edge.estimation: A list of data frames for group 0 and group 1 (only for two phenotypes) containing the edge information of the estimated network
#'   \item node.estimation: A data frame containing the node information of all the nodes shown in the estimated network(s)
#' }
#' @export
HARMONIES = function(count.matrix,
                     phenotype,
                     N.mcmc = 10000,
                     b = 1,
                     h = 20,
                     sparsity.cutoff = 0.5,
                     beta.stars = 0.05,
                     n.rep = 20,
                     bayes.fdr = 0.05,
                     seed = 123
){
  # check input
  if(! all(count.matrix == floor(count.matrix)) | any(count.matrix < 0)){
    stop("Elements in the input count matrix must be nonnegative integers")
  }
  phenotype = check.phenotype(phenotype)
  seed = check.seed(seed)
  if(b <= 0 | h <= 0){
    stop("parameter 'b' and 'h' must be strictly positive")
  }
  if(N.mcmc != floor(N.mcmc) ){
    stop("N.mcmc should be a large positive integer to ensure convergence (suggested value: >= 10000)")
  }
  if(sparsity.cutoff > 1 | sparsity.cutoff < 0){
    stop("sparsity.cutoff must be between 0 and 1 (suggested value: 0.5)")
  }
  if(beta.stars < 0 | beta.stars > 1){
    stop("The instability parameter beta.stars must be between 0-1 (suggested value: 0.05)")
  }
  if(n.rep <0 | n.rep != floor(n.rep) ){
    stop("Number of subsamples in StARS must be a positive integer (suggested value: 20)")
  }
  if(bayes.fdr < 0 | bayes.fdr > 0.1){
    stop("bayes.fdr controls the FDR of zero imputation (suggested value: 0-0.1)")
  }


  # load libraries
  # if (!require(qgraph)) {install.packages("qgraph", version = "1.6.3")}
  # if (!require(huge)) {install.packages("huge", version = "1.3.2")}
  # if (!require(pulsar)) {install.packages("pulsar", version = "0.3.5")}
  # if (!require(Rcpp)) {install.packages("Rcpp", version = "1.0.2")}


  # fixed parameters
  # minimun number of nonzero observations needed for each taxon for fitting the ZINB-DPP model
  COUNT.MIN = 2;


  # run MCMC
  mcmc.output = run.ZINBDPP(count.matrix = count.matrix,
                            phenotype = phenotype,
                            N.mcmc = N.mcmc,
                            b = b,
                            h = h,
                            count.min = COUNT.MIN,
                            seed = seed)


  # infer network
  if(is.null(colnames(count.matrix))){
    taxa.names = paste0("taxon", seq(1, ncol(count.matrix)))
  }else{
    taxa.names = colnames(count.matrix)
  }
  HARMONIES.res.list = get.network(mcmc.output = mcmc.output,
                                   count.matrix = count.matrix,
                                   phenotype = phenotype,
                                   sparsity.cutoff = sparsity.cutoff,
                                   taxa.names = taxa.names,
                                   beta.stars = beta.stars,
                                   n.rep = n.rep,
                                   bayes.fdr = bayes.fdr,
                                   seed = seed)
  return(HARMONIES.res.list)
}



#' Data Normalization by ZINB-DPP
#'
#' Description
#' run MCMC algorithm to fit the zero-inflated negative binomial (ZINB) model with Dirichlet process prior (DPP) for the input n-by-p count matrix

# Dependencies:
#' @import Rcpp

# Parameters:
#' @param count.matrix An \eqn{n*p} sequencing count matrix.
#'                     Columns represent the taxonomic features and rows represents the samples.
#'                     Set the column name to be the name of taxa, and all the taxa here should be on the same taxonomic level.
#' @param phenotype A phenotype indicator vector for all the samples (should be a vector of 0s if only have
#'                  one phenotype, or a vector of 0 and 1 if have 2 phenotypes. The current method can handle
#'                  at most 2 phenotypes)
#' @param N.mcmc Number of MCMC iterations
#' @param b Shape hyper-parameter for the variance term
#' @param h Scale hyper-parameter for the variance term
#' @param count.min: minimum number of none zero counts required for a taxon to fit the ZINB-DPP model
#' @param seed: Random seed

# Output:
#' @return A S3 object of class"mcmc.zinbdpp":
#' \itemize{
#'   \item FoldChange Fold change (group 1 over group 0) of the normalized abundance for all p taxa
#'   \item remove.idx A binary vector of length p that indicates if a taxon is dropped from fitting the ZINB-DPP model
#'   \item size.factor A matrix that stores the MCMC output (after burn-in the first half) for the size factor of all the sample
#'   \item alpha.matrix A matrix of the posterior mean of {alpha_ij}, i = 1, ..., n, and j = 1, ..., p
#'   \item phi A vector of the posterior mean of the dispersion parameter
#'   \item pi An n-by-p matrix of the posterior mean of the missing probability pi for each count in the input matrix
#'   \item H.ppi An n-by-p matrix of the posterior probability of inclusion (PPI) of being a missing value
#'         for each count in the input matrix
#'   \item gamma.ppi A vector of the posterior probability of inclusion (PPI) for each of p taxon to be discriminating between
#'         patient phenotypes (if we have 2 groups)
#'   \item gamma.accept.rate Acceptance rate for updating gamma in the Metropolis–Hastings algorithm
#'   \item phi.accept.rate Acceptance rate for updating phi in the Metropolis–Hastings algorithm
#'   \item alpha.accept.rate Acceptance rate for updating alpha in the Metropolis–Hastings algorithm
#' }

run.ZINBDPP = function(count.matrix = count.matrix,
                       phenotype = phenotype,
                       N.mcmc = N.mcmc,
                       b = b,
                       h = h,
                       count.min = COUNT.MIN,
                       seed = seed){
  # check input
  if(! all(count.matrix == floor(count.matrix)) | any(count.matrix < 0)){
    stop("Elements in the input count matrix must be nonnegative integers")
  }
  if(!is.matrix(count.matrix)){
    count.matrix = as.matrix(count.matrix)
  }
  phenotype = check.phenotype(phenotype)
  seed = check.seed(seed)
  if(b <= 0 | h <= 0){
    stop("parameter 'b' and 'h' must be strictly positive")
  }
  if(N.mcmc != floor(N.mcmc) ){
    stop("N.mcmc should be a large positive integer to ensure convergence (suggested value: >= 10000)")
  }


  # load package
  #if (!require(Rcpp)) {install.packages("Rcpp", version = "1.0.2")}


  # set initial values for the MCMC
  s0 = rep(1, nrow(count.matrix))
  S0 = matrix(1, 1, 1)


  # fit ZINB-DPP model
  set.seed(seed)
  cat("Start fitting the ZINB-DPP model. This may take a while... \n")
  MCMC.output  = fitZINBDPP(Y = count.matrix,
                            z = phenotype,
                            s = s0,
                            iter = N.mcmc,
                            DPP = TRUE,
                            S = S0,
                            aggregate = FALSE,
                            b = b,
                            h = h,
                            mincount = count.min,
                            MRF = FALSE,
                            G = S0)
  class(MCMC.output) = append(class(MCMC.output), "mcmc.zinbdpp")
  return(MCMC.output)
}

#' Estimate the sparse precision matrix by Glasso
#'
#' Description
#' Run graphical lasso on the normalized abundances given by the ZINB-DPP model. The tuning parameter is selected by StARS.

# Parameters:
#' @param mcmc.output MCMC output from function "run.MCMC.R"
#' @param count.matrix A count matrix from metagenomic shotgun sequencing or 16SrRNA sequencing technologies.
#'                  Columns represent the taxonomic features and rows represents the samples.
#' @param phenotype A phenotype indicator vector (can have at most 2 phenotypes)
#' @param sparsity.cutoff A threshold between 0-1. Taxa with proportions of zeros larger than the threshold
#'                        will be dropped for the network inference
#' @param taxa.names: Name of p taxa in the input count matrix for the ZINB-DPP model (default is NULL if using the simulated data)
#' @param beta.stars Stability parameter used in StARS method
#' @param n.rep Number of subsamples used in StARS method
#' @param bayes.fdr Bayesian false discovery rate controled in the normalization step
#' @param seed: Random seed

# Outputs:
#' @return A list with the following components:
#' \itemize{
#'  \item partial.corr: A list of the partical correlation estimations for group 0 and group 1(if have two phenotypes) by HARMONIES
#'  \item edge.estimation: A list of data frames for group 0 and group 1(if have two phenotypes) containing the edge information of the estimated network
#'  \item node.estimation: A data frame containing the node information of all the nodes shown in the estimated network(s)
#'    }

get.network = function(mcmc.output,
                       count.matrix,
                       phenotype,
                       sparsity.cutoff = 0.5,
                       taxa.names = NULL,
                       beta.stars = 0.1,
                       n.rep = 20,
                       bayes.fdr = 0.01,
                       seed = 123){
  # check input
  if(! inherits(mcmc.output, "mcmc.zinbdpp")){
    stop("mcmc.output must be the output from function 'run.MCMC.R' ")
  }
  if(! all(count.matrix == floor(count.matrix)) | any(count.matrix < 0)){
    stop("Elements in the input count matrix must be nonnegative integers")
  }
  phenotype = check.phenotype(phenotype)
  seed = check.seed(seed)
  if(sparsity.cutoff > 1 | sparsity.cutoff < 0){
    stop("sparsity.cutoff must be between 0 and 1 (suggested value: 0.5)")
  }
  if(beta.stars < 0 | beta.stars > 1){
    stop("The instability parameter beta.stars must be between 0-1 (suggested value: 0.05)")
  }
  if(n.rep <0 | n.rep != floor(n.rep) ){
    stop("Number of subsamples in StARS must be a positive integer (suggested value: 20)")
  }


  # load packages
  # if (!require(qgraph)) {install.packages("qgraph", version = "1.6.3")}
  # if (!require(huge)) {install.packages("huge", version = "1.3.2")}
  # if (!require(pulsar)) {install.packages("pulsar", version = "0.3.5")}
  #

  # load utility functions


  # step1. summarize the MCMC output from the ZINB-DPP model
  mcmc.summary = summarize.zinbdpp(mcmc.output = mcmc.output,
                                   count.matrix = count.matrix,
                                   phenotype = phenotype,
                                   bayes.fdr = bayes.fdr)
  rm(mcmc.output)

  # get the binary vector indicating the removed taxa (0: keep / 1: remove)
  taxa.remove = mcmc.summary$taxa.remove
  count.matrix = count.matrix[, taxa.remove == 0]


  # step2: filter the taxa by names and sparsity
  if(is.null(colnames(count.matrix))){
    taxa.names = paste0("taxon", seq(1, ncol(count.matrix)))
  }else{
    taxa.names = colnames(count.matrix)
  }
  taxa.names = taxa.names[taxa.remove == 0]
  rm.char = c("unclassified", "noname", "Virus", "virus")

  # filter by name
  rm.idx = sapply(taxa.names, function(x){any(sapply(rm.char, grepl,x))})
  rm.idx = which(rm.idx)

  # filter by sparsity
  sparse.idx = filter.sparsetaxa(count.matrix[phenotype == 0,],
                                 sparsity.cutoff = sparsity.cutoff)
  filter.idx = intersect(rm.idx, sparse.idx)
  keep.idx = seq(1, sum(taxa.remove == 0))
  if(length(filter.idx) != 0){
    keep.idx = keep.idx[-filter.idx]
  }

  # if have the additional phenotype
  if(length(unique(phenotype)) == 2){
    sparse.idx = filter.sparsetaxa(count.matrix[phenotype == 1,],
                                   sparsity.cutoff = sparsity.cutoff)
    filter.idx1 = intersect(rm.idx, sparse.idx)
    keep.idx1 = seq(1, sum(taxa.remove == 0))
    if(length(filter.idx1) != 0){
      keep.idx1 = keep.idx1[-filter.idx]
    }
  }


  # step3: estimated the partial correlation
  # (1) the default group
  pcorr.grp0 = est.pcor(alpha.matrix = as.matrix(mcmc.summary$nm.alpha0[, keep.idx]),
                        beta.stars = beta.stars,
                        n.rep = n.rep,
                        seed = seed)
  edge.grp0 = summarize.edge(pcorr = pcorr.grp0,
                             taxa.name = taxa.names[keep.idx]  )
  node.names0 = c()
  if(nrow(edge.grp0) == 0){
    cat('No detected edge for group labeled with 0 under the current parameter setting. \n')
    abundance.grp0 = data.frame(Taxon = character(),Abundance = numeric(), stringsAsFactors = F)
  }else{
    # abundance of all the taxa in group 0 (the default group)
    abundance.ref0 = get.abundance(taxa.names = taxa.names,
                                   alpha.matrix = mcmc.summary$nm.alpha0,
                                   count.matrix = count.matrix[phenotype == 0, ])
    node.names0 = unique(c(edge.grp0[, 1], edge.grp0[, 2]))

    # extract the abundance for the taxa in the network for group 0
    name.idx = sapply(node.names0, function(x){which(abundance.ref0[, 1] == x)})
    abundance.grp0 = abundance.ref0[name.idx, ]
  }

  if(length(unique(phenotype)) == 1){
    # the inference stops here since we only have a single group
    res.list = list(partial.corr = pcorr.grp0,
                    edge.estimation = edge.grp0,
                    node.estimation = abundance.grp0)
  }else{
    node.names1 = c()
    # (2) the additional group
    pcorr.grp1 = est.pcor(alpha.matrix = as.matrix(mcmc.summary$nm.alpha1[, keep.idx1]),
                          beta.stars = beta.stars,
                          n.rep = n.rep,
                          seed = seed)
    edge.grp1 = summarize.edge(pcorr = pcorr.grp1,
                               taxa.name = taxa.names[keep.idx1]  )
    node.names1 = unique(c(edge.grp1[, 1], edge.grp1[, 2]))

    # additional differential abundance information
    gamma.ppi = mcmc.summary$gamma.ppi
    fold.change = mcmc.summary$fold.change
    names(gamma.ppi) = names(fold.change) = taxa.names
    # abundance of all the taxa in group 1 (the additional group)
    abundance.ref1 = get.abundance(taxa.names = taxa.names,
                                   alpha.matrix = mcmc.summary$nm.alpha1,
                                   count.matrix = count.matrix[phenotype == 1, ])

    # create the output for the edge the node information
    if(nrow(edge.grp1) == 0){
      cat('No detected edge for group labeled with 1 under the current parameter setting. \n')
      abundance.grp1 = data.frame(Taxon = character(), Abundance = numeric(), stringsAsFactors = F)

      # update the data frame of the node infomation for group 0
      if(nrow(edge.grp0) == 0){
        # there is no egde in either networks, and the node information is empty
        abundance.df = data.frame(Taxon = character(), AbundanceGroup0 = numeric(), AbundanceGroup1 = numeric(), FoldChange = numeric(), PPI = numeric(), stringsAsFactors = F)
        res.list = list(partial.corr = list(Group0 = pcorr.grp0,
                                            Group1 = pcorr.grp1),
                        edge.estimation = list(Group0 = edge.grp0,
                                               Group1 = edge.grp1),
                        node.estimation = abundance.df)

      }else{
        # edge(s) in group 0 only
        abundance.df = merge.abundance(taxa =  node.names0,
                                       abundance.df0 = abundance.ref0,
                                       abundance.df1 = abundance.ref1,
                                       fold.change = fold.change,
                                       gamma.ppi = gamma.ppi)
        res.list = list(partial.corr = list(Group0 = pcorr.grp0,
                                            Group1 = pcorr.grp1),
                        edge.estimation = list(Group0 = edge.grp0,
                                               Group1 = edge.grp1),
                        node.estimation = abundance.df)
      }
    }else{
      # edge(s) in group 1 or both group 0 and 1
      if(nrow(edge.grp0) == 0){
        taxa.combine = node.names1
      }else{
        taxa.combine = unique( node.names0, node.names1)
      }
      abundance.df = merge.abundance(taxa = taxa.combine,
                                     abundance.df0 = abundance.ref0,
                                     abundance.df1 = abundance.ref1,
                                     fold.change = fold.change,
                                     gamma.ppi = gamma.ppi)
      res.list = list(partial.corr = list(Group0 = pcorr.grp0,
                                          Group1 = pcorr.grp1),
                      edge.estimation = list(Group0 = edge.grp0,
                                             Group1 = edge.grp1),
                      node.estimation = abundance.df)

    }
  }
  return(res.list)

}


