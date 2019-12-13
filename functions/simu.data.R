##########################################################################################
#  functions to simulate network data and network visualization used for HARMONIES
#  maintaniner: Shuang Jiang, <shuangj@smu.edu>
##########################################################################################

# Description:
#' Generate the network with adjacent matrix from the Erdős–Rényi model

# Dependency:
#' @import igraph

# Parameters:
#' @param p Number of nodes in the network
#' @param K Number of the major clusters in the network (K <= p).
#' @param edge.prob Probability of generating an edge between the two nodes
#' @param seed Random seed

# Output:
#' @return return.list A list of the binary adjacent matrix of the p nodes and the 
#'         table of the source and destination nodes

get.adjmatrix = function(p,
                         K = 2, 
                         edge.prob = 0.1,
                         seed = 123){
  # check input #
  if(p <= 0 | p!=floor(p)){
    stop("Number of nodes p must be a positive integer")
  }
  if(K > p){
    stop("Number of clusters K should not be larger than the number of nodes p")
  }else if(K <= 0 | K != floor(K)){
    stop("Number of clusters K must be a positive integer")
  }
  if(edge.prob <=0 | edge.prob>= 1){
    stop("edge.prob should be a number between 0-1")
  }
  seed.max = .Machine$integer.max
  if(seed < 0 | seed != floor(seed) | seed > seed.max){
    stop("random seed should be a positive integer with range of 1-",seed.max)
  }
  
  # load package 
  if (!require(igraph)) {install.packages("igraph", version = "1.2.4.1")}
  
  # initialization
  set.seed(seed)
  grp.index = floor(seq(0, p, length.out = K + 1));
  
  # generate the list of edge pairs from the erdos renyi model
  for (k in 1:K) {
    module = erdos.renyi.game(grp.index[k + 1] - grp.index[k], p.or.m = edge.prob, directed = FALSE);
    V(module)$name = as.character((grp.index[k] + 1):grp.index[k + 1]);
    if (k == 1) {
      g = module;
    } else {
      g = igraph::union(g, module, byname = "auto");
    }
  }
  if (K > 1) {
    for (k in 1:(K - 1)) {
      for (kk in (k + 1):K) {
        g = g + edge(as.character(sample((grp.index[k] + 1):grp.index[k + 1], 1)), as.character(sample((grp.index[kk] + 1):grp.index[kk + 1], 1)))
      }
    }
  }
  edge.l = data.frame(as_edgelist(g, names = TRUE), stringsAsFactors = FALSE);
  edge.l = rbind(edge.l, data.frame(cbind(1:p, 1:p), stringsAsFactors = FALSE));
  names(edge.l) = c("Source", "Target")
  
  # create the adjacent matrix 
  nodes.l = sort(unique(c(edge.l$Source, edge.l$Target)));
  
  # remove self-connected nodes
  sc.idx = which(edge.l$Source == edge.l$Target); 
  if (length(sc.idx) > 0) {
    edge.l = edge.l[-sc.idx,]; 
  }
  
  # Build the network
  G = make_empty_graph(directed = FALSE) + vertices(nodes.l);
  for (i in 1:nrow(edge.l)) {
    G = add_edges(G, c(as.character(edge.l[i, 1]), as.character(edge.l[i, 2])));
  }
  
  # Create the adjacent matrix
  adj = as.matrix(as_adjacency_matrix(G, type = "both"));
  node.order = as.numeric(rownames(adj))
  adj = adj[order(node.order), order(node.order)]
  return.list = list(edge.list = edge.l,
                     adj.matrix = adj)
  return (return.list)
}


# Description:
#' Generate the precision matrix and partial correlation matrix of a given network

# Parameters:
#' @param adj.mat A binary (0, 1) adjacent matrix. Can assign node names to the colunm or row names of \code{adj.mat}
#' @param seed Random seed

# Output:
#' @return Sigma Variance-covariance matrix
#' @return Omega Precision matrix
#' @return Corr Correlation matrix
#' @return Pcorr Partial correlation matrix

simu.network = function(adj.matrix,
                        seed = 123){
  # check input 
  if(ncol(adj.matrix) != ncol(adj.matrix)){
    stop("Input adj must be a square matrix ")
  }
  if(length(unique(as.vector(adj.matrix))) != 2 | any(unique(as.vector(adj.matrix)) != c(0, 1))){
    stop("Input adj must be a binary matrix with value 0 or 1")
  }
  seed.max = .Machine$integer.max
  if(seed < 0 | seed != floor(seed) | seed > seed.max){
    stop("random seed should be a positive integer with range of 1-",seed.max)
  }
  
  # initialization
  diag(adj.matrix) = 0
  set.seed(seed)
  p = nrow(adj.matrix)
  if(is.null(rownames(adj.matrix)) & is.null(colnames(adj.matrix))){
    names = seq(1, p)
  }else{
    name.idx = which(c(!is.null(rownames(adj.matrix)), !is.null(colnames(adj.matrix))))[1]
    if(name.idx == 1){
      names = rownames(adj.matrix)
    }else{
      names = colnames(adj.matrix)
    }
  }
  Omega = matrix(0, nrow = p, ncol = p);
  if (p >= 1000) {
    # set a larger value for the diagnal for a high dimension case
    diag(Omega) = 2;
  } else {
    diag(Omega) = 1;
  }
  
  # simulate a valid precision matrix 
  # Step 1
  edge.idx = which(adj.matrix == 1, arr.ind = TRUE)
  edge.idx = edge.idx[which(edge.idx[, 1] < edge.idx[, 2]),]
  edge.sign = sample(c(-1, 1), nrow(edge.idx), replace = TRUE) * runif(nrow(edge.idx), 0, 0.1);
  Omega[edge.idx] = edge.sign;
  Omega[cbind(edge.idx[, 2], edge.idx[, 1])] = edge.sign;
  
  # Step 2
  A = matrix(0, nrow = p, ncol = p);
  diag(A) = diag(Omega);
  for (j in 1:p) {
    if (sum(adj.matrix[j,]) > 0) {
      A[j, -j] = Omega[j, -j]/1.5/sum(abs(Omega[j, -j]))
    }
  }
  A = (A + t(A))/2;
  
  # Step 3
  temp = which(A > 0 & A < 0.1);
  if (length(temp) > 0) {
    A[temp] = 0.1;
  }
  temp = which(A < 0 & A > -0.1);
  if (length(temp) > 0) {
    A[temp] = -0.1;
  }
  rm(temp)
  
  # Create the covariance, correlation and partial correlation matrics
  Ar = solve(A);
  Sigma = Ar; 
  Par = matrix(0, nrow = p, ncol = p); 
  Cor = matrix(0, nrow = p, ncol = p); 
  for (j in 1:p) {
    for (jj in 1:p) {
      Cor[j, jj] = Ar[j, jj]/sqrt(Ar[j, j]*Ar[jj, jj]);
      Par[j, jj] = A[j, jj]/sqrt(A[j, j]*A[jj, jj]);
    }
  }
  
  
  rownames(Sigma) = names;
  colnames(Sigma) = names;
  rownames(A) = names;
  colnames(A) = names;
  rownames(Cor) = names;
  colnames(Cor) = names;
  rownames(Par) = names;
  colnames(Par) = names;
  
  Corr = Cor; 
  diag(Corr) = 1
  Pcorr = -Par; 
  diag(Pcorr) = 1
  
  return(list(Sigma = Sigma, 
              Omega = A, 
              Corr = Corr, 
              Pcorr = Pcorr));
  
}

# Description:
#' Visualize the simulated network in D3

# Dependency:
#' @import igraph
#' @importFrom networkD3 forceNetwork
#' @import plyr

# Parameters:
#' @param pcor.mat Simulated partial correlation matrix

# Output:
#' @return A D3 plot showing the simulated network with edge colored by the sign of the corresponding 
#'         nonzero partical correlation coefficient.

visualize.networkD3 = function(pcor.mat){
  # check input 
  if(!is.matrix(pcor.mat)){
    stop("pcor.mat must be a matrix")
  }
  if(!isSymmetric(pcor.mat)){
    stop("pcor.mat must be a symmetric matrix")
  }
  if(is.null(rownames(pcor.mat)) & is.null(colnames(pcor.mat))){
    nodes = paste0("Taxon", seq(1, nrow(pcor.mat)) )
  }else{
    name.idx = which(c(!is.null(rownames(pcor.mat)), !is.null(colnames(pcor.mat))))[1]
    if(name.idx == 1){
      nodes = rownames(pcor.mat)
    }else{
      nodes = colnames(pcor.mat)
    }
  }
  rownames(pcor.mat) = colnames(pcor.mat) = nodes
  
  # load package 
  if (!require(igraph)) {install.packages("networkD3", version = "1.2.4.1")}
  if (!require(networkD3)) {install.packages("networkD3", version = "0.4")}
  if (!require(plyr)) {install.packages("plyr", version = "1.8.4")}
  
  # get the partical correlation coefficients for the edge pairs
  p = nrow(pcor.mat)
  pcor.df = NULL;
  for (j in 1:(p - 1)) {
    for (jj in (j + 1):p) {
      if (abs(pcor.mat[j, jj]) > 0) {
        pcor.df = rbind(pcor.df, c(nodes[j], nodes[jj], pcor.mat[j, jj]))
      }
    }
  }
  edgeList = as.data.frame(pcor.df, stringsAsFactors = F)
  colnames(edgeList) <- c("SourceName", "TargetName", "Weight")
  edgeList$Weight = as.numeric(edgeList$Weight)
  edges.col <- ifelse(edgeList$Weight < 0, "FF0000", "00CC00")
  
  # Create a graph
  gD <- igraph::simplify(igraph::graph.data.frame(edgeList, directed=FALSE))
  
  # Create a data frame of information about nodes
  nodeList <- data.frame(ID = c(0:(igraph::vcount(gD) - 1)),
                         nName = igraph::V(gD)$name)
  nodeList$Size = 2
  nodeList$Group = 1
  
  # Map node names from the edge list to node IDs
  getNodeID <- function(x){which(x == igraph::V(gD)$name) - 1}
  edgeList <- plyr::ddply(edgeList, .variables = c("SourceName", "TargetName", "Weight"), 
                          function (x) data.frame(SourceID = getNodeID(x$SourceName), 
                                                  TargetID = getNodeID(x$TargetName)))
  
  # rescale the edge weight to a range of 2-6 for a better visualization 
  edge.weight = abs(as.numeric(edgeList$Weight))
  edge.weight = (edge.weight - min(edge.weight)) / max(edge.weight - min(edge.weight)) * (6 - 2) + 2
  edgeList$Weight = edge.weight
  
  D3_network_LM <- networkD3::forceNetwork(Links = edgeList, 
                                           Nodes = nodeList, 
                                           Source = "SourceID", 
                                           Target = "TargetID", 
                                           Value = "Weight",
                                           NodeID = "nName",
                                           Nodesize = "Size", 
                                           Group = "Group",  
                                           height = 500, 
                                           width = 1000, 
                                           fontSize = 15,
                                           opacity = 0.85,
                                           zoom = F, 
                                           opacityNoHover = 0.8, 
                                           bounded = T,
                                           linkColour = edges.col) 
  
  # Plot network
  D3_network_LM 
}


# Description:
#' Generate the network with adjacent matrix from the Erdős–Rényi model

# Dependency:
#' @import mvnfast

# Parameters:
#' @param n Number of random sample
#' @param Sigma Variance-covariance matrix of the underlying multivariate normal distribution
#' @param mu The global mean of the underlying multivariate normal distribution
#' @param zi Zero-inflation proportion, i.e. additional porportion of zeros added to the count matrix
#' @param seed Random seed

# Output:
#' @return Y Observed count matrix 
#' @return Y.full The count matrix without adding the inflated zeros
#' @return Mu.mat The true underlying multivariate normal data 
#' @return Sigma Variance-covariance matrix of the underlying multivariate normal distribution

get.countmat = function(n, 
                        Sigma, 
                        mu = 7,
                        zi = 0.1, 
                        seed = 123) {
  # check input 
  source("check.input.R")
  if(n != floor(n) | n<=0 ){
    stop("sample size n must be a psoitive integer")
  }
  if(!is.matrix(Sigma)){
    stop("Sigma should be a matrix")
  }
  if(!isSymmetric(Sigma)){
    stop("Sigma must be symmetric")
  }
  if( mu < 5 | mu > 10){
    warning("recommended value for mu is 5 - 10 to mimic the real microbiome data")
  }
  if(zi < 0 | zi > 1){
    stop("zi should be a proportion between 0-1")
  }
  seed = check.seed(seed)
  
  
  # load package 
  if (!require(mvnfast)) {install.packages("mvnfast", version = "0.2.5")}
  
  
  set.seed(seed);
  p = dim(Sigma)[1];
  # Fixed integers 
  # 1) range of total counts across samples
  N.LOWER = 50000
  N.UPPER = 100000
  N = floor(runif(n, N.LOWER, N.UPPER))
  # 2) half of the interval length to generate the mean vector of the underlying multivariate normal distribution
  RANGE = 5
  mu = runif(p, mu - RANGE, mu + RANGE);
  
  # generate the true underlying multivariateb normal data
  A = matrix(0, nrow = n, ncol = p);
  for(i in 1:n) {
    A[i,]= mvnfast::rmvn(1, mu, Sigma, ncores = 4);
  }
  A.true = A 
  A = exp(A);
  
  # generate Phi
  Phi = matrix(NA, nrow = n, ncol = p);
  for (i in 1:n) {
    Phi[i,] = rdirichlet(A[i,]);
  }
  
  # generate the observable counts
  Y = matrix(NA, nrow = n, ncol = p);
  for (i in 1:n) {
    Y[i,] = rmultinom(1, N[i], Phi[i,]);
  }
  
  # add the inflated zeros
  Y.obs = Y
  for(i in 1:n){
    for(j in 1:p){
      zi.idx = rbinom(1, 1, zi)
      Y.obs[i, j] = ifelse(zi.idx == 0, Y.obs[i, j], 0)
    }
  }

  return(list(Y = Y.obs, 
              Y.full = Y,
              Mu.mat = A.true, 
              Sigma = Sigma))
}

# Description:
#' Generate random sample from a Dirichlet distribution with parameter vector alpha

# Parameters:
#' @param alpha Concentration parameter vector for the Dirichlet distribution

# Output:
#' @return A random sample from a Dirichlet distribution with parameter vector alpha


rdirichlet = function(alpha) {
  # check input
  if(length(alpha) < 2 | any(alpha < 0)){
    stop("Concentration parameter of the Dirichlet distribution is invalid")
  }
  
  # generate random sample
  K = length(alpha)
  temp = rep(NA, K)
  for (j in 1:K) {
    temp[j] <- rgamma(1, alpha[j], 1)
  }
  temp <- temp/sum(temp)
  return (temp)
}
