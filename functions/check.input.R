##########################################################################################
#  Functions for checking input used for HARMONIES 
#  maintaniner: Shuang Jiang, <shuangj@smu.edu>
##########################################################################################

# Description: 
#' Check the phenotype indicator vector input

# Parameters:
#' @param phenotype A phenotype indicator vector for all the samples (can have at most 2 types of phenotype)

# Output:
#' @return A vector of 0 if there is only one phenotype, or a binary (0, 1) vector if there are 2 phenotypes

check.phenotype = function(phenotype){
  if(length(unique(phenotype)) > 2 ){
    stop("Can handle at most 2 phenotypes")
  }
  if(length(unique(phenotype)) == 1 && unique(phenotype) != 0 ){
    phenotype = rep(0, length(phenotype))
  }
  if(length(unique(phenotype)) == 2 & any(unique(phenotype) != c(0, 1)) ){
    phenotype.input = phenotype
    phenotype = factor(phenotype, levels = c(unique(phenotype)[1], unique(phenotype)[2]), labels = c(0, 1) )
    phenotype = as.numeric(levels(phenotype))[phenotype]
    cat("samples labeled as", unique(phenotype.input)[1], "is converted to 0, and samples labeled as", unique(phenotype.input)[2], "is converted to 1. \n" )
  }
  return(phenotype)
}
  

# Description: 
#' Check the random seed input

# Parameters:
#' @param seed A random seed (positive integer)

# Output:
#' @return A random seed

check.seed = function(seed){
  seed.max = .Machine$integer.max
  if(seed < 0 | seed != floor(seed) | seed > seed.max){
    stop("random seed should be a positive integer with range of 1-",seed.max)
  }
  return(seed)
}

  
  