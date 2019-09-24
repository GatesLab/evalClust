#' Strength Preserving Correlation Matrix Randomizer
#'
#'This function implements the Masuda, Kojaku & Sano (2018) "Configuration model for correlation matrices preserving the node strength" algorithm for 
#'generating correlation matrices that have the same strength distribution as the original matrix.
#'
#' @param org_corr A \eqn{p \times p}  positive definite covariance (which will be transformed into a correlation matrix) or a correlation matrix
#' @param sampleSize The sample size for a correlation matrix, lower values for more sampling error
#' @param n Number of random correlations matrices to return
#' @param tol Tolerance for configuration model convergence. Default to .0001
#' @param stepSize Step size for configuration model NR solver. Increase to increase convergence speed. Default to .001
#' @param verbose Print convergence information to screen.

#' @return A list containing:  A \eqn{p \times p \times n} array of rewired correlation matrices and a \eqn{p \times p} matrix of the configuration model.
#' @export
#'
#' @examples
#' 
#' org_corr = matrix(0, 20, 20)
#' 
#' org_corr[1:10, 1:10] = .5
#' org_corr[11:20, 11:20] = .5
#' 
#' noise = t(sapply(as.list(1:100),FUN= function(x){return(rnorm(20,0,1))}))
#' 
#' org_corr = cov2cor(org_corr + cov(noise))
#' 
#' rand_corr = corrmat_rand(org_corr, 10000, 1)[[1]]
#' 
#' org_strength = colSums(org_corr)
#' rand_strength = colSums(rand_corr)
#' 
#' cor(org_strength, rand_strength)
corrmat_rand = function(org_corr, sampleSize, n, tol = .0001, stepSize = .001, verbose = F){
  
  org_corr = stats::cov2cor(org_corr)
  
  config_mod_list = config_mod(org_corr, tol, stepSize, verbose)
  
  config_mod_mat = stats::cov2cor(config_mod_list[[1]])
  
  toReturn = stats::rWishart(n,sampleSize,config_mod_mat)

  toReturn = apply(toReturn, FUN = stats::cov2cor, MARGIN = 3)
  toReturn = list(array(toReturn, c(dim(org_corr), n)), config_mod_mat)
  return(toReturn)
    
  
}

