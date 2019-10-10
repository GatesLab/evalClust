#' @name perturbCorr
#' @aliases evalCorr
#' @title Perturb networks and evaluate subgroup structures of correlation matrices.
#' @description This function implements the Masuda, Kojaku & Sano (2018) "Configuration model for correlation matrices preserving the node strength" algorithm for 
#' generating correlation matrices that have the same strength distribution as the original matrix.  Returns modularity of each randomized correlation matrix.
#' @param sym.matrix A symmetric, correlation or covariance matrix object
#' @param sampleSize The sample size for a correlation matrix, lower values for more sampling error
#' @param n Number of random correlations matrices to return
#' @param tol Tolerance for configuration model convergence. Default to .0001
#' @param stepSize Step size for configuration model NR solver. Increase to increase convergence speed. Default to .001
#' @param verbose Logical, Print convergence information to screen.  Defaults to FALSE.
#' @param plot Logical, defaults to TRUE

#' @export 
#' @examples 
#' perturbCorr(examplecorr, sampleSize=100, n=25, plot=FALSE)


perturbCorr <- evalCorr <- function(sym.matrix, 
                                   plot = TRUE, 
                                   sampleSize = NULL,
                                   n = NULL,
                                   tol = .001,
                                   stepSize = .001,
                                   verbose = FALSE
                                   ){
  
  if (!isSymmetric(unname(sym.matrix))){ 
    # only recommended for count graphs; 
    # for correlation or those with (-) read in matrix
    #g <- sym.matrix 
    #sym.matrixg <- as.sym.matrix( get.adjacency(g, attr = "weight"))
    
    stop(paste("Symmetric matrix required."))
    
    
  } else{
    
    sym.matrixg                  <- as.data.frame(sym.matrix)
    g                            <- graph.adjacency(as.matrix(sym.matrixg), mode = "undirected", weighted = TRUE)
    
  }

    truemembership <- walktrap.community(g, weights = E(g)$weight, steps = 4)$membership

  # now randomly perturb and rewire
  
  n.elements       <- length(sym.matrix[,1])*(length(sym.matrix[,1])-1)/2 #number of unique elements; symmetric
  modularity.rando <- matrix(nrow = n)
  modularity.original <- modularity(walktrap.community(g, weights = E(g)$weight, steps = 4))
  
  randmats <- corrmat_rand(sym.matrix, sampleSize = sampleSize, n=n, tol=tol, stepSize=stepSize, verbose=verbose)
  randmats <- randmats[[1]]
  
  for(p in 1:n){ 
      new.g                 <- graph.adjacency(as.matrix(randmats[,,p]), mode = "undirected", weighted = TRUE)
      clust.sol             <- walktrap.community(new.g, weights = E(new.g)$weight)
      membership            <- clust.sol$membership
      modularity.rando[p] <-  modularity(clust.sol)
  }
  
 
  if (plot == TRUE){
    d <- density(modularity.rando)
    plot(d, main = "Comparison of original modularity to modularity of randomized graphs", xlab = "Modularity", ylab = "Density")
    abline(v=modularity.original, col='red')
}
  
  distribution <- sort(as.matrix(modularity.rando[,length(n)]))
  cutoff <- distribution[round(length(distribution)*.95)]
  
  res <- list(
    randmats = randmats,
    modularity = modularity.original,
    modularity.rando = modularity.rando
  )
  return(res)
}
