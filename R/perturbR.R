#' @name perturbR
#' @aliases evalClust
#' @title Perturb networks and evaluate subgroup structures.
#' @description Randomly rewires networks in increasing degrees of
#'    perturbation to evaluate stability of community solutions obtained from Walktrap.
#' @param sym.matrix A symmetric, weighted matrix object
#' @param plot Logical, defaults to TRUE
#' @param resolution The percentage of edges to iteratively alter. One percent is default, increase to go quicker. 
#' @param reps The number of repititions to do for each level of perturbation. Decrease to make it go quicker. 
#' @export perturbR
#' @examples 
#' perturbR(exampledata, plot=FALSE, resolution=0.10, reps=1)


perturbR <- function(sym.matrix, plot = TRUE, resolution = 0.01, reps = 100){
  
  if (!isSymmetric(unname(sym.matrix))){ 
    # only recommended for count graphs; 
    # for correlation or those with (-) read in matrix
    #g <- sym.matrix 
    #sym.matrixg <- as.sym.matrix( get.adjacency(g, attr = "weight"))
    
    stop(paste("Symmetric matrix required."))
    
  } else{
    
    sym.matrixg                  <- as.data.frame(sym.matrix)
    diag(sym.matrixg)            <- 0
    sym.matrixg[sym.matrixg < 0] <- 0
    g                            <- graph.adjacency(as.matrix(sym.matrixg), mode = "undirected", weighted = TRUE)
    
  }
  
  truemembership <- walktrap.community(g, weights = E(g)$weight, steps = 4)$membership
  
  # now randomly perturb and rewire
  
  n.elements       <- length(sym.matrix[,1])*(length(sym.matrix[,1])-1)/2 #number of unique elements; symmetric
  percent          <- seq(from=0, to = n.elements, by=max(round(resolution*(n.elements)), 1)) #perturb 1% at a time
  VI               <- matrix(, nrow = reps, ncol = length(percent))
  ARI              <- matrix(, nrow = reps,ncol = length(percent))
  modularity.value <- matrix(,nrow = reps, ncol = length(percent))
  VI[,1]           <- 0 # when 0 edges are perturbed no variation of information
  ARI[,1]          <- 1 # when 0 edges are perturbed ARI = 1
  iters            <- seq(1,length(percent))

  for(p in 1:reps){ # for each degree of perturbation run 100 times
    
    for(k in 2:length(iters)) {
      
      new.v                 <- as.matrix(rewirematrix(sym.matrix, percent[iters[k]]))
      diag(new.v)           <- 0
      new.v[new.v< 0]       <- 0
      new.g                 <- graph.adjacency(as.matrix(new.v), mode = "undirected", weighted = TRUE)
      clust.sol             <- walktrap.community(new.g, weights = E(new.g)$weight)
      membership            <- clust.sol$membership
      modularity.value[p,k] <- modularity(clust.sol)
      VI[p,k]               <- vi.dist(membership, truemembership) # ARI and VI for comparison of this.g and membership in original
      ARI[p,k]              <- arandi(membership,  truemembership) #since ARI can be negative
      
    }
    
  }
  
  # explicitly change 10 and 20 percent of community affiliations to add lines to graph
  comms <- unique(truemembership)
  lengths <- NULL
  
  for (p in 1:length(comms)){
    
    lengths[p] <- length(which(truemembership == comms[p]))
    
  }
  
  max       <- which(lengths == max(lengths))
  perc10    <- round(.10*length(truemembership))
  perc20    <- round(.20*length(truemembership))
  tochange  <- which(truemembership == comms[max[1]])
  changed10 <- truemembership
  
  mincomm <- which(lengths == min(lengths))
  if (length(mincomm)>1)
    mincomm <- mincomm[2]
  
  for (p in 1:perc10){
    
    changed10[tochange[p]] <- comms[mincomm]
    
  }
  
  rep10ari  <- matrix(arandi(changed10, truemembership), 1, length(percent))
  rep10vi   <- matrix(vi.dist(changed10, truemembership), 1, length(percent))
  changed20 <- truemembership
  
  for (p in 1:perc20){
    
    changed20[tochange[p]] <- comms[mincomm]
    
  }
  
  # create random graph, rewire, compare to the above in a figure.
  rep20ari         <- matrix(arandi(changed20, truemembership), 1, length(percent))
  rep20vi          <- matrix(vi.dist(changed20, truemembership), 1, length(percent))
  rando            <- new.v
  VI.rando         <- matrix(,nrow = reps, ncol = length(percent))
  ARI.rando        <- matrix(,nrow = reps,ncol = length(percent))
  modularity.rando <- matrix(,nrow = reps, ncol = length(percent))
  plotVI           <- NULL
  plotARI          <- NULL
  percentlab       <- NULL
  
  if (plot == TRUE){

    VI.rando[,1] <- 0 #when 0 edges are disrupted no variation of information
    ARI.rando[,1]<- 1 #when 0 edges are disrupted ARI = 1
    
    for(k in 2 :length(percent)) {
      
      for(p in 1:reps){ 
        
        new.v <-  as.matrix(rewirematrix(rando, percent[k]))
        diag(new.v)      <- 0
        new.v[new.v< 0]  <- 0
        new.g            <- graph.adjacency(as.matrix(new.v), mode = "undirected", weighted = TRUE)
        clust.sol        <-  walktrap.community(new.g, weights = E(new.g)$weight, steps = 4)
        membership       <- clust.sol$membership
        modularity.rando[p,k] <- modularity(clust.sol)
        VI.rando[p,k]         <- vi.dist(membership, truemembership)        ## ARI and VI for comparison of this.g and membership in original
        ARI.rando[p,k]        <- max(arandi(membership, truemembership), 0) ## since ARI can be negative
        
      }
    }
    
    # plots of ARI and VI compared to original
    percentlab <- percent/n.elements
  
    plotARI <- plot(percentlab, colMeans(ARI), col = "black", main = "Comparison of original result against perturbed graphs: ARI", xlab = "Proportion Perturbed", ylab = "Mean ARI")
    plotARI <-  points(percentlab, colMeans(ARI.rando), col = "red") + lines(percentlab, rep10ari) + lines(percentlab, rep20ari)
    
    plotVI <- plot(percentlab, colMeans(VI.rando), col = "red", main = "Comparison of original result against perturbed graphs: VI", xlab = "Proportion Perturbed", ylab = "Mean VI")
    plotVI <- plotVI + points(percentlab, colMeans(VI), col = "black") + lines(percentlab, rep10vi) + lines(percentlab, rep20vi)
  }
  
  res <- list(
    VI  = VI, # only one column if Plot == FALSE; column is the 20% perturb point
    ARI = ARI,
    modularity = modularity.value,
    VI.rando = VI.rando, 
    ARI.rando = ARI.rando,
    modularity.rando = modularity.rando,
    percent = percentlab,
    ari10mark = rep10ari[1],
    vi10mark = rep10vi[1],
    ari20mark = rep20ari[1],
    vi20mark = rep20vi[1],
    plotVI = plotVI,
    plotARI = plotARI
  )
  return(res)
}
