#' @name perturbR
#' @aliases evalClust
#' @title Perturb networks and evaluate subgroup structures.
#' @description Randomly rewires networks in increasing degrees of
#'    perturbation to evaluate stability of community solutions obtained from Walktrap.
#' @param sym.matrix A symmetric, sparse count matrix object
#' @param plot Logical, defaults to TRUE
#' @param resolution The percentage of edges to iteratively alter. One percent is default, increase to go quicker. 
#' @param reps The number of repititions to do for each level of perturbation. Decrease to make it go quicker. 
#' @param confirm_cluster Dataframe. Option to provide confirmatory cluster labels contained in the dataframe. Dataframe has 2 columns,
#' the first referring to node number, and the second an integer variable referring to cluster label/assignment.
#' @param errbars Logical, defaults to FALSE. Option to add error bars of one standard deviation above and below the mean for each point.
#' @export perturbR
#' @examples 
#' perturbR(exampledata, plot=FALSE, resolution=0.10, reps=1, cluster_assign = NULL, errbars = FALSE)


perturbR <- evalClust <- function(sym.matrix, plot = TRUE, resolution = 0.01, reps = 100, cluster_assign = NULL, errbars = FALSE){
  
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
  if(is.null(cluster_assign)){
  truemembership <- walktrap.community(g, weights = E(g)$weight, steps = 4)$membership
  
  } else {
    truemembership <- cluster_assign[order(cluster_assign[,1]),]
    truemembership <- truemembership[,2]
  }
  # now randomly perturb and rewire
  
  n.elements       <- length(sym.matrix[,1])*(length(sym.matrix[,1])-1)/2 #number of unique elements; symmetric
  percent          <- seq(from=0, to = n.elements, by=max(round(resolution*(n.elements)), 1)) #perturb 1% at a time
  VI               <- matrix(, nrow = reps, ncol = length(percent))
  ARI              <- matrix(, nrow = reps,ncol = length(percent))
  modularity.value <- matrix(,nrow = reps, ncol = length(percent))
  modularity.value[,1] <- modularity(walktrap.community(g, weights = E(g)$weight, steps = 4))
  VI[,1]           <- 0 # when 0 edges are perturbed no variation of information
  ARI[,1]          <- 1 # when 0 edges are perturbed ARI = 1
  iters            <- seq(1,length(percent))

  for(p in 1:reps){ # for each degree of perturbation run 100 times
    
    for(k in 2:length(iters)) {
      
      new.v                 <- as.matrix(rewireR(sym.matrix, percent[iters[k]]))
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
  
  for (p in 1:length(comms))
    lengths[p] <- length(which(truemembership == comms[p]))

  
#  max       <- which(lengths == max(lengths))
  comms <- unique(truemembership)
  perc10    <- round(.10*length(truemembership))
  perc20    <- round(.20*length(truemembership))
#  tochange  <- which(truemembership == comms[max[1]])
  tochange10 <- sample(seq(1,length(sym.matrix[,1])), perc10)
  tochange20 <- sample(seq(1,length(sym.matrix[,1])), perc20)
  
  changed10 <- truemembership

  for (p in 1:perc10){
    commchange <- comms[-truemembership[tochange10[p]]]
    changed10[tochange10[p]] <- commchange[sample(length(commchange), 1)] #randomly select which community it gets assigned
    }
  
  rep10ari  <- matrix(arandi(changed10, truemembership), 1, length(percent))
  rep10vi   <- matrix(vi.dist(changed10, truemembership), 1, length(percent))
  changed20 <- truemembership
  
  for (p in 1:perc20){
    commchange <- comms[-truemembership[tochange20[p]]]
    changed20[tochange20[p]] <- commchange[sample(length(commchange), 1)] #randomly select which community it gets assigned
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
  
  diag(new.v)      <- 0
  new.v[new.v< 0]  <- 0
  new.g            <- graph.adjacency(as.matrix(new.v), mode = "undirected", weighted = TRUE)
  modularity.rando[,1]       <-  modularity(walktrap.community(new.g, weights = E(new.g)$weight, steps = 4))
  
  if (plot == TRUE){

    VI.rando[,1] <- 0 #when 0 edges are disrupted no variation of information
    ARI.rando[,1]<- 1 #when 0 edges are disrupted ARI = 1
    
    for(k in 2 :length(percent)) {
      
      for(p in 1:reps){ 
        
        new.v <-  as.matrix(rewireR(rando, percent[k]))
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
  
    plotARI <- plot(percentlab, colMeans(ARI),  pch=19, col = "black", main = "Comparison of original result against perturbed graphs: ARI", xlab = "Proportion Perturbed", ylab = "Mean ARI")
    if(errbars == TRUE)
    arrows(percentlab, (colMeans(ARI)-apply(ARI, 2, sd)), percentlab, (colMeans(ARI)+apply(ARI, 2, sd)), length=0.05, angle=90, code=3)
    plotARI <-  points(percentlab, colMeans(ARI.rando), pch=19, col = "red") + lines(percentlab, rep10ari) + lines(percentlab, rep20ari)
    if(errbars == TRUE)
      arrows(percentlab, (colMeans(ARI.rando)-apply(ARI.rando, 2, sd)), percentlab, (colMeans(ARI.rando)+apply(ARI.rando, 2, sd)), length=0.05, angle=90, code=3)
    
    plotVI <- plot(percentlab, colMeans(VI.rando), ylim=range(c(0, colMeans(VI.rando)+apply(VI.rando, 2, sd))), pch=19,
                   col = "red", main = "Comparison of original result against perturbed graphs: VI", xlab = "Proportion Perturbed", ylab = "Mean VI")
    if(errbars == TRUE)
      arrows(percentlab, (colMeans(VI.rando)-apply(VI.rando, 2, sd)), percentlab, (colMeans(VI.rando)+apply(VI.rando, 2, sd)), length=0.05, angle=90, code=3)
    plotVI <- plotVI + points(percentlab, colMeans(VI), pch=19, col = "black") + lines(percentlab, rep10vi) + lines(percentlab, rep20vi)
    if(errbars == TRUE)
      arrows(percentlab, (colMeans(VI)-apply(VI, 2, sd)), percentlab, (colMeans(VI)+apply(VI, 2, sd)), length=0.05, angle=90, code=3)
    
    }
  
  distribution <- sort(as.matrix(modularity.rando[,length(percent)]))
  cutoff <- distribution[round(length(distribution)*.95)]
  
  res <- list(
    VI  = VI, # only one column if Plot == FALSE; column is the 20% perturb point
    ARI = ARI,
    VI.rando = VI.rando, 
    ARI.rando = ARI.rando,
    modularity = modularity.value,
    modularity.rando = modularity.rando,
    percent = percentlab,
    ari10mark = rep10ari[1],
    vi10mark = rep10vi[1],
    ari20mark = rep20ari[1],
    vi20mark = rep20vi[1],
    plotVI = plotVI,
    plotARI = plotARI,
    cutoff = cutoff
  )
  return(res)
}
