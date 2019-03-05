#' @name perturbR
#' @aliases evalClust
#' @title Perturb networks and evaluate subgroup structures.
#' @description Randomly rewires networks in increasing degrees of
#'    perturbation to evaluate stability of community solutions obtained from Walktrap.
#' @param sym.matrix A symmetric, sparse count matrix object
#' @param plot Logical, defaults to TRUE
#' @param resolution The percentage of edges to iteratively alter. One percent is default, increase to go quicker. 
#' @param reps The number of repititions to do for each level of perturbation. Decrease to make it go quicker. 
#' @param errbars Logical, defaults to FALSE. Option to add error bars of one standard deviation above and below the mean for each point.
#' @export perturbR 
#' @examples 
#' perturbR(exampledata, plot=FALSE, resolution=0.10, reps=1, errbars = FALSE)


perturbR <- evalClust <- function( sym.matrix, 
                                   plot = TRUE, 
                                   resolution = 0.01, 
                                   reps = 100,
                                   errbars = FALSE ){
  
  cluster_assign = NULL # removed option for now to have a priori clusters
  dist = "NegBinom"     # removed argument for now to have neg binom or diff distribution
  
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
      
      new.v                 <- as.matrix(rewireR(sym.matrix, nperturb = percent[iters[k]], dist = dist))
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
  
  # randomly order vector to see modularity distribution
  modularity.randclust <- matrix(NA, reps, 1)
  
  for (p in 1:reps)
    modularity.randclust[p] <-modularity(g, sample(truemembership))
  
  # explicitly change 10 and 20 percent of community affiliations to add lines to graph
  comms <- unique(truemembership)
  lengths <- NULL
  
  for (p in 1:length(comms))
    lengths[p] <- length(which(truemembership == comms[p]))
  
  #  max       <- which(lengths == max(lengths))
  perc10    <- round(.10*length(truemembership))
  perc20    <- round(.20*length(truemembership))
  
  rep10arim <- matrix(,100, 1)
  rep10vim <- matrix(,100, 1)
  for (k in 1:100){
    changed10 <- truemembership
    #  tochange  <- which(truemembership == comms[max[1]])
    tochange10 <- sample(seq(1,length(sym.matrix[,1])), perc10)
    
    for (p in 1:perc10){
      commchange <- comms[-which(comms == truemembership[tochange10[p]])]
      changed10[tochange10[p]] <- commchange[sample(length(commchange), 1)] #randomly select which community it gets assigned
    }
    
    rep10arim[k]  <- arandi(changed10, truemembership)
    rep10vim[k]   <- vi.dist(changed10, truemembership)
  }
  rep10ari  <- matrix(mean(rep10arim), 1, length(percent))
  rep10vi   <- matrix(mean(rep10vim), 1, length(percent))
  
  rep20arim <- matrix(,100, 1)
  rep20vim <- matrix(,100, 1)
  for (k in 1:100){
    changed20 <- truemembership
    tochange20 <- sample(seq(1,length(sym.matrix[,1])), perc20)
    
    for (p in 1:perc20){
      commchange <- comms[-which(comms == truemembership[tochange20[p]])]
      changed20[tochange20[p]] <- commchange[sample(length(commchange), 1)] #randomly select which community it gets assigned
    }
    
    rep20arim[k]  <- arandi(changed20, truemembership)
    rep20vim[k]   <- vi.dist(changed20, truemembership)
  }
  rep20ari  <- matrix(mean(rep20arim), 1, length(percent))
  rep20vi   <- matrix(mean(rep20vim), 1, length(percent))
  
  rando            <- new.v
  VI.rando         <- matrix(,nrow = reps, ncol = length(percent))
  ARI.rando        <- matrix(,nrow = reps,ncol = length(percent))
  modularity.rando <- matrix(,nrow = reps, ncol = length(percent))
  plotVI           <- NULL
  plotARI          <- NULL

  diag(new.v)      <- 0
  new.v[new.v< 0]  <- 0
  new.g            <- graph.adjacency(as.matrix(new.v), mode = "undirected", weighted = TRUE)
  modularity.rando[,1]       <-  modularity(walktrap.community(new.g, weights = E(new.g)$weight, steps = 4))
  
  percentlab <- percent/n.elements
  
  if (plot == TRUE){
    
    VI.rando[,1] <- 0 #when 0 edges are disrupted no variation of information
    ARI.rando[,1]<- 1 #when 0 edges are disrupted ARI = 1
    
    for(k in 2 :length(percent)) {
      
      for(p in 1:reps){ 
        
        new.v <-  as.matrix(rewireR(rando, nperturb = percent[k], dist = dist))
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
    
    plotARI <- plot(percentlab, colMeans(ARI),  pch=16, col = "black", main = "Comparison of original result against perturbed graphs: ARI", xlab = "Proportion Perturbed", ylab = "Mean ARI")
    if(errbars == TRUE)
      graphics::arrows(percentlab, (colMeans(ARI)-apply(ARI, 2, stats::sd)), percentlab, (colMeans(ARI)+apply(ARI, 2, stats::sd)), length=0.05, angle=90, code=3)
    plotARI <-  points(percentlab, colMeans(ARI.rando), pch=18, col = "red") + lines(percentlab, rep10ari) + lines(percentlab, rep20ari)
    if(errbars == TRUE)
      graphics::arrows(percentlab, (colMeans(ARI.rando)-apply(ARI.rando, 2, stats::sd)), percentlab, (colMeans(ARI.rando)+apply(ARI.rando, 2, stats::sd)), length=0.05, angle=90, code=3)
    
    plotVI <- plot(percentlab, colMeans(VI.rando), ylim=range(c(0, colMeans(VI.rando)+apply(VI.rando, 2, stats::sd))), pch=18,
                   col = "red", main = "Comparison of original result against perturbed graphs: VI", xlab = "Proportion Perturbed", ylab = "Mean VI")
    if(errbars == TRUE)
      graphics::arrows(percentlab, (colMeans(VI.rando)-apply(VI.rando, 2, stats::sd)), percentlab, (colMeans(VI.rando)+apply(VI.rando, 2, stats::sd)), length=0.05, angle=90, code=3)
    plotVI <- plotVI + points(percentlab, colMeans(VI), pch=16, col = "black") + lines(percentlab, rep10vi) + lines(percentlab, rep20vi)
    if(errbars == TRUE)
      graphics::arrows(percentlab, (colMeans(VI)-apply(VI, 2, stats::sd)), percentlab, (colMeans(VI)+apply(VI, 2, stats::sd)), length=0.05, angle=90, code=3)
    
  }
  
  distribution <- sort(as.matrix(modularity.rando[,length(percent)]))
  cutoff <- distribution[round(length(distribution)*.95)]
  
  res <- list(
    comm.assign = truemembership,
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
