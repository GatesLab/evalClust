#' #' @name perturb
#' #' @title Perturb networks and evaluate subgroup structures.
#' #' @description Randomly rewires networks in increasing degrees of
#'    perturbation to evaluate stability of cluster solutions.
#' #' @param sym.matrix A symmetric, weighted sym.matrix or weighted graph
#' #'
#'
library(mcclust) # move this to somewhere else?
library(dils)
library(igraph)
library(ggplot2)

evalClust <- function(sym.matrix){

  if (is.igraph(sym.matrix))
  { g <- sym.matrix
  sym.matrix <- as.sym.matrix( get.adjacency(g, attr = "weight"))
  } else
    g                   <- graph.adjacency(as.matrix(sym.matrix), mode = "undirected", weighted = TRUE)

  diag(sym.matrix)        <- 0

  truemembership          <- walktrap.community(g, weights = E(g)$weight, steps = 4)$membership

  # now randomly perturb and rewire
  n.elements <- length(sym.matrix[,1])*(length(sym.matrix[,1])-1)/2
  percent <-seq(from=0, to = n.elements, by=max(round(0.01*(n.elements)), 1)) #disrupt 1% at a time
  VI<-matrix(,nrow = 100, ncol = length(percent))
  ARI<-matrix(,nrow = 100,ncol = length(percent))
  modularity.value <-matrix(,nrow = 100, ncol = length(percent))

  VI[,1] <- 0 #when 0 edges are disrupted no variation of information
  ARI[,1]<- 1 #when 0 edges are disrupted ARI = 1

  for(k in 1:length(percent)) {
    for(p in 1:100){ # for each degree of perturbation run 100 times
      new.v <- as.matrix(rewirematrix(sym.matrix, percent[k]))
      new.g                  <- graph.adjacency(new.v, mode = "undirected", weighted = TRUE)
      membership       <- walktrap.community(new.g, steps = 4)$membership
      modularity.value[p,k] <- modularity(walktrap.community(new.g, steps = 4))
      ## ARI and VI for comparison of this.g and membership in original
      VI[p,k] <- vi.dist(membership, truemembership)
      ARI[p,k] <-max(arandi(membership, truemembership), 0) #since ARI can be negative
    }
  }

  # create random graph, rewire, compare to the above in a figure.
  rando <- rewirematrix(sym.matrix, percent[length(percent)]) #create randomized network with same overall weight

  VI.rando<-matrix(,nrow = 100, ncol = length(percent))
  ARI.rando<-matrix(,nrow = 100,ncol = length(percent))
  modularity.rando <-matrix(,nrow = 100, ncol = length(percent))

  VI.rando[,1] <- 0 #when 0 edges are disrupted no variation of information
  ARI.rando[,1]<- 1 #when 0 edges are disrupted ARI = 1

  for(k in 2:length(percent)) {
    for(p in 1:100){ # for each degree of perturbation run 100 times
      new.v <-  as.matrix(rewirematrix(rando, percent[k]))
      new.g                  <- graph.adjacency(new.v, mode = "undirected", weighted = TRUE)
      membership       <- walktrap.community(new.g, steps = 4)$membership
      modularity.rando[p,k] <- modularity(walktrap.community(new.g, steps = 4))
      ## ARI and VI for comparison of this.g and membership in original
      VI.rando[p,k] <- vi.dist(membership, truemembership)
      ARI.rando[p,k] <-max(arandi(membership, truemembership), 0) #since ARI can be negative
    }
  }

  # plots of ARI and VI compared to original
  percentlab <- percent/n.elements

  # explicitly change 10 and 20 percent of community affiliations to add lines to graph
  comms <- unique(truemembership)
  lengths <- NULL
  for (p in 1:length(comms))
    lengths[p] <- length(which(truemembership == comms[p]))
  max <- which(lengths == max(lengths))
  perc10 <- round(.10*length(truemembership))
  perc20 <- round(.20*length(truemembership))
  tochange <- which(truemembership == comms[max])
  changed10 <- truemembership
   for (p in 1:perc10)
    changed10[tochange[p]] <- comms[sample(which(lengths == min(lengths)[1]), 1)]
  rep10ari <- matrix(arandi(changed10, truemembership), 1, length(percent))
  rep10vi <- matrix(vi.dist(changed10, truemembership), 1, length(percent))
  changed20 <- truemembership
  for (p in 1:perc20)
    changed20[tochange[p]] <- comms[sample(which(lengths == min(lengths)[1]), 1)]
  rep20ari <- matrix(arandi(changed20, truemembership), 1, length(percent))
  rep20vi <- matrix(vi.dist(changed20, truemembership), 1, length(percent))

  plotARI <- plot(percentlab, colMeans(ARI.rando), col = "red", main = "Comparison of original result against perturbed graphs: ARI", xlab = "Proportion Perturbed", ylab = "Mean ARI")
  plotARI <-  points(percentlab, colMeans(ARI), col = "black") + lines(percentlab, rep10ari) + lines(percentlab, rep20ari)

  plotVI <- plot(percentlab, colMeans(VI.rando), col = "red", main = "Comparison of original result against perturbed graphs: VI", xlab = "Proportion Perturbed", ylab = "Mean VI")
  plotVI <- plotVI + points(percentlab, colMeans(VI), col = "black") + lines(percentlab, rep10vi) + lines(percentlab, rep20vi)

  res <- list(
    VI = VI,
    ARI = ARI,
    modularity = modularity,
    VI.rando = VI.rando,
    ARI.rando = ARI.rando,
    modularity.rando = modularity.rando,
    percent = percentlab,
    plotVI = plotVI,
    plotARI = plotARI
  )
  return(res)
}
