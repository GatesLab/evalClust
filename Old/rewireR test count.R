## this function follows from the work describe in "The weighted random graph model" by Garlaschelli,
## New Journal of Physics, 11, (2009), 073005.
# only undirected, weighted graphs are considered here.
# nperturb is the number of edges to perturb.
# matrix is the symmetric graph to rewire.
rewirematrixC <- function(sym.matrix, nperturb, type)
{
  eligable <-
    which(lower.tri(sym.matrix), arr.ind = T)#ensure that diagonal isn't considered, and all potential edges are
  toalter <- sample(1:length(eligable[, 1]), nperturb)
  new.v <- sym.matrix
  if (type == "count") {
    # W <-
    #   sum(colSums(sym.matrix)) / 2 # sum of strengths; divide by two since undirected
    # # Optimal parameter choice - retain overall weight, not edge number; Garlaschelli 2009
    # pstar <- (2 * W) / (length(sym.matrix[, 1]) * (length(sym.matrix[, 1]) - 1) +
    #                       2 * W)
    for (l in 1:length(toalter)) {
      # new.v[eligable[toalter[l],1],eligable[toalter[l],2]]<- 0
      W <-
        (sum(sym.matrix[eligable[toalter[l],1],]) + sum(sym.matrix[,eligable[toalter[l],2]])) / 2 # sum of strengths; divide by two since undirected
             # # Optimal parameter choice - retain overall weight, not edge number; Garlaschelli 2009
             pstar <- (2 * W) / (2*(length(sym.matrix[,1])) + 2 * W)
             bernoulli <- 1
             toadd <- 0
             while (bernoulli == 1) {
               bernoulli <- rbinom(1, 1, pstar)
               toadd <- bernoulli + toadd
             }
             new.v[eligable[toalter[l], 1], eligable[toalter[l], 2]] <-
               toadd
    }
  } else if (type == "cov") { #need to update - P.J. Brown 1994 has covariance analogue 
    if (length(toalter)>1)
    {
      randomized <- sample(toalter)
    } else {
      randomized <- toalter
    }
    for (l in 1:length(randomized))
      new.v[eligable[toalter[l],1],eligable[toalter[l],2]]<- new.v[[eligable[randomized[l],1],eligable[randomized[l],2]]] 
  } else if (type == "cor") {
    for (l in 1:length(toalter))
      new.v[eligable[toalter[l],1], eligable[toalter[l],2]] <- rjm(sym.matrix, 1, row = eligable[toalter[l],1], col = eligable[toalter[l],2])
  }
  for (l in 1:length(toalter))
    new.v[eligable[toalter[l],2], eligable[toalter[l],1]] <- new.v[eligable[toalter[l],1], eligable[toalter[l],2]] 
  return(as.data.frame(new.v))
}
