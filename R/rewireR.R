## this function follows from the work describe in "The weighted random graph model" by Garlaschelli,
## New Journal of Physics, 11, (2009), 073005.
# only undirected, weighted graphs are considered here.
# nperturb is the number of edges to perturb.
# matrix is the symmetric graph to rewire.
rewirematrix <- function(matrix, nperturb, type)
{
  diag(matrix) <- 0
  matrix[which(matrix < 0)] = 0 # set negative values to zero; not considered in most comm. det. approaches
  eligable <-
    which(lower.tri(matrix), arr.ind = T)#ensure that diagonal isn't considered, and all potential edges are
  toalter <- sample(1:length(eligable[, 1]), nperturb)
  new.v <- matrix
  if (type == "count") {
    W <-
      sum(colSums(matrix)) / 2 # sum of strengths; divide by two since undirected
    # Optimal parameter choice - retain overall weight, not edge number; Garlaschelli 2009
    pstar <- (2 * W) / (length(matrix[, 1]) * (length(matrix[, 1]) - 1) +
                          2 * W)
    for (l in 1:length(toalter)) {
      # new.v[eligable[toalter[l],1],eligable[toalter[l],2]]<- 0
      bernoulli <- 1
      toadd <- 0
      while (bernoulli == 1) {
        bernoulli <- rbinom(1, 1, pstar)
        toadd <- bernoulli + toadd
      }
      new.v[eligable[toalter[l], 1], eligable[toalter[l], 2]] <-
        toadd
      }
    elseif (type == "cov")
    {
      if (length(toalter)>1)
      {
        randomized <- sample(toalter)
      } else {
        randomized <- toalter
      }
      for (l in 1:length(randomized))
        new.v[eligable[toalter[l],1],eligable[toalter[l],2]]<- new.v[[eligable[randomized[l],1],eligable[randomized[l],2]]] 
    }
  }
  new.v[upper.tri(new.v)] <-
    t(new.v)[upper.tri(new.v)] # maintain symmetry
  return(as.data.frame(new.v))
}
