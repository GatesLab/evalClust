## this function follows from the work describe in "The weighted random graph model" by Garlaschelli,
## New Journal of Physics, 11, (2009), 073005.
# only undirected, weighted graphs are considered here.
# nperturb is the number of edges to perturb.
# matrix is the symmetric graph to rewire.

rewirematrix <- function(sym.matrix, nperturb)
{
  eligable <-
    which(lower.tri(sym.matrix), arr.ind = T)#ensure that diagonal isn't considered, and all potential edges are
  toalter <- sample(1:length(eligable[, 1]), nperturb)
  new.v <- sym.matrix
#  if (type == "count") {
    W <-
      sum(colSums(sym.matrix)) / 2 # sum of strengths; divide by two since undirected
    # Optimal parameter choice - retain overall weight, not edge number; Garlaschelli 2009
    pstar <- (2 * W) / (length(sym.matrix[, 1]) * (length(sym.matrix[, 1]) - 1) +
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
  for (l in 1:length(toalter)) #make symmetric again
    new.v[eligable[toalter[l],2], eligable[toalter[l],1]] <- new.v[eligable[toalter[l],1], eligable[toalter[l],2]] 
  return(as.data.frame(new.v))
}
