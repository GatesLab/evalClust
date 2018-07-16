#' @name rewireR
#' @title Rewire graph by randomly assigning new values for a given degree of perturbation.
#' @description Randomly rewires graphs by altering a specific number of edges using Bernoulli 
#' trials as described in "The weighted random graph model" by Garlaschelli,
#' New Journal of Physics, 11, (2009), 073005. Only undirected, weigghted count matrices 
#' are considered here.
#' @param sym.matrix  A symmetric, sparse count matrix object.
#' @param nperturb The number of edges to randomly alter.
#' @param dist Option to rewire in a manner that retains overall graph weight regardless of distribution of edge weights. 
#' This option is invoked by putting any text into this field. Defaults to "NegBinom" for negative binomial.
#' @export rewireR
#' @examples 
#' rewireR(exampledata, nperturb=40, dist = "Normal")


rewireR <- function(sym.matrix, nperturb, dist)
{
  eligable <-
    which(lower.tri(sym.matrix), arr.ind = T)#ensure that diagonal isn't considered, and all potential edges are
  toalter <- sample(1:length(eligable[, 1]), nperturb)
  new.v <- sym.matrix
  maxEdge <- max(sym.matrix)
  #  if (type == "count") {
  if (dist == "NegBinom"){
    W <-
      sum(sym.matrix[lower.tri(sym.matrix)]) # sum of strengths; 
    # Optimal parameter choice - retain overall graph weight, not edge number; Garlaschelli 2009
    pstar <- (2*W) / (length(sym.matrix[, 1]) * (length(sym.matrix[, 1]) - 1) +
                         2*W)
    for (l in 1:length(toalter)) {
      # new.v[eligable[toalter[l],1],eligable[toalter[l],2]]<- 0
      stopcri <- maxEdge + 1
      while (stopcri > maxEdge){
      bernoulli <- 1
      toadd <- 0
      while (bernoulli == 1) {
        bernoulli <- rbinom(1, 1, pstar)
        toadd <- bernoulli + toadd
      }
      stopcri <- toadd
      }
      new.v[eligable[toalter[l], 1], eligable[toalter[l], 2]] <-
        toadd
    } }else {
        if (length(toalter)>1){
        randomized <- sample(toalter)
        } else {
          randomized <- toalter
        }
      for (l in 1:length(randomized))
        new.v[eligable[toalter[l],1],eligable[toalter[l],2]]<- new.v[[eligable[randomized[l],1],eligable[randomized[l],2]]] 
      }
  for (l in 1:length(toalter)) #make symmetric again
    new.v[eligable[toalter[l],2], eligable[toalter[l],1]] <- new.v[eligable[toalter[l],1], eligable[toalter[l],2]] 
  return(as.data.frame(new.v))
  
}
