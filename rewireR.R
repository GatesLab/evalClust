## this function follows from the work describe in "The weighted random graph model" by Garlaschelli,
## New Journal of Physics, 11, (2009), 073005.
# only undirected, weighted graphs are considered here.
# nperturb is the number of edges to perturb. 
# matrix is the symmetric graph to rewire. 
rewirematrix <- function(matrix, nperturb)
{
  diag(matrix) <- 0
  matrix[which(matrix < 0)] =0 # set negative values to zero; not considered in most comm. det. approaches
  cov <- NULL
  if (max(matrix) < 1 & max(matrix)>0){  #is it a correlation matrix? 
    matrix = round(matrix *100) # scale up
    cov <- 1}
  W <- sum(colSums(matrix))/2 # sum of strengths; divide by two since undirected
  # Optimal parameter choice - retain overall weight, not edge number; Garlaschelli 2009
  pstar <- (2*W)/(length(matrix[,1])*(length(matrix[,1])-1)+2*W)
  eligable <- which(lower.tri(matrix),arr.ind = T)#ensure that diagonal isn't considered, and all potential edges are
  toalter <- sample(1:length(eligable[,1]), nperturb)
  new.v <- matrix
  for (l in 1:length(toalter)){
    new.v[eligable[toalter[l],1],eligable[toalter[l],2]]<- 0
    bernoulli <- 1
    while (bernoulli == 1){
      bernoulli <- rbinom(1,1,pstar)
      new.v[eligable[toalter[l],1],eligable[toalter[l],2]]<- new.v[eligable[toalter[l],1],eligable[toalter[l],2]] + bernoulli
    }
    # maintain symmetry:
    new.v[eligable[toalter[l],2],eligable[toalter[l],1]]<- new.v[eligable[toalter[l],1],eligable[toalter[l],2]]
  }
  if (!is.null(cov))
    new.v <- new.v/100
  return(as.data.frame(new.v))
}
