# Reference: 
# Joe H (2006). Generating random correlation matrices based on partial
# correlations. J. Mult. Anal. Vol. 97, 2177--2189.

# Generate random correlation matrix based on random partial correlations
# choice of alpha parameter lead to invariance to index order.
# d = dimension, 
# alphad = alpha parameter for partial of 1,d given 2,...,d-1
# default value alphad = 1 leads to random matrix uniform over
#      space of positive definite correlation matrices
# Other each correlation has a Beta(a,a) distribution on (-1,1) where
#   a=alphad+(d-2)/2

# rsub is the symmetrix matrix, 
# generate the correlation for the (row,col) component
#
rjm2 <-function(rsub, row, col)
{ b<-nrow(rsub)     # b is dimension 
alp <-  1 
ii<-1:b
ii <- ii[-c(row, col)]
r1<-rsub[ii,row]
r3<-rsub[ii,col]
# ii <- 2:(b-1)
# r1 <- rsub[ii,1]
# r3 <- rsub[ii,b]
R2<-rsub[ii,ii]
Ri<-solve(R2)       #partial correlation matrix
rcond<-2*rbeta(1,alp,alp)-1      #draw from beta parameter
tem13<-t(r1)%*%Ri%*%r3  
tem11<-t(r1)%*%Ri%*%r1
tem33<-t(r3)%*%Ri%*%r3
res<-tem13+rcond*sqrt((1-tem11)*(1-tem33))
return(res)
}