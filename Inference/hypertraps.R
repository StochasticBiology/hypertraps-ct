library(Rcpp)
sourceCpp("hypertraps-r.cpp")
  
test.mat = matrix(c(0,0,1,
                    0,1,1,
                    1,1,1), nrow=3, ncol=3)
my.post = HyperTraPS(test.mat)

c4.dat = read.csv("../Data/c4-curated.csv", header=FALSE)
c4.mat = as.matrix(c4.dat)
c4.res = HyperTraPS(c4.mat, outputinput_arg = 1)
matrix(c4.res$best, nrow=len, ncol=len)
