library(Rcpp)
sourceCpp("hypertraps-r.cpp")
  
my.post = HyperTraPS( c(0,0,0,0,0,1,
                            0,0,0,0,1,1,
                            0,0,0,1,1,1), 3, 3, kernel_index_arg = 4, outputinput_arg = 1)
matrix(my.post, nrow=3, ncol=3)

c4.dat = read.csv("../Data/c4-curated.csv", header=FALSE)
c4.mat = cbind(0*as.matrix(c4.dat), as.matrix(c4.dat))
c4.vec = c(apply(c4.mat, 1, c))
len = ncol(c4.dat)
c4.res = HyperTraPS(c4.vec, len_arg = len, ntarg_arg = nrow(c4.dat),
                    crosssectional_arg = 1)
matrix(c4.res, nrow=len, ncol=len)
