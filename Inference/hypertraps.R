library(Rcpp)
sourceCpp("hypertraps-r.cpp")
HyperTraPS( c(0,0,0,0,0,1,
              0,0,0,0,1,1,
	      0,0,0,1,1,1), 3, 3)
	      