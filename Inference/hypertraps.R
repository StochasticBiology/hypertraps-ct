library(Rcpp)
sourceCpp("hypertraps-r.cpp")
my.post = HyperTraPS( rep(c(0,0,0,0,0,1,
              0,0,0,0,1,1,
	            0,0,0,1,1,1),10), 3, 30, kernel_index_arg = 4, outputinput_arg = 1)
	      
my.post = HyperTraPS( c(0,0,0,0,0,1,
                            0,0,0,0,1,1,
                            0,0,0,1,1,1), 3, 3, kernel_index_arg = 4, outputinput_arg = 1)
matrix(my.post, nrow=3, ncol=3)
