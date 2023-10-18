library(Rcpp)
library(ggplot2)

sourceCpp("hypertraps-r.cpp")

test.mat = as.matrix(read.table("../VerifyData/synth-cross-samples-1.txt"))
my.post = HyperTraPS(test.mat)

my.post.out = PosteriorAnalysis(my.post)

ggplot(my.post.out$Bubbles, aes(x=Time, y=OriginalIndex, size=Probability)) + 
  geom_point()


test.mat = matrix(c(0,0,1,
                    0,1,1,
                    1,1,1), nrow=3, ncol=3)
my.post = HyperTraPS(test.mat)

c4.dat = read.csv("../Data/c4-curated.csv", header=FALSE)
c4.mat = as.matrix(c4.dat)
c4.res = HyperTraPS(c4.mat, outputinput_arg = 1)
matrix(c4.res$best, nrow=len, ncol=len)
