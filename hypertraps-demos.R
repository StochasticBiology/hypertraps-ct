#################

# if we just want the HyperTraPS code without other dependencies, we only need Rcpp
# library(Rcpp)
# sourceCpp("hypertraps-r.cpp")

# if we want various helper functions and ggplot functions in addition, use this
source("hypertraps.R")

m.1 = matrix(rep(c(0,0,0,0,0,
               1,0,0,0,0,
               1,1,0,0,0,
               1,1,1,0,0,
               1,1,1,1,0,
               0,0,0,0,0,
               0,0,0,0,1,
               0,0,0,1,1,
               0,0,1,1,1,
               0,1,1,1,1),5), byrow=TRUE, ncol=5)
m.2 = matrix(rep(c(1,0,0,0,0,
               1,1,0,0,0,
               1,1,1,0,0,
               1,1,1,1,0,
               1,1,1,1,1,
               0,0,0,0,1,
               0,0,0,1,1,
               0,0,1,1,1,
               0,1,1,1,1,
               1,1,1,1,1),5), byrow=TRUE, ncol=5)
times = rep(c(0.1, 0.2, 0.3, 0.4, 0.5), 10)

### simple demo
my.post = HyperTraPS(m.2, initialstates_arg = m.1, starttimes_arg = times, featurenames_arg = c("A", "B", "C", "D", "E")); 
plotHypercube.summary(my.post)

# write output to files
writeHyperinf(my.post, "simpledemo", my.post$L, postlabel = "simpledemo", fulloutput=TRUE)

# retrieve output from files
my.post.r = readHyperinf("simpledemo", postlabel = "simpledemo", fulloutput=TRUE)
plotHypercube.summary(my.post.r)

# run an example with fewer walkers
my.post.sparse = HyperTraPS(m.2, initialstates_arg = m.1, starttimes_arg = times, featurenames_arg = c("A", "B", "C", "D", "E"), walkers_arg = 2); 

# q-gram distance
qgramdist(my.post, my.post.sparse)

### various other demos
# other plots
plotHypercube.motifs(my.post)
plotHypercube.timeseries(my.post)
plotHypercube.sampledgraph(my.post)

# regularisation
my.post.regularise = HyperTraPS(m.2, initialstates_arg = m.1, regularise_arg = 1, walkers_arg = 20)
plotHypercube.regularisation(my.post.regularise)
plotHypercube.summary(my.post.regularise)

# simulated annealing output
my.post.sa = HyperTraPS(m.2, initialstates_arg = m.1, sa_arg = 1)
plotHypercube.summary(my.post.sa)

# phenotypic landscape inference
my.post.pli = HyperTraPS(m.2, initialstates_arg = m.1, PLI_arg = 1)
plotHypercube.summary(my.post.pli)

# start with every edge parameterised, then regularise
my.post.bigmodel.regularise = HyperTraPS(m.2, initialstates_arg = m.1, model_arg = -1, regularise_arg = 1, walkers_arg = 20)
plotHypercube.regularisation(my.post.bigmodel.regularise)

# tool use paper reproduction
test.mat = as.matrix(read.table("RawData/total-observations.txt-trans.txt"))
my.names = as.vector(read.table("RawData/tools-names.txt"))[[1]]
starts = test.mat[seq(from=1, to=nrow(test.mat), by=2),]
ends = test.mat[seq(from=2, to=nrow(test.mat), by=2),]

my.post.tools = HyperTraPS(ends, initialstates_arg = starts, 
                           length_index_arg = 4, outputinput= 1, 
                           featurenames_arg = my.names) 
ggarrange(plotHypercube.lik.trace(my.post.tools), plotHypercube.bubbles(my.post.tools, reorder=TRUE), nrow=2)
plotHypercube.sampledgraph(my.post.tools, max=100)
