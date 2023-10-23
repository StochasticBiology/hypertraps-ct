#################
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
my.post.r = readHyperinf("simpledemo", my.post$L, postlabel = "simpledemo", fulloutput=TRUE)
plotHypercube.summary(my.post)

### various other demos
# regularisation
my.post.regularise = HyperTraPS(m.2, initialstates_arg = m.1, regularise_arg = 1, walkers_arg = 20)
plotHypercube.regularisation(my.post.regularise)

# simulated annealing output
my.post.sa = HyperTraPS(test.mat, sa_arg = 1)
plotHypercube.summary(my.post.sa)

# phenotypic landscape inference
my.post.pli = HyperTraPS(test.mat, PLI_arg = 1)
plotHypercube.summary(my.post.pli)

# start with every edge parameterised, then regularise
my.post.bigmodel.regularise = HyperTraPS(test.mat, model_arg = -1, regularise_arg = 1, walkers_arg = 20)
plotHypercube.regularisation(my.post.bigmodel.regularise)

ggarrange( ggplot(my.post))
plot(my.post$regularisation$reg.process$params, my.post$regularisation$reg.process$AIC)

# continuous time demo
test.mat = as.matrix(read.table("../VerifyData/synth-cross-samples-1.txt"))
test.start.times = rep(0.1,nrow(test.mat))
test.end.times = rep(0.2,nrow(test.mat))
my.post = HyperTraPS(test.mat, 
                     starttimes_arg = test.start.times, 
                     endtimes_arg = test.end.times, 
                     outputinput_arg = 1)
my.post.out = PosteriorAnalysis(my.post)
plotHypercube(my.post, my.post.out)

# tool use paper reproduction
test.mat = as.matrix(read.table("../Data/total-observations.txt-trans.txt"))
my.names = as.vector(read.table("../Data/tools-names.txt"))[[1]]
starts = test.mat[seq(from=1, to=nrow(test.mat), by=2),]
ends = test.mat[seq(from=2, to=nrow(test.mat), by=2),]

my.post.tools = HyperTraPS(ends, initialstates_arg = starts, 
                           length_index_arg = 4, outputinput= 1, 
                           featurenames_arg = my.names) 
plotHypercube(my.post.tools)

g.lik.trace = ggplot(my.post.tools$lik.traces) + geom_line(aes(x=sample.times, y=l.samples.1)) +
  geom_line(aes(x=sample.times, y=l.samples.2))
g.bubbles = ggplot(my.post.tools$Bubbles, aes(x=Time, y=factor(Name, levels=unique(my.post.tools$Bubbles$Name)), size=Probability)) +
  geom_point() 
ggarrange(g.lik.trace, g.bubbles, nrow=2)
