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
my.post = HyperTraPS(m.2, initialstates_arg = m.1, 
                     starttimes_arg = times, endtimes_arg = times, 
                     length_index_arg = 4,
                     featurenames_arg = c("A", "B", "C", "D", "E")); 
plotHypercube.summary(my.post)
plotHypercube.sampledgraph2(my.post, thresh=0.1, use.arc=FALSE, edge.label.size=3) + 
  theme(legend.position="none") + expand_limits(x = c(-0.1, 1.1))

# write output to files
writeHyperinf(my.post, "simpledemo", my.post$L, postlabel = "simpledemo", fulloutput=TRUE)

# retrieve output from files
my.post.r = readHyperinf("simpledemo", postlabel = "simpledemo", fulloutput=TRUE)
plotHypercube.summary(my.post.r)

# run an example with fewer walkers
my.post.sparse = HyperTraPS(m.2, initialstates_arg = m.1, 
                            starttimes_arg = times, endtimes_arg = times,
                            featurenames_arg = c("A", "B", "C", "D", "E"), walkers_arg = 2); 
plotHypercube.summary(my.post.sparse)

# direct time run (no time window specified)
my.post.dt = HyperTraPS(m.2, initialstates_arg = m.1, featurenames_arg = c("A", "B", "C", "D", "E")); 
plotHypercube.summary(my.post.dt)

### various other demos
# other plots
plotHypercube.motifs(my.post)
plotHypercube.timeseries(my.post)
plotHypercube.graph(my.post)

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

#### short-form examples from past studies -- should run in a few minutes and give approximations to the original results

# ovarian cancer case study reproduction
# traits are chromosomal aberrations, observations are independent patient samples
cgh.mat = readLines("RawData/ovarian.txt")
cgh.mat = do.call(rbind, lapply(strsplit(cgh.mat, ""), as.numeric))
cgh.names = as.vector(read.table("RawData/ovarian-names.txt", sep=","))[[1]]

my.post.cgh = HyperTraPS(cgh.mat, 
                        length_index_arg = 4, outputinput= 1, 
                        featurenames_arg = cgh.names) 
ggarrange(plotHypercube.lik.trace(my.post.cgh), 
          plotHypercube.bubbles(my.post.cgh, reorder=TRUE), 
          plotHypercube.sampledgraph2(my.post.cgh, no.times=TRUE), nrow=3)
plotHypercube.sampledgraph2(my.post.cgh, no.times=TRUE)

# C4 paper reproduction
# traits are physical/genetic features associated with C4, observations are (incomplete) phylogenetically independent intermediate species
c4.mat = as.matrix(read.table("RawData/c4-curated.csv", sep=","))
c4.names = as.vector(read.table("RawData/c4-trait-names.txt", sep=","))[[1]]

my.post.c4 = HyperTraPS(c4.mat, 
                        length_index_arg = 4, outputinput= 1, 
                        losses_arg = 1,
                        featurenames_arg = c4.names) 
ggarrange(plotHypercube.lik.trace(my.post.c4), plotHypercube.bubbles(my.post.c4, reorder=TRUE), nrow=2)

# malaria paper reproduction
# traits are clinical features, observations are (incomplete) independent patient presentations
malaria.df = read.csv("RawData/jallow_dataset_binary_with2s.csv")
malaria.mat = as.matrix(malaria.df[,2:ncol(malaria.df)])
malaria.names = as.vector(read.table("RawData/malaria-names.txt", sep=","))[[1]]

my.post.malaria = HyperTraPS(malaria.mat, 
                        length_index_arg = 3, outputinput= 1, 
                        kernel_index_arg = 2,
                        walkers_arg = 20,
                        featurenames_arg = malaria.names) 
ggarrange(plotHypercube.lik.trace(my.post.malaria), plotHypercube.bubbles(my.post.malaria, reorder=TRUE, transpose=TRUE), nrow=2)

# tool use paper reproduction
# traits are modes of tool use, observations are phylogenetically coupled species observations (phylogeny has been accounted for, giving transition pairs)
tools.mat = as.matrix(read.table("RawData/total-observations.txt-trans.txt"))
tools.names = as.vector(read.table("RawData/tools-names.txt"))[[1]]
tools.starts = tools.mat[seq(from=1, to=nrow(tools.mat), by=2),]
tools.ends = tools.mat[seq(from=2, to=nrow(tools.mat), by=2),]

my.post.tools = HyperTraPS(tools.ends, initialstates_arg = tools.starts, 
                           length_index_arg = 4, outputinput= 1, 
                           featurenames_arg = tools.names) 
ggarrange(plotHypercube.lik.trace(my.post.tools), plotHypercube.bubbles(my.post.tools, reorder=TRUE, transpose=TRUE), nrow=2)
plotHypercube.sampledgraph(my.post.tools, max=100)

sf = 2
png("demo-science-plots.png", width=1200*sf, height=600*sf, res=72*sf)
ggarrange( plotHypercube.bubbles(my.post.cgh, reorder=TRUE) + xlab("Order") + ylab("Aberration") + ggtitle("Ovarian cancer progression"),
           plotHypercube.bubbles(my.post.c4, reorder=TRUE) + xlab("Order") + ylab("C4 feature") + ggtitle("C4 photosynthesis"),
           plotHypercube.bubbles(my.post.malaria, reorder=TRUE, transpose=TRUE) + xlab("Symptom") + ylab("Order") + ggtitle("Severe malaria progression"),
           plotHypercube.bubbles(my.post.tools, reorder=TRUE, transpose=TRUE) + xlab("Tool use mode") + ylab("Order") + ggtitle("Tool use emergence"))
dev.off()
