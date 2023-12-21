# load the essentials
library(ggplot2)
library(ggpubr)
library(readxl)

setwd("..")
source("hypertraps.R")

### Cancer case study 1

# pull the AML data object from TreeMHN
load("RawData/AML_tree_obj.RData")

# helper function for shifting bits in a character string
new.state = function(old.state, locus) {
  tmp = old.state
  substr(tmp, locus, locus) = "1"
  return(tmp)
}

null.str = paste(rep("0", AML[[1]]), collapse="")

# read the data from the AML structure
trans.df = data.frame()
for(aml.list in 1:length(AML[[5]])) {
  # pull a specific cancer phylogeny from the dataset
  sub.ex =  AML$tree_df[AML$tree_df$Tree_ID == match(AML[[5]][aml.list], AML$tree_labels),]
  print(c(aml.list, nrow(sub.ex) ))
  # initialise data structure for this phylogeny
  sub.ex$state = NA
  sub.ex$state[sub.ex$Node_ID==1] = null.str
  # initialise to-do list for tree traversal
  to.do = sub.ex$Node_ID[sub.ex$Parent_ID==1 & sub.ex$Node_ID != 1]
  while(length(to.do) > 0) { 
    # go through tree, following descendant paths
    new.to.do = c()
    for(i in 1:length(to.do)) {
      # get before and after state for this edge
      b.state = sub.ex$state[sub.ex$Node_ID==sub.ex$Parent_ID[to.do[i]]]
      a.state = new.state(b.state, sub.ex$Mutation_ID[to.do[i]])  
      sub.ex$state[to.do[i]] = a.state
      # add to growing data structure and update to-do list to this node's daughters
      trans.df = rbind(trans.df, data.frame(before = b.state, after = a.state))
      new.to.do = c(new.to.do, sub.ex$Node_ID[sub.ex$Parent_ID==to.do[i] & sub.ex$Node_ID != 1])
    }
    to.do = unique(new.to.do)
  }
}

# put data into a form we can work with
befores = matrix(unlist(lapply(strsplit(trans.df$before, ""), as.numeric)), ncol=AML[[1]], byrow=TRUE)
afters = matrix(unlist(lapply(strsplit(trans.df$after, ""), as.numeric)), ncol=AML[[1]], byrow=TRUE)

if(run.simulations == TRUE) {
  # run HyperTraPS and regularise
  cancer.post = HyperTraPS(afters, initialstates = befores, length = 4, kernel = 3)
  cancer.post.autoreg = HyperTraPS(afters, initialstates = befores, length = 4, kernel = 3, autoregularise = 1)
  cancer.post.reg = HyperTraPS(afters, initialstates = befores, length = 4, kernel = 3, walkers=20, regularise = 1)
  cancer.post.sa = HyperTraPS(afters, initialstates = befores, sa = 1, length = 3, kernel = 3)
  cancer.post.sa.autoreg = HyperTraPS(afters, initialstates = befores, length = 3, kernel = 3, autoregularise = 1)

  # save output
  writeHyperinf(cancer.post, "cancer.post", postlabel="cancer.post", fulloutput = FALSE, regularised = FALSE)
  writeHyperinf(cancer.post.reg, "cancer.post.reg", postlabel="cancer.post.reg", fulloutput = FALSE, regularised = TRUE)
  writeHyperinf(cancer.post.autoreg, "cancer.post.autoreg", postlabel="cancer.post.autoreg", fulloutput = FALSE, regularised = FALSE)
  writeHyperinf(cancer.post.sa, "cancer.post.sa", postlabel="cancer.post.sa", fulloutput = FALSE, regularised = FALSE)
  writeHyperinf(cancer.post.sa.autoreg, "cancer.post.sa.autoreg", postlabel="cancer.post.sa.autoreg", fulloutput = FALSE, regularised = FALSE)
}

# example subsequent read
cancer.post = readHyperinf("cancer.post", postlabel="cancer.post", fulloutput = FALSE, regularised = FALSE)
cancer.post.reg = readHyperinf("cancer.post.reg", postlabel="cancer.post.reg", fulloutput = FALSE, regularised = TRUE)
cancer.post.autoreg = readHyperinf("cancer.post.autoreg", postlabel="cancer.post.autoreg", fulloutput = FALSE, regularised = FALSE)
cancer.post.sa = readHyperinf("cancer.post.sa", postlabel="cancer.post.sa", fulloutput = FALSE, regularised = FALSE)
cancer.post.sa.autoreg = readHyperinf("cancer.post.sa.autoreg", postlabel="cancer.post.sa.autoreg", fulloutput = FALSE, regularised = FALSE)

plotHypercube.influences(cancer.post, feature.names = AML[[4]], 
                         upper.right = TRUE, reorder = TRUE)
plotHypercube.influences(cancer.post.reg, feature.names = AML[[4]], 
                         use.regularised = TRUE, upper.right = TRUE, reorder = TRUE)
plotHypercube.influences(cancer.post.autoreg, feature.names = AML[[4]], 
                         upper.right = TRUE, reorder = TRUE)
plotHypercube.influences(cancer.post.sa, feature.names = AML[[4]], 
                         upper.right = TRUE, reorder = TRUE)
plotHypercube.influences(cancer.post.sa.autoreg, feature.names = AML[[4]], 
                         upper.right = TRUE, reorder = TRUE)

ggarrange(plotHypercube.influences(cancer.post.reg, feature.names = AML[[4]], 
                                   use.regularised = TRUE, upper.right = TRUE, reorder = TRUE),
          plotHypercube.influences(cancer.post.sa.autoreg, feature.names = AML[[4]], 
                                   upper.right = TRUE, reorder = TRUE),
          nrow=1, ncol=2)
### do this just with best parameterisation?

# plot influences between genes
plotHypercube.influences(cancer.post.reg, feature.names = AML[[4]], 
                         use.regularised = TRUE, upper.right = TRUE, reorder = TRUE)

# more samples from the regularised parameterisation
cancer.more.samples = PosteriorAnalysis(cancer.post.reg, use_regularised = TRUE, samples_per_row = 1000)

# various plots
plotHypercube.sampledgraph2(cancer.post.reg, use.arc = FALSE, node.labels = FALSE, no.times = TRUE)
plotHypercube.sampledgraph3(cancer.post.reg, use.arc = FALSE, node.labels = FALSE, 
                            no.times = TRUE, truncate = 6, edge.label.size=4, thresh = 0.01,
                            feature.names = labels)
plotHypercube.sampledgraph3(cancer.more.samples, use.arc = FALSE, node.labels = FALSE, 
                            no.times = TRUE, truncate = 6, edge.label.size=2, thresh = 0.02,
                            feature.names = labels)


# put parameter loss details in regularisation output; add use.regularised to all plot functions
# consider clustering for influence maps

### Cancer case study 2

# read a big list of mutations across samples
big.c.df = as.data.frame(read_excel("RawData/Supplementary Table 14.Driver.Events.By.Mutation.Type.01052015.v2.xlsx", sheet=4))
genes = unique(big.c.df$Gene)
samples = unique(big.c.df$Sample)

# initialise a dataframe to store binary states
data.df <- data.frame(matrix(0, ncol = length(genes), nrow = length(samples)))
colnames(data.df) <- genes
data.df$Sample = samples
# populate this dataframe when we find a given mutation for a given sample
for(i in 1:nrow(big.c.df)) {
  sref = which(samples == big.c.df$Sample[i])
  gref = which(genes == big.c.df$Gene[i])
  data.df[sref, gref] = 1
}

# put into matrix form and run HyperTraPS
big.c.m = as.matrix(data.df[1:(ncol(data.df)-1)])
big.c.post = HyperTraPS(big.c.m, kernel = 3)
writeHyperinf(big.c.post, "big-c-post")

# plot influences between genes
plotHypercube.influences(big.c.post, feature.names = genes, 
                         upper.right = FALSE, reorder = TRUE)

# more sampling from posterior
big.c.more = PosteriorAnalysis(big.c.post, samples_per_row = 100)

# followup plots
plotHypercube.sampledgraph3(big.c.post, use.arc = FALSE, node.labels = FALSE, feature.names = genes, 
                            no.times = TRUE, thresh=0.02, truncate = 3)

plotHypercube.sampledgraph3(big.c.more, use.arc = FALSE, node.labels = FALSE, feature.names = genes, 
                            no.times = TRUE, thresh=0.005, truncate = 3)

