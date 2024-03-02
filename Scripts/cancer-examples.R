# load the essentials
library(ggplot2)
library(ggpubr)
library(readxl)

setwd("..")
source("hypertraps.R")
setwd("Scripts/")

### Cancer case study 1

# pull the AML data object from TreeMHN
load("../Data/AML_tree_obj.RData")

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
srcs = matrix(unlist(lapply(strsplit(trans.df$before, ""), as.numeric)), ncol=AML[[1]], byrow=TRUE)
dests = matrix(unlist(lapply(strsplit(trans.df$after, ""), as.numeric)), ncol=AML[[1]], byrow=TRUE)

parallel.fn = function(fork, srcs, dests) {
  if(fork == 1) { return(HyperTraPS(dests, initialstates = srcs, seed = 1, samplegap = 10, length = 4, kernel = 3))}
  if(fork == 2) { return(HyperTraPS(dests, initialstates = srcs, seed = 2, samplegap = 10, length = 4, kernel = 3))}
  if(fork == 3) { return(HyperTraPS(dests, initialstates = srcs, penalty = 1, seed = 1, samplegap = 10, length = 4, kernel = 3))}
  if(fork == 4) { return(HyperTraPS(dests, initialstates = srcs, penalty = 1, seed = 2, samplegap = 10, length = 4, kernel = 3))}
  if(fork == 5) { return(HyperTraPS(dests, initialstates = srcs, penalty = 1, seed = 1, samplegap = 10, model = 3, length = 4, kernel = 3))}
  if(fork == 6) { return(HyperTraPS(dests, initialstates = srcs, penalty = 1, seed = 2, samplegap = 10, model = 3, length = 4, kernel = 3))}
  if(fork == 7) { return(HyperTraPS(dests, initialstates = srcs, seed = 1, regularise = 1, samplegap = 10, length = 4, kernel = 3))}
  if(fork == 8) { return(HyperTraPS(dests, initialstates = srcs, seed = 2, regularise = 1, samplegap = 10, length = 4, kernel = 3))}
}


# run HyperTraPS and regularise via penalised likelihood

# FIXME: a general concern about having, by default, seed = 1.
#  I think this is different from the way most code does this in R,
#  where the seed is not fixed, so different runs lead to possibly different
#  results, unless the user explicitly sets the seed.
#  It might also be worth thinking/illustrating how to do this if running in
#  parallel (getting an appropriate seed from R)

require(parallel)
# Create the data frame of options, print to check it is what we want
# and pass to mcmapply

nfork = 8
expt.names = c("seed 1", "seed 2", "pen seed 1", "pen seed 2", "pen mod 3 seed 1", "pen mod 3 seed 2", "reg seed 1", "reg seed 2")
parallelised.runs <- mcmapply(parallel.fn,
                              fork = 1:nfork,
                              MoreArgs = list(srcs=srcs,
                                              dests=dests),
                              SIMPLIFY = FALSE,
                              mc.cores = min(detectCores(), nfork))

cancer.post = parallelised.runs[[1]]
cancer.post.1 = parallelised.runs[[2]]
cancer.post.autoreg = parallelised.runs[[3]]
cancer.post.autoreg.1 = parallelised.runs[[4]]
cancer.post.3.autoreg = parallelised.runs[[5]]
cancer.post.3.autoreg.1 = parallelised.runs[[6]]

writeHyperinf(cancer.post, "cancer.post", postlabel="cancer.post", fulloutput = FALSE, regularised = FALSE)
writeHyperinf(cancer.post.autoreg, "cancer.post.autoreg", postlabel="cancer.post.autoreg", fulloutput = FALSE, regularised = FALSE)
writeHyperinf(cancer.post.3.autoreg, "cancer.post.3.autoreg", postlabel="cancer.post.3.autoreg", fulloutput = FALSE, regularised = FALSE)

# uncomment if we're reading previously computed experiments from file
#  cancer.post = readHyperinf("cancer.post", postlabel="cancer.post", fulloutput = FALSE, regularised = FALSE)
#  cancer.post.autoreg = readHyperinf("cancer.post.autoreg", postlabel="cancer.post.autoreg", fulloutput = FALSE, regularised = FALSE)
#  cancer.post.3.autoreg = readHyperinf("cancer.post.3.autoreg", postlabel="cancer.post.3.autoreg", fulloutput = FALSE, regularised = FALSE)

# hypercube without timings
g.cancer.graph2 = plotHypercube.sampledgraph2(cancer.post.autoreg, use.arc = FALSE, featurenames = AML[[4]], 
                            edge.label.size=3, edge.label.angle = "along", node.labels=FALSE,
                            no.times=TRUE, small.times=FALSE, thresh=0.006, truncate=6,
                            use.timediffs = FALSE, edge.check.overlap = FALSE) +
  theme(legend.position="none") + coord_flip() + scale_y_reverse()

# hypercube with timings (messier)
g.cancer.graph2t = plotHypercube.sampledgraph2(cancer.post.autoreg, use.arc = FALSE, featurenames = AML[[4]], 
                                               edge.label.size=3, edge.label.angle = "none", node.labels=FALSE,
                                               thresh=0.004, truncate=6,
                                               use.timediffs = FALSE) + 
  theme(legend.position="none") + coord_flip() + scale_y_reverse()

# create plots of influences under different regularisation protocols
plot.base = plotHypercube.influences(cancer.post, featurenames = AML[[4]], 
                                     use.final = TRUE, upper.right = TRUE, reorder = TRUE) +
  guides(alpha=FALSE)
plot.autoreg = plotHypercube.influences(cancer.post.autoreg, featurenames = AML[[4]], 
                                        upper.right = TRUE, reorder = TRUE) +
  guides(alpha = FALSE)

plot.base.pruned = plotHypercube.influences(cancer.post, featurenames = AML[[4]], 
                                            upper.right = TRUE, reorder = TRUE, cv.thresh = 0.5) +
  guides(alpha = FALSE)

plot.autoreg.pruned = plotHypercube.influences(cancer.post.autoreg, cv.thresh = 2, featurenames = AML[[4]], 
                                               upper.right = TRUE, reorder = TRUE) +
  guides(alpha = FALSE) 

# compare with and without penalised likelihood
ggarrange(plot.base.pruned, plot.autoreg.pruned)


sf = 3
png("cancer-post-autoreg.png", width=800*sf, height=800*sf, res=72*sf)
print(ggarrange(g.cancer.graph2 + theme(legend.position = "none"),
                plot.autoreg.pruned, nrow=2, labels=c("A", "B")))
dev.off()

png("cancer-post-si.png", width=400*sf, height=400*sf, res=72*sf)
print(plotHypercube.influencegraph(cancer.post.autoreg, 
                                   cv.thresh = 2, thresh = 1, featurenames = AML[[4]])
)
dev.off()

g.cancer.motif = plotHypercube.motifs(cancer.post.autoreg, 
                                      label.size = 3,
                                      featurenames = AML[[4]],
                                      label.scheme = "sparse") + theme(legend.position="none")

png("cancer-post-v3.png", width=800*sf, height=800*sf, res=72*sf)
ggarrange(g.cancer.graph2,
                     ggarrange(plot.autoreg.pruned, g.cancer.motif, nrow=2, labels=c("B", "C")),
                     labels=c("A", ""), nrow=1, widths=c(1,1.5))
dev.off()

png("cancer-post-v2.png", width=1200*sf, height=900*sf, res=72*sf)
print(ggarrange( ggarrange(g.cancer.graph2 + theme(legend.position = "none"),
                           plot.autoreg.pruned, nrow=1, labels=c("A", "B")),
                 plotHypercube.motifs(cancer.post.autoreg, featurenames = AML[[4]]) + theme(legend.position="none"),
                 labels=c("", "C"), nrow=2))
dev.off()


png("cancer-post-all.png", width=1200*sf, height=800*sf, res=72*sf)
print(ggarrange(plot.base, plot.autoreg, 
                plot.base.pruned, plot.autoreg.pruned, labels = c("A", "B", "C", "D")))
dev.off()


# put parameter loss details in regularisation output; add use.regularised to all plot functions
# consider clustering for influence maps

