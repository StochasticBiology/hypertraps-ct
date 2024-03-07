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
  if(fork == 1) { return(HyperTraPS(dests, initialstates = srcs, walkers = 20, seed = 1, samplegap = 10, length = 5, kernel = 3))}
  if(fork == 2) { return(HyperTraPS(dests, initialstates = srcs, walkers = 20, seed = 2, samplegap = 10, length = 5, kernel = 3))}
  if(fork == 3) { return(HyperTraPS(dests, initialstates = srcs, walkers = 20, penalty = 1, seed = 1, samplegap = 10, length = 5, kernel = 3))}
  if(fork == 4) { return(HyperTraPS(dests, initialstates = srcs, walkers = 20, penalty = 1, seed = 2, samplegap = 10, length = 5, kernel = 3))}
  if(fork == 5) { return(HyperTraPS(dests, initialstates = srcs, walkers = 20, penalty = 1, seed = 1, samplegap = 10, model = 3, length = 5, kernel = 3))}
  if(fork == 6) { return(HyperTraPS(dests, initialstates = srcs, walkers = 20, penalty = 1, seed = 2, samplegap = 10, model = 3, length = 5, kernel = 3))}
  if(fork == 7) { return(HyperTraPS(dests, initialstates = srcs, seed = 1, regularise = 1, samplegap = 10, length = 4, kernel = 3))}
  if(fork == 8) { return(HyperTraPS(dests, initialstates = srcs, seed = 2, regularise = 1, samplegap = 10, length = 4, kernel = 3))}
}


# run HyperTraPS and regularise via penalised likelihood

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

# various checks for consistency across random number seeds
ggarrange(plotHypercube.influences(parallelised.runs[[3]], cv.thresh = 2, featurenames = AML[[4]], 
                                   upper.right = TRUE, reorder = TRUE),
          plotHypercube.influences(parallelised.runs[[4]], cv.thresh = 2, featurenames = AML[[4]], 
                                   upper.right = TRUE, reorder = TRUE)
)

ggarrange(plotHypercube.influences(parallelised.runs[[1]], cv.thresh = 2, featurenames = AML[[4]], 
                                   upper.right = TRUE, reorder = TRUE),
          plotHypercube.influences(parallelised.runs[[2]], cv.thresh = 2, featurenames = AML[[4]], 
                                   upper.right = TRUE, reorder = TRUE)
)

##### plots for manuscript

# hypercube without timings
g.cancer.graph2 = plotHypercube.sampledgraph2(parallelised.runs[[4]], use.arc = FALSE, featurenames = AML[[4]], 
                                              edge.label.size=3, edge.label.angle = "along", node.labels=FALSE,
                                              no.times=TRUE, small.times=FALSE, thresh=0.008, truncate=5,
                                              use.timediffs = FALSE, edge.check.overlap = FALSE) +
  theme(legend.position="none") + coord_flip() + scale_y_reverse()

# hypercube with timings (messier)
g.cancer.graph2t = plotHypercube.sampledgraph2(parallelised.runs[[4]], use.arc = FALSE, featurenames = AML[[4]], 
                                               edge.label.size=3, edge.label.angle = "none", node.labels=FALSE,
                                               no.times=TRUE, small.times=TRUE,
                                               thresh=0.004, truncate=6,
                                               use.timediffs = FALSE) + 
  theme(legend.position="none") + coord_flip() + scale_y_reverse()

g.influences = plotHypercube.influences(parallelised.runs[[4]], cv.thresh = 2, featurenames = AML[[4]], 
                                        upper.right = TRUE, reorder = TRUE)

g.influence.graph = plotHypercube.influencegraph(parallelised.runs[[4]], 
                                                 cv.thresh = 1, thresh = 1, featurenames = AML[[4]])


g.motif = plotHypercube.motifs(parallelised.runs[[4]], 
                               label.size = 3,
                               featurenames = AML[[4]],
                               label.scheme = "sparse") + theme(legend.position="none")

sf = 3

png("cancer-post-si.png", width=400*sf, height=400*sf, res=72*sf)
print(g.influence.graph)
dev.off()


png("cancer-post-main-text.png", width=900*sf, height=700*sf, res=72*sf)
print(ggarrange(g.cancer.graph2,
                ggarrange(g.influences, g.motif, nrow=2, labels=c("B", "C")),
                labels=c("A", ""), nrow=1, widths=c(1,1.5)))
dev.off()

