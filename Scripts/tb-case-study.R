setwd("..")
source("hypertraps.R")
setwd("Scripts")

require(parallel)

# read and curate data on transitions for TB case study
tb.set = curate.tree("../Data/ng.2878-S2.txt", "../Data/tuberculosis-v5-header-19-29.csv")
tb.srcs = tb.set$srcs
tb.dests = tb.set$dests
tb.times = tb.set$times*1000

# for parallelisation -- with different control parameter "fork" this function will return different HyperTraPS experiments
parallel.fn = function(fork, srcs, dests, times) {
  if(fork == 1) { return(HyperTraPS(dests, initialstates = srcs, samplegap = 10, seed = 1, length = 4, kernel = 3))}
  if(fork == 2) { return(HyperTraPS(dests, initialstates = srcs, samplegap = 10, seed = 2, length = 4, kernel = 3))}
  if(fork == 3) { return(HyperTraPS(dests, initialstates = srcs, samplegap = 10, starttimes = times, endtimes = times, penalty = 1, seed = 1, length = 4, kernel = 3))}
  if(fork == 4) { return(HyperTraPS(dests, initialstates = srcs, samplegap = 10, starttimes = times, endtimes = times, penalty = 1, seed = 2, length = 4, kernel = 3))}
  if(fork == 5) { return(HyperTraPS(dests, initialstates = srcs, samplegap = 10, endtimes = times, seed = 1, length = 4, kernel = 3))}
  if(fork == 6) { return(HyperTraPS(dests, initialstates = srcs, samplegap = 10, endtimes = times, seed = 2, length = 4, kernel = 3))}
  if(fork == 7) { return(HyperTraPS(dests, initialstates = srcs, samplegap = 10, endtimes = times, penalty = 1, seed = 1, length = 4, kernel = 3))}
  if(fork == 8) { return(HyperTraPS(dests, initialstates = srcs, samplegap = 10, endtimes = times, penalty = 1, seed = 2, length = 4, kernel = 3))}
  if(fork == 9) { return(HyperTraPS(dests, initialstates = srcs, samplegap = 10, endtimes = times, penalty = 1, seed = 1, model = 3, length = 4, kernel = 3))}
  if(fork == 10) { return(HyperTraPS(dests, initialstates = srcs, samplegap = 10, endtimes = times, penalty = 1, seed = 2, model = 3, length = 4, kernel = 3))}
}

# run these experiments in parallel
n.fork = 10
expt.names = c("DT, seed 1", "DT, seed 2", "CT precise, penalty, seed 1", "CT precise, penalty, seed 2", 
               "CT wide, seed 1", "CT wide, seed 2", "CT wide, penalty, seed 1", "CT wide, penalty, seed 2",
               "CT wide, pen mod 3 seed 1", "CT wide, pen mod 3 seed 2")

parallelised.runs2 <- mcmapply(parallel.fn, fork=1:n.fork,
                               MoreArgs = list(src = tb.srcs,
                                               dests = tb.dests,
                                               times = tb.times),
                               SIMPLIFY = FALSE,
                               mc.cores = min(detectCores(), n.fork))

# labels for experiments and features

tb.names = readLines("../RawData/tb-labels.txt")

# pull the "bubble" plot data together for comparison
bdf = thdf = data.frame()
for(i in 1:n.fork) {
  this.post = parallelised.runs2[[i]]
  tmpdf = this.post$bubbles
  tmpdf$Expt=i
  bdf = rbind(bdf, tmpdf)
  tmpdf = this.post$timehists
  tmpdf$Expt=i
  thdf = rbind(thdf, tmpdf)
}
bdf$ExptLabel = expt.names[bdf$Expt]

# bubble plot comparing experiments
g.tb.bubbles = ggplot(bdf, aes(x=Time+Expt/10, y=OriginalIndex, size=Probability, color=ExptLabel)) +
  geom_point(alpha=0.5) + labs(x = "Ordering", y = "Feature", color = "Experiment") +
  theme_light()

sf = 2
png("plot-tb-bubbles-compare.png", width=600*sf, height=400*sf, res=72*sf)
print(g.tb.bubbles)
dev.off()

# look at likelihood traces
ggarrange(plotlist = lapply(parallelised.runs2, plotHypercube.lik.trace))

# influence graphs for all the L^2 cases
ggarrange(plotlist = lapply(parallelised.runs2[1:8], plotHypercube.influencegraph, thresh = 0, cv.thresh = 0.3))

# influence graphs for the L^3 cases
ggarrange(plotlist = lapply(parallelised.runs2[9:10], plotHypercube.influencegraph, thresh = 0, cv.thresh = 1))

# various plot collections
ggarrange(plotlist = lapply(parallelised.runs2, plotHypercube.influencegraph, thresh = 0, cv.thresh = 1.5))
#ggarrange(plotlist = lapply(parallelised.runs2, plotHypercube.influences, cv.thresh = 1.5, featurenames = tb.names))

# produce summary plots for each of the experiments
cv.ts = c(0.3, 0.3, 1.5, 1.5, 0.3, 0.3, 1.5, 1.5, 1.5, 1.5)
for(expt in 1:n.fork) {
  # focus on one penalised likelihood, continuous time case
  tb.post.1 = parallelised.runs2[[expt]]
  
  g.tb.graph = plotHypercube.sampledgraph2(tb.post.1, use.arc = FALSE, featurenames = tb.names, 
                                            edge.label.size=3, edge.label.angle = "none", node.labels=FALSE,
                                            no.times=TRUE, small.times=TRUE)
  g.tb.thist = plotHypercube.timehists(tb.post.1, featurenames = tb.names)
  g.tb.motifs = plotHypercube.motifs(tb.post.1, featurenames=tb.names)
  g.tb.influences = plotHypercube.influences(tb.post.1, featurenames=tb.names, reorder=TRUE,
                                             cv.thresh = cv.ts[expt])
  g.tb.tseries = plotHypercube.timeseries(tb.post.1, featurenames=tb.names)
  
  g.tb.influencegraph = plotHypercube.influencegraph(tb.post.1, featurenames=tb.names,
                                                     cv.thresh = cv.ts[expt])
  g.tb.motifseries = plotHypercube.motifseries(tb.post.1, t.set=c(0.001, 0.01, 0.1, 0.5, 1, 5, 10))
  g.blank = ggplot(data.frame()) + geom_blank()
  
  png(paste0("plot-tb-summary-", expt, ".png"), width=800*sf, height=800*sf, res=72*sf)
  print(ggarrange(ggarrange(g.blank, g.tb.graph + theme(legend.position="none"), widths=c(1,2), nrow=1), 
                  ggarrange(g.tb.influencegraph, g.tb.motifseries, widths=c(1,2), nrow=1), nrow=2))
  dev.off()
  
  png(paste0("plot-tb-summary-si-", expt, ".png"), width=1000*sf, height=400*sf, res=72*sf)
  print(ggarrange(g.tb.influences, g.tb.motifs, g.tb.tseries, nrow=1))
  dev.off()
}
