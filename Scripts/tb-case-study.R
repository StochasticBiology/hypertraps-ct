setwd("..")
source("hypertraps.R")
setwd("Scripts")

require(parallel)

# read prepared data on transitions for TB case study
tb.tab = as.matrix(read.table("../Data/ng.2878-S2.txt-cooked.txt-data.txt"))
tb.srcs = tb.tab[seq(1, nrow(tb.tab), by = 2),]
tb.dests = tb.tab[seq(2, nrow(tb.tab), by = 2),]

tb.times = as.matrix(read.table("../Data/ng.2878-S2.txt-cooked.txt-datatime.txt"))*1000

# for parallelisation -- with different control parameter "fork" this function will return different HyperTraPS experiments
parallel.fn = function(fork, srcs, dests, times) {
  if(fork == 1) { return(HyperTraPS(dests, initialstates = srcs, seed = 1, length = 4, kernel = 4))}
  if(fork == 2) { return(HyperTraPS(dests, initialstates = srcs, seed = 2, length = 4, kernel = 4))}
  if(fork == 3) { return(HyperTraPS(dests, initialstates = srcs, penalty = 1, seed = 1, length = 4, kernel = 4))}
  if(fork == 4) { return(HyperTraPS(dests, initialstates = srcs, penalty = 1, seed = 2, length = 4, kernel = 4))}
  if(fork == 5) { return(HyperTraPS(dests, initialstates = srcs, endtimes = times, seed = 1, length = 4, kernel = 4))}
  if(fork == 6) { return(HyperTraPS(dests, initialstates = srcs, endtimes = times, seed = 2, length = 4, kernel = 4))}
  if(fork == 7) { return(HyperTraPS(dests, initialstates = srcs, endtimes = times, penalty = 1, seed = 1, length = 4, kernel = 4))}
  if(fork == 8) { return(HyperTraPS(dests, initialstates = srcs, endtimes = times, penalty = 1, seed = 2, length = 4, kernel = 4))}
}

# run these experiments in parallel
n.fork = 8
parallelised.runs <- mcmapply(parallel.fn, fork=1:n.fork,
                              MoreArgs = list(src = tb.srcs,
                                              dests = tb.dests,
                                              times = tb.times),
                              SIMPLIFY = FALSE,
                              mc.cores = min(detectCores(), n.fork))

# labels for experiments and features
expt.names = c("DT, seed 1", "DT, seed 2", "DT, penalty, seed 1", "DT, penalty, seed 2", "CT, seed 1", "CT, seed 2", "CT, penalty, seed 1", "CT, penalty, seed 2")
tb.names = readLines("../RawData/tb-labels.txt")

# pull the "bubble" plot data together for comparison
bdf = thdf = data.frame()
for(i in 1:n.fork) {
  this.post = parallelised.runs[[i]]
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

# focus on one penalised likelihood, continuous time case
tb.post.1 = parallelised.runs[[7]]

g.tb.graph = plotHypercube.sampledgraph2(tb.post.1, use.arc = FALSE, featurenames = tb.names, edge.label.size=3)
g.tb.graph2 = plotHypercube.sampledgraph2(tb.post.1, use.arc = FALSE, featurenames = tb.names, 
                                          edge.label.size=3, edge.label.angle = "none", node.labels=FALSE,
                                          no.times=TRUE, small.times=TRUE)
g.tb.thist = plotHypercube.timehists(tb.post.1, featurenames = tb.names)

png("plot-tb-summary.png", width=2000*sf, height=1000*sf, res=72*sf)
ggarrange(g.tb.graph2, g.tb.thist)
dev.off()

g.tb.motifs = plotHypercube.motifs(tb.post.1, featurenames=tb.names)
g.tb.influences = plotHypercube.influences(tb.post.1, featurenames=tb.names, reorder=TRUE)
g.tb.tseries = plotHypercube.timeseries(tb.post.1, featurenames=tb.names)

png("plot-tb-summary-si.png", width=800*sf, height=300*sf, res=72*sf)
ggarrange(g.tb.motifs, g.tb.influences)
dev.off()

png("plot-tb-alt.png", width=800*sf, height=1000*sf, res=72*sf)
ggarrange(g.tb.graph2 + theme(legend.position="none"), ggarrange(g.tb.motifs, g.tb.thist, nrow=1, labels=c("B", "C")), labels=c("A", ""), heights=c(2,1), nrow=2)
dev.off()

png("plot-tb-alt-2.png", width=800*sf, height=1000*sf, res=72*sf)
ggarrange(g.tb.graph2 + theme(legend.position="none"), ggarrange(g.tb.motifs, g.tb.influences, nrow=1, labels=c("B", "C")), labels=c("A", ""), 
          heights=c(1.5,1), nrow=2)
dev.off()

png("plot-tb-alt-si.png", width=800*sf, height=1000*sf, res=72*sf)
ggarrange(g.tb.bubbles, ggarrange(g.tb.influences, g.tb.tseries, nrow=1, labels=c("B", "C")), labels=c("A", ""), heights=c(2,1), nrow=2)
dev.off()

