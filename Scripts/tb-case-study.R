setwd("..")
source("hypertraps.R")
setwd("Scripts")

require(parallel)

# read and curate data on transitions for TB case study
tb.set = curate.tree("../Data/ng.2878-S2.txt", "../Data/tuberculosis-v5-header-19-29.csv")
tb.srcs = tb.set$srcs
tb.dests = tb.set$dests
tb.times = tb.set$times*1000
g.curated.tree = plotHypercube.curated.tree(tb.set)

# for parallelisation -- with different control parameter "fork" this function will return different HyperTraPS experiments
parallel.fn = function(fork, srcs, dests, times) {
  if(fork == 1) { return(HyperTraPS(dests, initialstates = srcs, seed = 1, length = 5, kernel = 3))}
  if(fork == 2) { return(HyperTraPS(dests, initialstates = srcs, seed = 2, length = 5, kernel = 3))}
  if(fork == 3) { return(HyperTraPS(dests, initialstates = srcs, starttimes = times, endtimes = times, penalty = 1, seed = 1, length = 5, kernel = 3))}
  if(fork == 4) { return(HyperTraPS(dests, initialstates = srcs, starttimes = times, endtimes = times, penalty = 1, seed = 2, length = 5, kernel = 3))}
  if(fork == 5) { return(HyperTraPS(dests, initialstates = srcs, starttimes=times*0.75, endtimes = times*1.25, penalty = 1, seed = 1, length = 5, kernel = 3))}
  if(fork == 6) { return(HyperTraPS(dests, initialstates = srcs, starttimes=times*0.75, endtimes = times*1.25, penalty = 1, seed = 2, length = 5, kernel = 3))}
  if(fork == 7) { return(HyperTraPS(dests, initialstates = srcs, endtimes = times, penalty = 1, seed = 1, length = 5, kernel = 3))}
  if(fork == 8) { return(HyperTraPS(dests, initialstates = srcs, endtimes = times, penalty = 1, seed = 2, length = 5, kernel = 3))}
  if(fork == 9) { return(HyperTraPS(dests, initialstates = srcs, starttimes=times, endtimes = times, penalty = 1, seed = 1, model = 3, length = 5, kernel = 3))}
  if(fork == 10) { return(HyperTraPS(dests, initialstates = srcs, starttimes=times, endtimes = times, penalty = 1, seed = 2, model = 3, length = 5, kernel = 3))}
  if(fork == 11) { return(HyperTraPS(dests, initialstates = srcs, starttimes=times, endtimes = times, penalty = 3, seed = 1, model = 3, length = 5, kernel = 3))}
  if(fork == 12) { return(HyperTraPS(dests, initialstates = srcs, starttimes=times, endtimes = times, penalty = 3, seed = 2, model = 3, length = 5, kernel = 3))}
}

# run these experiments in parallel
n.fork = 12

# labels for experiments and features
expt.names = c("DT 1", "DT 2", "CT 1", "CT 2", "CTu 1", "CTu 2", "CTw 1", "CTw 2", "CTm3 1", "CT m3 2", "CT m3p 1", "CT m3p 2")
tb.names = colnames(tb.set$data[2:ncol(tb.set$data)])
parallelised.runs2 <- mcmapply(parallel.fn, fork=1:n.fork,
                               MoreArgs = list(src = tb.srcs,
                                               dests = tb.dests,
                                               times = tb.times),
                               SIMPLIFY = FALSE,
                               mc.cores = min(detectCores(), n.fork))


####### some checks for consistency across different random seeds
ggarrange(plotHypercube.influencegraph(parallelised.runs2[[3]], featurenames=tb.names, cv.thresh = cv.thresh),
          plotHypercube.influencegraph(parallelised.runs2[[4]], featurenames=tb.names, cv.thresh = cv.thresh) )

ggarrange(plotHypercube.influencegraph(parallelised.runs2[[11]], featurenames=tb.names, cv.thresh = cv.thresh),
          plotHypercube.influencegraph(parallelised.runs2[[12]], featurenames=tb.names, cv.thresh = cv.thresh) )

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
g.tb.bubbles = ggplot(bdf, aes(x=Time+Expt/20, y=OriginalIndex, size=Probability, color=ExptLabel)) +
  geom_point(alpha=0.5) + labs(x = "Ordering", y = "Feature", color = "Experiment") +
  theme_light()

sf = 2
png("plot-tb-bubbles-compare.png", width=600*sf, height=400*sf, res=72*sf)
print(g.tb.bubbles)
dev.off()

# look at likelihood traces
ggarrange(plotlist = lapply(parallelised.runs2, plotHypercube.lik.trace))

##########

# main text figure: penalised likelihood, precisely specified timings
png("plot-tb-main-text.png", width=800*sf, height=800*sf, res=72*sf)
ggarrange(
  g.curated.tree,
  plotHypercube.sampledgraph2(parallelised.runs2[[3]], use.arc = FALSE, featurenames = tb.names, 
                              edge.label.size=3, edge.label.angle = "none", node.labels=FALSE,
                              no.times=TRUE, small.times=TRUE) + theme(legend.position="none"),
  plotHypercube.influencegraph(parallelised.runs2[[3]], cv.thresh=0.3, featurenames=tb.names),
  plotHypercube.motifseries(parallelised.runs2[[4]], t.set=c(1, 5, 10, 15, 20, 50), label.size=2.5, thresh = 0.02),
  widths=c(1,2,1,2),
  labels=c("A", "B", "C", "D")
  
)
dev.off()

# SI figure: alternative plots, discrete time, uncertain time specification, tripletwise model

png("plot-tb-si.png", width=800*sf, height=1000*sf, res=72*sf)
ggarrange( 
  ggarrange(plotHypercube.influences(parallelised.runs2[[3]], featurenames=tb.names, reorder=TRUE,
                                     cv.thresh = 0.3),
            plotHypercube.motifs(parallelised.runs2[[3]], featurenames=tb.names) + theme(legend.position="none"),
            plotHypercube.timeseries(parallelised.runs2[[3]], featurenames=tb.names), 
            nrow=1,
            labels=c("A", "B", "C")),
  
  ggarrange(
    plotHypercube.sampledgraph2(parallelised.runs2[[1]], use.arc = FALSE, featurenames = tb.names, 
                                edge.label.size=3, edge.label.angle = "none", node.labels=FALSE,
                                no.times=TRUE, small.times=FALSE) + theme(legend.position="none"),
    plotHypercube.influencegraph(parallelised.runs2[[7]], cv.thresh=0.3, featurenames=tb.names), 
    widths=c(2,1),
    nrow=1,
    labels=c("D", "E")),
  
  ggarrange(
    plotHypercube.sampledgraph2(parallelised.runs2[[7]], use.arc = FALSE, featurenames = tb.names, 
                                edge.label.size=3, edge.label.angle = "none", node.labels=FALSE,
                                no.times=TRUE, small.times=TRUE) + theme(legend.position="none"),
    plotHypercube.influencegraph(parallelised.runs2[[11]], cv.thresh=0.3, featurenames=tb.names),
    widths=c(2,1),
    nrow=1,
    labels=c("F", "G")),
  nrow=3
)
dev.off()


# other plots for curiosity

# influence graphs for all the L^2 cases
ggarrange(plotlist = lapply(parallelised.runs2[1:8], plotHypercube.influencegraph, thresh = 0, cv.thresh = 0.3))

# influence graphs for the L^3 cases
ggarrange(plotlist = lapply(parallelised.runs2[9:10], plotHypercube.influencegraph, thresh = 0, cv.thresh = 1))

# various plot collections
ggarrange(plotlist = lapply(parallelised.runs2, plotHypercube.influencegraph, thresh = 0, cv.thresh = 1.5))
#ggarrange(plotlist = lapply(parallelised.runs2, plotHypercube.influences, cv.thresh = 1.5, featurenames = tb.names))

# produce summary plots for each of the experiments
cv.ts = rep(0.4, n.fork)
for(expt in 1:n.fork) {
  # focus on one penalised likelihood, continuous time case
  tb.post.1 = parallelised.runs2[[expt]]
  
  g.tb.graph = plotHypercube.sampledgraph2(tb.post.1, use.arc = FALSE, featurenames = tb.names, 
                                           edge.label.size=3, edge.label.angle = "none", node.labels=FALSE,
                                           no.times=TRUE, small.times=TRUE)
  g.tb.thist = plotHypercube.timehists(tb.post.1, featurenames = tb.names)
  g.tb.motifs = plotHypercube.motifs(tb.post.1, featurenames=tb.names)
  
  g.tb.tseries = plotHypercube.timeseries(tb.post.1, featurenames=tb.names)
  
  g.tb.influencegraph = plotHypercube.influencegraph(tb.post.1, featurenames=tb.names,
                                                     cv.thresh = cv.ts[expt])
  if(expt <= 4 | expt >= 11) {
    this.t.set=c(0.001, 0.01, 0.1, 0.5, 1, 5, 10)
  } else {
    this.t.set=c(1e-5,1e-4,2e-4,5e-4,1e-3,1e-2)
  }
  g.tb.motifseries = plotHypercube.motifseries(tb.post.1, t.set=this.t.set)
  
  png(paste0("plot-tb-summary-", expt, ".png"), width=800*sf, height=800*sf, res=72*sf)
  print(ggarrange(ggarrange(g.curated.tree, g.tb.graph + theme(legend.position="none"), widths=c(1,2), nrow=1, labels=c("A", "B")), 
                  ggarrange(g.tb.influencegraph, g.tb.motifseries, widths=c(1,2), nrow=1, labels=c("C", "D")), nrow=2))
  dev.off()
  
  if(tb.post.1$model == 2) {
    g.tb.influences = plotHypercube.influences(tb.post.1, featurenames=tb.names, reorder=TRUE,
                                               cv.thresh = cv.ts[expt])
    png(paste0("plot-tb-summary-si-", expt, ".png"), width=1000*sf, height=400*sf, res=72*sf)
    print(ggarrange(g.tb.influences, g.tb.motifs, g.tb.tseries, nrow=1, labels=c("A", "B", "C")))
    dev.off()
  }
}

summary.list = lapply(parallelised.runs2, plotHypercube.lik.trace)
ggarrange(plotlist=summary.list)




