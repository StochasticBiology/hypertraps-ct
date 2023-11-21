setwd("..")
source("hypertraps.R")
setwd("Scripts")

#### bubble and hypercube plots for TB experiments

sims.completed = FALSE

tb.names = readLines("../RawData/tb-labels.txt")
fname = c("tb-dt-1", "tb-dt-2", "tb-ct-1", "tb-ct-2")
bdf = thdf = data.frame()
for(i in 1:length(fname)) {
  this.name = paste(c("../Data/", fname[i]), collapse="")
  if(sims.completed == FALSE) {
    this.post = readHyperinf(this.name)
    this.post = PosteriorAnalysis(this.post,featurenames = tb.names)
  } else {
    this.post = readHyperinf(this.name, postlabel = this.name)
  }
  tmpdf = this.post$bubbles
  tmpdf$Expt=i
  bdf = rbind(bdf, tmpdf)
  tmpdf = this.post$timehists
  tmpdf$Expt=i
  thdf = rbind(thdf, tmpdf)
}
g.tb.bubbles = ggplot(bdf, aes(x=Time+Expt/10, y=OriginalIndex, size=Probability, color=factor(Expt))) +
  geom_point() +
  theme_light()

tb.post.1 = readHyperinf("../Data/tb-ct-1")
tb.post.1 = PosteriorAnalysis(tb.post.1, featurenames=tb.names)

plotHypercube.sampledgraph2(tb.post.1, use.arc = FALSE)
g.tb.graph = plotHypercube.sampledgraph2(tb.post.1, use.arc = FALSE, feature.names = tb.names, edge.label.size=3)
g.tb.thist = plotHypercube.timehists(tb.post.1, feature.names = tb.names)

png("plot-tb-summary.png", width=2000*sf, height=1000*sf, res=72*sf)
ggarrange(g.tb.graph, g.tb.thist)
dev.off()

g.tb.motifs = plotHypercube.motifs(tb.post.1, feature.names=tb.names)
g.tb.influences = plotHypercube.influences(tb.post.1, feature.names=tb.names)
g.tb.tseries = plotHypercube.timeseries(tb.post.1, feature.names=tb.names)

png("plot-tb-summary-si.png", width=800*sf, height=300*sf, res=72*sf)
ggarrange(g.tb.motifs, g.tb.influences)
dev.off()

#### bubble and hypercube plots for MRO experiments

sims.completed = FALSE

mro.names = readLines("../RawData/mro-labels.txt")
fname = c("mro-3", "mro-4") 
bdf = thdf = data.frame()
for(i in 1:length(fname)) {
  this.name = paste(c("../Data/", fname[i]), collapse="")
  if(sims.completed == FALSE) {
    this.post = readHyperinf(this.name)
    this.post = PosteriorAnalysis(this.post,featurenames = tb.names)
  } else {
    this.post = readHyperinf(this.name, postlabel = this.name)
  }
  tmpdf = this.post$bubbles
  tmpdf$Expt=i
  bdf = rbind(bdf, tmpdf)
  tmpdf = this.post$timehists
  tmpdf$Expt=i
  thdf = rbind(thdf, tmpdf)
}
g.mro.bubbles = ggplot(bdf, aes(x=Time+Expt/10, y=OriginalIndex, size=Probability, color=factor(Expt))) +
  geom_point() +
  theme_light()

mro.post.1 = readHyperinf("../Data/mro-3")
mro.post.1 = PosteriorAnalysis(mro.post.1, featurenames=mro.names)

plotHypercube.sampledgraph2(mro.post.1, use.arc = FALSE)
g.mro.graph = plotHypercube.sampledgraph2(mro.post.1, use.arc = FALSE, feature.names = mro.names, edge.label.size = 3)
g.mro.thist = plotHypercube.timehists(mro.post.1, feature.names = mro.names, t.thresh = 2)

png("plot-mro-summary.png", width=2000*sf, height=1000*sf, res=72*sf)
ggarrange(g.mro.graph, g.mro.thist)
dev.off()

g.mro.motifs = plotHypercube.motifs(mro.post.1, feature.names=mro.names)
g.mro.influences = plotHypercube.influences(mro.post.1, feature.names=mro.names)
g.mro.tseries = plotHypercube.timeseries(mro.post.1, feature.names=mro.names)

png("plot-mro-summary-si.png", width=800*sf, height=300*sf, res=72*sf)
ggarrange(g.mro.motifs, g.mro.influences)
dev.off()

png("plot-tb-alt.png", width=800*sf, height=1000*sf, res=72*sf)
ggarrange(g.tb.graph + theme(legend.position="none"), ggarrange(g.tb.motifs, g.tb.thist, nrow=1, labels=c("B", "C")), labels=c("A", ""), heights=c(2,1), nrow=2)
dev.off()

png("plot-tb-alt-si.png", width=800*sf, height=1000*sf, res=72*sf)
ggarrange(g.tb.bubbles, ggarrange(g.tb.influences, g.tb.tseries, nrow=1, labels=c("B", "C")), labels=c("A", ""), heights=c(2,1), nrow=2)
dev.off()

png("plot-mro-alt.png", width=800*sf, height=1000*sf, res=72*sf)
ggarrange(g.mro.graph + theme(legend.position="none"), ggarrange(g.mro.motifs, g.mro.thist, nrow=1, labels=c("B", "C")), labels=c("A", ""), heights=c(2,1), nrow=2)
dev.off()

png("plot-mro-alt-si.png", width=800*sf, height=1000*sf, res=72*sf)
ggarrange(g.mro.bubbles, ggarrange(g.mro.influences, g.mro.tseries, nrow=1, labels=c("B", "C")), labels=c("A", ""), heights=c(2,1), nrow=2)
dev.off()

