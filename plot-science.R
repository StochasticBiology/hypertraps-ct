library(ggplot2)
library(gridExtra)
source("plot-trans.R")

#### bubble and hypercube plots for TB experiments

fname = c("tb-dt-1", "tb-dt-2", "tb-ct-1", "tb-ct-2")
bdf = thdf = data.frame()
for(i in 1:length(fname)) {
  bubble.name = paste(c("Data/", fname[i], "-posterior.txt-bubbles.csv"), collapse="")
  tmpdf = read.csv(bubble.name)
  tmpdf$Expt=i
  bdf = rbind(bdf, tmpdf)
  thist.name = paste(c("Data/", fname[i], "-posterior.txt-timehists.csv"), collapse="")
  tmpdf = read.csv(thist.name)
  tmpdf$Expt=i
  thdf = rbind(thdf, tmpdf)
}
g.tb.bubbles = ggplot(bdf, aes(x=Time+Expt/10, y=OriginalIndex, size=Probability, color=factor(Expt))) +
  geom_point() +
  theme_light()

g.tb.thist = ggplot(thdf, aes(x=Time, y=Probability, fill=factor(OriginalIndex))) + 
  geom_col(position="dodge") + xlim(-0.1,1.05) + facet_wrap(~Expt, nrow=4) +
  theme_light() #+ scale_x_continuous(trans="log10")

trans.1 = read.csv("Data/tb-dt-1-trans.txt", sep=" ")
trans.s.1 = read.csv("Data/tb-dt-1-states.txt", sep=" ")
g.tb.cube = plot.hypercube3(trans.1, statesdf=trans.s.1, 
                node.labels = FALSE, seg.labels = TRUE, threshold=5e-2)

g.tb.summary = grid.arrange(g.tb.bubbles, g.tb.hist, nrow=1)
sf = 2
png("plot-science-tb.png", width=800*sf, height=600*sf, res=72*sf)
grid.arrange(g.tb.cube, g.tb.summary, nrow=2)
dev.off()

#### bubble and hypercube plots for MRO experiments

fname = c("mro-1", "mro-2", "mro-3", "mro-4")
bdf = thdf = data.frame()
for(i in 1:length(fname)) {
  bubble.name = paste(c("Data/", fname[i], "-posterior.txt-bubbles.csv"), collapse="")
  tmpdf = read.csv(bubble.name)
  tmpdf$Expt=i
  bdf = rbind(bdf, tmpdf)
  thist.name = paste(c("Data/", fname[i], "-posterior.txt-timehists.csv"), collapse="")
  tmpdf = read.csv(thist.name)
  tmpdf$Expt=i
  thdf = rbind(thdf, tmpdf)
}
g.mro.bubbles = ggplot(bdf, aes(x=Time+Expt/10, y=OriginalIndex, size=Probability, color=factor(Expt))) +
  geom_point() +
  theme_light()

g.mro.thist = ggplot(thdf, aes(x=Time, y=Probability, fill=factor(OriginalIndex))) + 
  geom_col(position="dodge") + xlim(-0.1,1.05) + facet_wrap(~Expt, nrow=4) +
  theme_light() #+ scale_x_continuous(trans="log10")

trans.1 = read.csv("Data/mro-1-trans.txt", sep=" ")
trans.s.1 = read.csv("Data/mro-1-states.txt", sep=" ")
g.mro.cube = plot.hypercube3(trans.1, statesdf=trans.s.1, 
                node.labels = FALSE, seg.labels = TRUE, threshold=1e-5)

g.mro.summary = grid.arrange(g.mro.bubbles, g.mro.hist, nrow=1)
sf = 2
png("plot-science-mro.png", width=800*sf, height=600*sf, res=72*sf)
grid.arrange(g.mro.cube, g.mro.summary, nrow=2)
dev.off()