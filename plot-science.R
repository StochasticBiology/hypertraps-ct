#### x range, scale and log type to be edited

library(ggplot2)
library(gridExtra)
library(igraph)
library(ggraph)
library(stringr)
#source("plot-trans.R")

DecToBin <- function(x, len) {
  s = c()
  for(j in (len-1):0)
  {
    if(x >= 2**j) { s=c(s,1); x = x-2**j } else { s=c(s,0)}
  }
  return(paste(s, collapse=""))
}

#### bubble and hypercube plots for TB experiments

fname = c("tb-dt-1", "tb-dt-2", "tb-ct-1", "tb-ct-2")
bdf = thdf = data.frame()
for(i in 1:length(fname)) {
  bubble.name = paste(c("../Data/", fname[i], "-bubbles.csv"), collapse="")
  tmpdf = read.csv(bubble.name)
  tmpdf$Expt=i
  bdf = rbind(bdf, tmpdf)
  thist.name = paste(c("../Data/", fname[i], "-timehists.csv"), collapse="")
  tmpdf = read.csv(thist.name)
  tmpdf$Expt=i
  thdf = rbind(thdf, tmpdf)
}
g.tb.bubbles = ggplot(bdf, aes(x=Time+Expt/10, y=OriginalIndex, size=Probability, color=factor(Expt))) +
  geom_point() +
  theme_light()

thdfp = data.frame()
thresh = 20
for(i in c(3,4)) {
  for(j in unique(thdf$OriginalIndex)) {
    sub = thdf[thdf$Expt==i & thdf$OriginalIndex == j & thdf$Time < thresh,]
    sub1 = thdf[thdf$Expt==i & thdf$OriginalIndex == j & thdf$Time >= thresh,]
    thdfp = rbind(thdfp, sub)
    thdfp = rbind(thdfp, data.frame(Expt=i, OriginalIndex=j, Time=thresh, Probability=sum(sub1$Probability)))
  }
}

g.tb.thist = ggplot(thdfp[thdfp$Time < thresh,], aes(x=log(Time+1), y=Probability, color=factor(Expt))) + 
  #geom_col(position="dodge") + xlim(-0.1,thresh+0.5) + facet_wrap(~OriginalIndex, ncol=2, scales="free") +
  geom_line() + xlim(-0.1,log(thresh+1)) + facet_wrap(~OriginalIndex, nrow=2) +
  theme_light() #+ scale_x_continuous(trans="log10")

g.tb.thist2 = ggplot(thdfp[thdfp$Time == thresh,], aes(x=OriginalIndex, y=Probability, fill=factor(Expt))) + 
  #geom_col(position="dodge") + xlim(-0.1,thresh+0.5) + facet_wrap(~OriginalIndex, ncol=2, scales="free") +
  geom_col(position="dodge") + 
  theme_light() #+ scale_x_continuous(trans="log10")


#g.tb.thist = ggplot(thdfp, aes(x=Time, y=Probability, color=factor(Expt), fill=factor(Expt))) + 
  #geom_col(position="dodge") +
#  geom_line() + 
#  xlim(-0.1,10.05) + facet_wrap(~OriginalIndex, ncol=2) +
#  theme_light() #+ scale_x_continuous(trans="log10")

ggplot(thdfp[thdfp$Expt>2,], aes(x=Time, y=Probability, color=factor(Expt))) + 
  geom_line() + xlim(-0.1,11) +
  geom_col(data=thdfp[thdfp$Expt>2 && thdfp$Time==10,]) + facet_wrap(~OriginalIndex, ncol=2) +
  theme_light() #+ scale_x_continuous(trans="log10")

# read transition and state probabilities
trans.1 = read.csv("../Data/tb-dt-1-trans.txt", sep=" ")
trans.1 = trans.1[!is.nan(trans.1$Probability),]

trans.s.1 = read.csv("../Data/tb-dt-1-states.txt", sep=" ")
trans.s.1 = trans.s.1[!is.nan(trans.s.1$Probability),]

# set up metadata for ggraph plot
trans.1$Flux = trans.1$Probability*trans.s.1$Probability[trans.1$From+1]

bigL = 10
trans.p = trans.1[trans.1$Flux > 0.05,]
trans.g = graph_from_data_frame(trans.p)
bs = unlist(lapply(as.numeric(V(trans.g)$name), DecToBin, len=bigL))
V(trans.g)$binname = bs
layers = str_count(bs, "1")
g.tb.cube = ggraph(trans.g, layout="sugiyama", layers=layers) + geom_edge_link(aes(edge_width=Flux, edge_alpha=Flux)) + 
  geom_node_point() + geom_node_label(aes(label=binname),size=2) +
  scale_edge_width(limits=c(0,NA)) + scale_edge_alpha(limits=c(0,NA)) +
  theme_graph() #aes(label=bs)) + theme_graph() 


#g.tb.cube = plot.hypercube3(trans.1, statesdf=trans.s.1, 
#                node.labels = FALSE, seg.labels = TRUE, threshold=5e-2)

g.tb.summary = grid.arrange(g.tb.bubbles, g.tb.thist, g.tb.thist2, nrow=3)

######## routes analysis TB
routes = read.table("../Data/tb-ct-1-routes.txt")
routetimes = read.table("../Data/tb-ct-1-times.txt")

# motif plot
rdf = data.frame()
for(j in 1:ncol(routes)) {
  startprob = 0
  for(i in 0:max(routes)) {
    thisprob = length(which(routes[,j]==i))/nrow(routes)
    rdf = rbind(rdf, data.frame(Index=i, Time=j, Start=startprob, End=startprob+thisprob, Probability=thisprob))
    startprob = startprob+thisprob
  }
}
g.tb.motifs = ggplot(rdf) + geom_rect(aes(xmin=Time-0.5,xmax=Time+0.5,ymin=Start,ymax=End,fill=factor(Index))) +
  geom_text(aes(x=Time,y=(Start+End)/2,label=Index), color="#FFFFFF") + ylab("Probability") + theme_light()

# time series illustration
rtdf = data.frame()
for(i in 1:1000) { #nrow(routes)) {
  prevtime = 0
  for(j in 1:ncol(routes)) {
    rtdf = rbind(rtdf, data.frame(Run=i, Step=j, Index=routes[i,j], PrevTime=prevtime, Time=routetimes[i,j]))
    prevtime = routetimes[i,j]
  }
}
g.tb.ts = ggplot(rtdf) + geom_segment(aes(x=PrevTime,xend=Time,y=Step-1,yend=Step,color=factor(Index)), alpha=0.5) +
  scale_x_continuous(trans="log") + theme_light()

sf = 2
png("plot-science-tb.png", width=1200*sf, height=1200*sf, res=72*sf)
grid.arrange(g.tb.cube, g.tb.summary, g.tb.motifs, g.tb.ts, nrow=2)
dev.off()

#### bubble and hypercube plots for MRO experiments

fname = c("mro-3", "mro-4") #, "mro-5", "mro-6")
bdf = thdf = data.frame()
for(i in 1:length(fname)) {
  bubble.name = paste(c("../Data/", fname[i], "-bubbles.csv"), collapse="")
  tmpdf = read.csv(bubble.name)
  tmpdf$Expt=i
  bdf = rbind(bdf, tmpdf)
  thist.name = paste(c("../Data/", fname[i], "-timehists.csv"), collapse="")
  tmpdf = read.csv(thist.name)
  tmpdf$Expt=i
  thdf = rbind(thdf, tmpdf)
}
g.mro.bubbles = ggplot(bdf, aes(x=Time+Expt/10, y=OriginalIndex, size=Probability, color=factor(Expt))) +
  geom_point() +
  theme_light()

thdfp = data.frame()
thresh = 10
for(i in c(3,4)) {
  for(j in unique(thdf$OriginalIndex)) {
    sub = thdf[thdf$Expt==i & thdf$OriginalIndex == j & thdf$Time < thresh,]
    sub1 = thdf[thdf$Expt==i & thdf$OriginalIndex == j & thdf$Time >= thresh,]
    thdfp = rbind(thdfp, sub)
    thdfp = rbind(thdfp, data.frame(Expt=i, OriginalIndex=j, Time=thresh, Probability=sum(sub1$Probability)))
  }
}

g.mro.thist = ggplot(thdfp[thdfp$Time < thresh,], aes(x=log(Time+1), y=Probability, fill=factor(Expt))) + 
  #geom_col(position="dodge") + xlim(-0.1,thresh+0.5) + facet_wrap(~OriginalIndex, ncol=2, scales="free") +
    geom_col(position="dodge") + xlim(-0.1,log(thresh+1)) + facet_wrap(~OriginalIndex, ncol=2) +
  theme_light() #+ scale_x_continuous(trans="log10")

g.mro.thist = ggplot(thdfp[thdfp$Time < thresh,], aes(x=log(Time+1), y=Probability, color=factor(Expt))) + 
  #geom_col(position="dodge") + xlim(-0.1,thresh+0.5) + facet_wrap(~OriginalIndex, ncol=2, scales="free") +
  geom_line() + xlim(-0.1,log(thresh+1)) + facet_wrap(~OriginalIndex, nrow=2) +
  theme_light() #+ scale_x_continuous(trans="log10")


g.mro.thist2 = ggplot(thdfp[thdfp$Time == thresh,], aes(x=OriginalIndex, y=Probability, fill=factor(Expt))) + 
  #geom_col(position="dodge") + xlim(-0.1,thresh+0.5) + facet_wrap(~OriginalIndex, ncol=2, scales="free") +
  geom_col(position="dodge") + 
  theme_light() #+ scale_x_continuous(trans="log10")


grid.arrange(g.mro.thist, g.mro.thist2, nrow=2)

trans.1 = read.csv("../Data/mro-1-trans.txt", sep=" ")
trans.s.1 = read.csv("../Data/mro-1-states.txt", sep=" ")
# set up metadata for ggraph plot
trans.1$Flux = trans.1$Probability*trans.s.1$Probability[trans.1$From+1]

bigL = 10
trans.p = trans.1[trans.1$Flux > 0.05,]
trans.g = graph_from_data_frame(trans.p)
bs = unlist(lapply(as.numeric(V(trans.g)$name), DecToBin, len=bigL))
V(trans.g)$binname = bs
layers = str_count(bs, "1")
g.mro.cube = ggraph(trans.g, layout="sugiyama", layers=layers) + geom_edge_link(aes(edge_width=Flux, edge_alpha=Flux)) + 
  geom_node_point() + geom_node_label(aes(label=binname),size=2) +
  scale_edge_width(limits=c(0,NA)) + scale_edge_alpha(limits=c(0,NA)) +
  theme_graph() #aes(label=bs)) + theme_graph() 

#g.mro.cube = plot.hypercube3(trans.1, statesdf=trans.s.1, 
#                node.labels = FALSE, seg.labels = TRUE, threshold=1e-5)

g.mro.summary = grid.arrange(g.mro.bubbles, g.mro.thist, g.mro.thist2, ncol=1)

#### routes analysis MRO
routes = read.table("../Data/mro-3-routes.txt")
routetimes = read.table("../Data/mro-3-times.txt")

# motif plot
rdf = data.frame()
for(j in 1:ncol(routes)) {
  startprob = 0
  for(i in 0:max(routes)) {
    thisprob = length(which(routes[,j]==i))/nrow(routes)
    rdf = rbind(rdf, data.frame(Index=i, Time=j, Start=startprob, End=startprob+thisprob, Probability=thisprob))
    startprob = startprob+thisprob
  }
}
g.mro.motifs = ggplot(rdf) + geom_rect(aes(xmin=Time-0.5,xmax=Time+0.5,ymin=Start,ymax=End,fill=factor(Index))) +
  geom_text(aes(x=Time,y=(Start+End)/2,label=Index), color="#FFFFFF") + ylab("Probability") + theme_light()

# time series illustration
rtdf = data.frame()
for(i in 1:1000) { #nrow(routes)) {
  prevtime = 0
  for(j in 1:ncol(routes)) {
    rtdf = rbind(rtdf, data.frame(Run=i, Step=j, Index=routes[i,j], PrevTime=prevtime, Time=routetimes[i,j]))
    prevtime = routetimes[i,j]
  }
}
g.mro.ts = ggplot(rtdf) + geom_segment(aes(x=PrevTime,xend=Time,y=Step-1,yend=Step,color=factor(Index)), alpha=0.5) +
  scale_x_continuous(trans="log") + theme_light()

sf = 2
png("plot-science-mro.png", width=1200*sf, height=1200*sf, res=72*sf)
grid.arrange(g.mro.cube, g.mro.summary, g.mro.motifs, g.mro.ts, nrow=2)
dev.off()

