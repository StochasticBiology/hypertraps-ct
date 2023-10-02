#### x range, scale and log type to be edited

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

trans.1 = read.csv("Data/tb-dt-1-trans.txt", sep=" ")
trans.s.1 = read.csv("Data/tb-dt-1-states.txt", sep=" ")
g.tb.cube = plot.hypercube3(trans.1, statesdf=trans.s.1, 
                node.labels = FALSE, seg.labels = TRUE, threshold=5e-2)

g.tb.summary = grid.arrange(g.tb.bubbles, g.tb.thist, g.tb.thist2, nrow=3)
sf = 2
png("plot-science-tb.png", width=600*sf, height=1200*sf, res=72*sf)
grid.arrange(g.tb.cube, g.tb.summary, nrow=2)
dev.off()

######## routes analysis TB
routes = read.table("Data/tb-ct-1-posterior.txt-routes.txt")
routetimes = read.table("Data/tb-ct-1-posterior.txt-times.txt")

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
ggplot(rdf) + geom_rect(aes(xmin=Time-0.5,xmax=Time+0.5,ymin=Start,ymax=End,fill=factor(Index))) +
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
ggplot(rtdf) + geom_segment(aes(x=PrevTime,xend=Time,y=Step-1,yend=Step,color=factor(Index)), alpha=0.5) +
  scale_x_continuous(trans="log") + theme_light()

thist = data.frame()
for(i in 0:9) {
  times = routetimes[routes==1]
  thist = rbind(thist, data.frame(Index=i, times=times))
  XXX
}

#### bubble and hypercube plots for MRO experiments

fname = c("mro-3", "mro-4", "mro-5", "mro-6")
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

trans.1 = read.csv("Data/mro-1-trans.txt", sep=" ")
trans.s.1 = read.csv("Data/mro-1-states.txt", sep=" ")
g.mro.cube = plot.hypercube3(trans.1, statesdf=trans.s.1, 
                node.labels = FALSE, seg.labels = TRUE, threshold=1e-5)

g.mro.summary = grid.arrange(g.mro.bubbles, g.mro.thist, g.mro.thist2, ncol=1)
sf = 2
png("plot-science-mro.png", width=600*sf, height=1200*sf, res=72*sf)
grid.arrange(g.mro.cube, g.mro.summary, nrow=2)
dev.off()



#### routes analysis MRO
routes = read.table("Data/mro-3-posterior.txt-routes.txt")
routetimes = read.table("Data/mro-3-posterior.txt-times.txt")

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
ggplot(rdf) + geom_rect(aes(xmin=Time-0.5,xmax=Time+0.5,ymin=Start,ymax=End,fill=factor(Index))) +
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
ggplot(rtdf) + geom_segment(aes(x=PrevTime,xend=Time,y=Step-1,yend=Step,color=factor(Index)), alpha=0.5) +
  scale_x_continuous(trans="log") + theme_light()


####### comparison to old code for debugging
########

thdf.old = read.csv("OldTH/synth-1-data.txt-posterior-1-1-3-3-0.txt.ctrec.process", sep=" ", header=FALSE)
colnames(thdf.old) = c("OriginalIndex", "Time", "Probability")
thdf.old$year="old"

thdf.new = read.csv("../hctdump-main-old/VerifyData/synth-1-data.txt-posterior-1-1-3-3-0.txt-timehists.csv")
#colnames(thdf) = c("OriginalIndex", "Time", "Probability")
thdf.new$year="new"

thdf=rbind(thdf.old,thdf.new)
ggplot() + geom_line(data=thdf, 
                     aes(x=Time, y=Probability, color=year)) +
  facet_wrap(~OriginalIndex) + xlim(0,3) +
  theme_light()

########

thdf.old = read.csv("OldTH/ng.2878-S2.txt-pruned.txt-labels.txt-data.txt-posterior-1-2-4-4-0.txt.ctrec.process", sep=" ", header=FALSE)
colnames(thdf.old) = c("OriginalIndex", "Time", "Probability")
thdf.old$year="old"

thdf.new = read.csv("Data/tb-ct-1-posterior.txt-timehists.csv")
#colnames(thdf) = c("OriginalIndex", "Time", "Probability")
thdf.new$year="new"

thdf = thdf.old
thdf=rbind(thdf.old,thdf.new)
ggplot() + geom_line(data=thdf, 
                     aes(x=Time, y=Probability, color=year)) +
  facet_wrap(~OriginalIndex) + xlim(0,10) +
  theme_light()

ggplot(thdf, aes(x=Time, y=Probability, fill=factor(OriginalIndex))) + 
  geom_col(position="dodge") + xlim(-0.1,1.05) + facet_wrap(~year) +
  theme_light() #+ scale_x_continuous(trans="log10")


g.tb.thist = ggplot(thdf, aes(x=Time, y=Probability, fill=factor(year))) + 
  geom_col(position="dodge") + xlim(-0.1,10.5) + facet_wrap(~OriginalIndex) +
  theme_light() #+ scale_x_continuous(trans="log10")
