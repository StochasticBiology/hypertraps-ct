library(Rcpp)
library(ggplot2)
library(ggpubr)
library(ggraph)
library(igraph)
library(stringr)

DecToBin <- function(x, len) {
  s = c()
  for(j in (len-1):0)
  {
    if(x >= 2**j) { s=c(s,1); x = x-2**j } else { s=c(s,0)}
  }
  return(paste(s, collapse=""))
}

plotHypercube = function(my.post, my.post.out, t.thresh = 20, f.thresh = 0.05) {
  ### likelihood traces
  g.lik.trace = ggplot(my.post$lik.traces) + geom_line(aes(x=sample.times, y=l.samples.1)) +
    geom_line(aes(x=sample.times, y=l.samples.2))
  
  ### bubble plot
  g.bubbles = ggplot(my.post.out$Bubbles, aes(x=Time, y=OriginalIndex, size=Probability)) +
    geom_point() 
  
  ### produce hypercube subgraph
  bigL = my.post$L
  trans.p = my.post$dynamics$Edges[my.post$dynamics$Edges$Flux > f.thresh,]
  trans.g = graph_from_data_frame(trans.p)
  bs = unlist(lapply(as.numeric(V(trans.g)$name), DecToBin, len=bigL))
  V(trans.g)$binname = bs
  layers = str_count(bs, "1")
  g.cube = ggraph(trans.g, layout="sugiyama", layers=layers) + geom_edge_link(aes(edge_width=Flux, edge_alpha=Flux)) + 
    geom_node_point() + geom_node_label(aes(label=binname),size=2) +
    scale_edge_width(limits=c(0,NA)) + scale_edge_alpha(limits=c(0,NA)) +
    theme_graph() #aes(label=bs)) + theme_graph() 
  
  thdfp = data.frame()
  for(i in c(3,4)) {
    for(j in unique(my.post.out$THist$OriginalIndex)) {
      sub = my.post.out$THist[my.post.out$THist$OriginalIndex == j & my.post.out$THist$Time < t.thresh,]
      sub1 = my.post.out$THist[my.post.out$THist$OriginalIndex == j & my.post.out$THist$Time >= t.thresh,]
      thdfp = rbind(thdfp, sub)
      thdfp = rbind(thdfp, data.frame(OriginalIndex=j, Time=t.thresh, Probability=sum(sub1$Probability)))
    }
  }
  
  g.thist = ggplot(thdfp[thdfp$Time < t.thresh,], aes(x=log(Time+1), y=Probability)) + 
    #geom_col(position="dodge") + xlim(-0.1,thresh+0.5) + facet_wrap(~OriginalIndex, ncol=2, scales="free") +
    geom_line() + xlim(-0.1,log(thresh+1)) + facet_wrap(~OriginalIndex, nrow=2) +
    theme_light() #+ scale_x_continuous(trans="log10")
  
  g.thist2 = ggplot(thdfp[thdfp$Time == t.thresh,], aes(x=OriginalIndex, y=Probability)) + 
    #geom_col(position="dodge") + xlim(-0.1,thresh+0.5) + facet_wrap(~OriginalIndex, ncol=2, scales="free") +
    geom_col(position="dodge") + 
    theme_light() #+ scale_x_continuous(trans="log10")
  
  return(ggarrange(g.lik.trace, g.bubbles, g.cube, g.thist, g.thist2))
}

sourceCpp("hypertraps-r.cpp")

test.mat = as.matrix(read.table("../VerifyData/synth-cross-samples-1.txt"))
my.post = HyperTraPS(test.mat, model = -1, regularise = 1, walkers_arg = 20)
plot(my.post$regularisation$reg.process$AIC)

test.mat = as.matrix(read.table("../VerifyData/synth-cross-samples-1.txt"))
test.start.times = rep(0.1,nrow(test.mat))
test.end.times = rep(0.2,nrow(test.mat))
my.post = HyperTraPS(test.mat, 
                     starttimes_arg = test.start.times, 
                     endtimes_arg = test.end.times, 
                     outputinput_arg = 1)
my.post.out = PosteriorAnalysis(my.post)

g.obj = plotHypercube(my.post, my.post.out)

test.mat = as.matrix(read.table("../Data/total-observations.txt-trans.txt"))
starts = test.mat[seq(from=1, to=nrow(test.mat), by=2),]
ends = test.mat[seq(from=2, to=nrow(test.mat), by=2),]
my.post = HyperTraPS(ends, initialstates_arg = starts, length_index_arg = 4, outputinput= 1) 
my.names = as.vector(read.table("../Data/tools-names.txt"))[[1]]
my.post.out = PosteriorAnalysis(my.post, featurenames_arg = my.names)

g.lik.trace = ggplot(my.post$lik.traces) + geom_line(aes(x=sample.times, y=l.samples.1)) +
  geom_line(aes(x=sample.times, y=l.samples.2))
g.bubbles = ggplot(my.post.out$Bubbles, aes(x=Time, y=factor(Name, levels=unique(my.post.out$Bubbles$Name)), size=Probability)) +
  geom_point() 
ggarrange(g.lik.trace, g.bubbles, nrow=2)
