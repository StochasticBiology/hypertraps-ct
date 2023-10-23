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

plotHypercube.lik.trace = function(my.post) {
  ### likelihood traces
  return(ggplot(my.post$lik.traces) + geom_line(aes(x=Step, y=LogLikelihood1)) +
    geom_line(aes(x=Step, y=LogLikelihood2)) + theme_light() )
}

plotHypercube.bubbles = function(my.post, reorder=FALSE) {
  if(reorder == TRUE) {
    toplot = my.post$bubbles
    toplot$Name = factor(toplot$Name, levels=unique(toplot$Name))
    return(ggplot(toplot, aes(x=Time, y=Name, size=Probability)) +
             geom_point() +theme_light() )  
  } else {
  ### bubble plot
  return(ggplot(my.post$bubbles, aes(x=Time, y=Name, size=Probability)) +
    geom_point() +theme_light() )
  }
}
  
plotHypercube.graph = function(my.post, f.thresh = 0.05) {
  ### produce hypercube subgraph
  bigL = my.post$L
  trans.p = my.post$dynamics$trans[my.post$dynamics$trans$Flux > f.thresh,]
  trans.g = graph_from_data_frame(trans.p)
  bs = unlist(lapply(as.numeric(V(trans.g)$name), DecToBin, len=bigL))
  V(trans.g)$binname = bs
  layers = str_count(bs, "1")
  return( ggraph(trans.g, layout="sugiyama", layers=layers) + geom_edge_link(aes(edge_width=Flux, edge_alpha=Flux)) + 
    geom_node_point() + geom_node_label(aes(label=binname),size=2) +
    scale_edge_width(limits=c(0,NA)) + scale_edge_alpha(limits=c(0,NA)) +
    theme_graph() #aes(label=bs)) + theme_graph() 
  )
}

plotHypercube.timehists = function(my.post, t.thresh = 20) {
  thdfp = data.frame()
  for(i in c(3,4)) {
    for(j in unique(my.post$timehists$OriginalIndex)) {
      sub = my.post$timehists[my.post$timehists$OriginalIndex == j & my.post$timehists$Time < t.thresh,]
      sub1 = my.post$timehists[my.post$timehists$OriginalIndex == j & my.post$timehists$Time >= t.thresh,]
      thdfp = rbind(thdfp, sub)
      thdfp = rbind(thdfp, data.frame(OriginalIndex=j, Time=t.thresh, Probability=sum(sub1$Probability)))
    }
  }
  
  g.thist = ggplot(thdfp[thdfp$Time < t.thresh,], aes(x=log(Time+1), y=Probability)) + 
    #geom_col(position="dodge") + xlim(-0.1,thresh+0.5) + facet_wrap(~OriginalIndex, ncol=2, scales="free") +
    geom_line() + xlim(-0.1,log(t.thresh+1)) + facet_wrap(~OriginalIndex, nrow=2) +
    theme_light() #+ scale_x_continuous(trans="log10")
  
  g.thist2 = ggplot(thdfp[thdfp$Time == t.thresh,], aes(x=OriginalIndex, y=Probability)) + 
    #geom_col(position="dodge") + xlim(-0.1,thresh+0.5) + facet_wrap(~OriginalIndex, ncol=2, scales="free") +
    geom_col(position="dodge") + 
    theme_light() #+ scale_x_continuous(trans="log10")
  
  return(ggarrange(g.thist, g.thist2))
}

plotHypercube.regularisation = function(my.post) {
  return(ggplot(my.post.regularise$regularisation$reg.process, 
         aes(x=params, y=AIC)) + geom_point() + theme_light() )
}

plotHypercube.motifs = function(my.post) {
# motif plot
rdf = data.frame()
for(j in 1:ncol(my.post$routes)) {
  startprob = 0
  for(i in 0:max(my.post$routes)) {
    thisprob = length(which(my.post$routes[,j]==i))/nrow(my.post$routes)
    rdf = rbind(rdf, data.frame(Index=i, Time=j, Start=startprob, End=startprob+thisprob, Probability=thisprob))
    startprob = startprob+thisprob
  }
}
return(ggplot(rdf) + geom_rect(aes(xmin=Time-0.5,xmax=Time+0.5,ymin=Start,ymax=End,fill=factor(Index))) +
  geom_text(aes(x=Time,y=(Start+End)/2,label=Index), color="#FFFFFF") + ylab("Probability") + theme_light())
}

plotHypercube.timeseries = function(my.post) {
# time series illustration
rtdf = data.frame()
for(i in 1:(min(nrow(my.post$routes),1000))) {
  prevtime = 0
  for(j in 1:ncol(my.post$routes)) {
    rtdf = rbind(rtdf, data.frame(Run=i, Step=j, Index=my.post$routes[i,j], PrevTime=prevtime, Time=my.post$times[i,j]))
    prevtime = my.post$times[i,j]
  }
}
return( ggplot(rtdf) + geom_segment(aes(x=PrevTime,xend=Time,y=Step-1,yend=Step,color=factor(Index)), alpha=0.5) +
  scale_x_continuous(trans="log") + theme_light())
}

plotHypercube.summary = function(my.post, f.thresh = 0.05, t.thresh = 20) {
  return (ggarrange(plotHypercube.lik.trace(my.post.r),
            plotHypercube.bubbles(my.post.r),
            plotHypercube.graph(my.post.r, f.thresh),
            plotHypercube.timehists(my.post.r, t.thresh), nrow=2, ncol=2) )
}

mylabel = function(label, suffix) {
  return(paste(c(label, suffix), collapse=""))
}

readHyperinf = function(label, L, postlabel = "", fulloutput=FALSE) {
  rL = list()
  rL$label = label
  rL$L = L
  rL$best = read.table(mylabel(label, "-best.txt"))
  rL$posterior.samples = read.table(mylabel(label, "-posterior.txt"))
  rL$lik.traces = read.csv(mylabel(label, "-lik.csv"))
  
  if(fulloutput == TRUE) {
  tmpL = list()
  tmpL$states = read.csv(mylabel(label, "-states.csv"))
  tmpL$trans = read.csv(mylabel(label, "-trans.csv"))
  rL$dynamics = tmpL
  }
  
  if(postlabel != "") {
  rL$bubbles = read.csv(mylabel(postlabel, "-bubbles.csv"))
  rL$timehists = read.csv(mylabel(postlabel, "-timehists.csv"))
  rL$routes = read.table(mylabel(postlabel, "-routes.txt"))
  rL$betas = read.table(mylabel(postlabel, "-betas.txt"))
  rL$times = read.table(mylabel(postlabel, "-times.txt")) 
  }
  
  return(rL)
}

writeHyperinf = function(wL, label, L, postlabel = "", fulloutput=FALSE) {
  write.table(t(wL$best), mylabel(label, "-best.txt"), row.names=FALSE, col.names=FALSE)
  write.table(wL$posterior.samples, mylabel(label, "-posterior.txt"), row.names=FALSE, col.names=FALSE)
  write.table(wL$lik.traces, mylabel(label, "-lik.csv"), row.names=FALSE, sep=",", quote=FALSE)
  
  if(fulloutput == TRUE) {
    write.table(wL$dynamics$states, mylabel(label, "-states.csv"), row.names=FALSE, sep=",", quote=FALSE)
    write.table(wL$dynamics$trans, mylabel(label, "-trans.csv"), row.names=FALSE, sep=",", quote=FALSE)
  }
  
  if(postlabel != "") {
    write.table(wL$bubbles, mylabel(postlabel, "-bubbles.csv"), row.names=FALSE, sep=",", quote=FALSE)
    write.table(wL$timehists, mylabel(postlabel, "-timehists.csv"), row.names=FALSE, sep=",", quote=FALSE)
    write.table(wL$routes, mylabel(postlabel, "-routes.txt"), row.names=FALSE, col.names = FALSE, sep=",", quote=FALSE)
    write.table(wL$betas, mylabel(postlabel, "-betas.txt"), row.names=FALSE, col.names = FALSE, sep=",", quote=FALSE)
    write.table(wL$times, mylabel(postlabel, "-times.txt"), row.names=FALSE, col.names = FALSE, sep=",", quote=FALSE)
  }
}

sourceCpp("hypertraps-r.cpp")

