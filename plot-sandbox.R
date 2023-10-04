library(ggplot2)
library(gridExtra)
source("plot-trans.R")
library(stringdist)
library(pheatmap)
library(igraph)
library(ggraph)
library(ggplotify)

#### various test bed experiments

for(expt in 1:4) {
  # cross
  if(expt == 1) {
    fname = c("test-cross-mod-1", "test-cross-mod-2", "test-cross-mod-3"); 
    oneshot="no"; bigL=5
  }
  # hi-order logic (normal inference)
  if(expt == 2) {
    fname = c("test-ho-mod--1", "test-ho-mod-2", "test-ho-mod-3"); 
    oneshot="no"; bigL=5
  }
  # hi-order logic (longer MCMC chain)
  if(expt == 3) {
    fname = c("test-ho-mod--1-l", "test-ho-mod-2-l", "test-ho-mod-3-l"); 
    oneshot="no"; bigL=5
  }
  # hi-order logic (simulated annealing)
  if(expt == 4) {
    fname = c("test-ho-mod--1-sa", "test-ho-mod-2-sa", "test-ho-mod-3-sa"); 
    oneshot="yes"; bigL=5
  }
  # hi-order logic (simulated annealing, regularised)
  if(expt == 5) {
    fname = c("test-ho-mod--1-sa", "test-ho-mod-2-sa", "test-ho-mod-3-sa"); 
    oneshot="regularised"; bigL=4
  }
  
  # initialise plot list
  g.bubbles = g.lik = g.cube = g.gcube = g.motif = g.pheatmap = list()
  # loop over files in output
  for(i in 1:length(fname)) {
    bdf = thdf = data.frame()
    # set up filenames for input
    if(oneshot=="regularised") {
      base.name = paste(c(fname[i], "-regularised"))
      post.name = paste(c(fname[i], "-regularised.txt"))
    } else if(oneshot == "yes") {
      base.name = fname[i]
      post.name = paste(c(fname[i], "-best.txt"))
    } else {
      base.name = fname[i]
      post.name = paste(c(fname[i], "-posterior.txt"))
    }
      lik.name = paste(c("VerifyData/", base.name, "-lik.txt"), collapse="")
      trans.name = paste(c("VerifyData/", base.name, "-trans.txt"), collapse="")
      states.name = paste(c("VerifyData/", base.name, "-states.txt"), collapse="")
      bubble.name = paste(c("VerifyData/", post.name, "-bubbles.csv"), collapse="")
      thist.name = paste(c("VerifyData/", post.name, "-timehists.csv"), collapse="")
      routes.name = paste(c("VerifyData/", post.name, "-routes.txt"), collapse="")
   
    # read likelihood trace and produce plot
    lik.df = read.csv(lik.name)
    g.lik[[i]] = ggplot(lik.df) + 
      geom_line(aes(x=Step,y=LogLikelihood1), color="#FF8888", alpha=0.3, linewidth=3) +
      geom_line(aes(x=Step,y=LogLikelihood2), color="#8888FF", alpha=0.3, linewidth=3) +
      theme_light()
    
    # read bubble output and produce plot
    tmpdf = read.csv(bubble.name)
    tmpdf$Expt=i
    bdf = rbind(bdf, tmpdf)
    g.bubbles[[i]] = ggplot(bdf, aes(x=Time, y=OriginalIndex, size=Probability)) +
      geom_point() +
      theme_light()
    
    tmpdf = read.csv(thist.name)
    tmpdf$Expt=i
    thdf = rbind(thdf, tmpdf)
    
    # read transition and state probabilities
    trans.1 = read.csv(trans.name, sep=" ")
    trans.1 = trans.1[!is.nan(trans.1$Probability),]
    trans.s.1 = read.csv(states.name, sep=" ")
    trans.s.1 = trans.s.1[!is.nan(trans.s.1$Probability),]
    # produce plot using custom code
    g.cube[[i]] = plot.hypercube3(trans.1, statesdf=trans.s.1, 
                                  node.labels = FALSE, seg.labels = TRUE, threshold=5e-2)
    
    # set up metadata for ggraph plot
    trans.1$Flux = trans.1$Probability*trans.s.1$Probability[trans.1$From+1]
    
    trans.p = trans.1[trans.1$Flux > 1e-2,]
    trans.g = graph_from_data_frame(trans.p)
    bs = unlist(lapply(as.numeric(V(trans.g)$name), DecToBin, len=bigL))
    V(trans.g)$binname = bs
    layers = str_count(bs, "1")
    g.gcube[[i]] = ggraph(trans.g, layout="sugiyama", layers=layers) + geom_edge_link(aes(edge_width=Flux, edge_alpha=Flux)) + 
      geom_node_point() + geom_node_label(aes(label=binname),size=2) +
      scale_edge_width(limits=c(0,NA)) + scale_edge_alpha(limits=c(0,NA)) +
      theme_graph() #aes(label=bs)) + theme_graph() 
    
    # read individual routes
    routes = read.table(routes.name)
    
    # process route probabilities and produce motif plot
    rdf = data.frame()
    for(j in 1:ncol(routes)) {
      startprob = 0
      for(k in 0:max(routes)) {
        thisprob = length(which(routes[,j]==k))/nrow(routes)
        rdf = rbind(rdf, data.frame(Index=k, Time=j, Start=startprob, End=startprob+thisprob, Probability=thisprob))
        startprob = startprob+thisprob
      }
    }
    g.motif[[i]] = ggplot(rdf) + geom_rect(aes(xmin=Time-0.5,xmax=Time+0.5,ymin=Start,ymax=End,fill=factor(Index))) +
      geom_text(aes(x=Time,y=(Start+End)/2,label=Index), color="#FFFFFF") + ylab("Probability") + theme_light()
    
    # cluster routes and produce heatmap
    routes.str = apply(routes, 1, paste, collapse="")
    routes.tab = table(routes.str)
    r.str.df = data.frame(routes=rownames(routes.tab), counts=as.vector(routes.tab))
    big.routes = r.str.df$routes[r.str.df$counts>2]
    routes.str.uniq = unique(big.routes)
    sdm = stringdistmatrix(routes.str.uniq, routes.str.uniq, method='jw',p=0.1)
    rownames(sdm) = routes.str.uniq
    colnames(sdm) = routes.str.uniq
    g.pheatmap[[i]] = as.ggplot(pheatmap(sdm))
    
  }
  
  # produce summary output
  out.name = paste(c("plot-sandbox-", expt, ".png"), collapse="")
  sf = 2
  png(out.name, width=2000*sf, height=800*sf, res=72*sf)
  grid.arrange(g.lik[[1]], g.bubbles[[1]], g.gcube[[1]], g.motif[[1]], g.pheatmap[[1]], 
               g.lik[[2]], g.bubbles[[2]], g.gcube[[2]], g.motif[[2]], g.pheatmap[[2]],
               g.lik[[3]], g.bubbles[[3]], g.gcube[[3]], g.motif[[3]], g.pheatmap[[3]], 
               #           g.lik[[4]], g.bubbles[[4]], g.gcube[[4]], g.motif[[4]], g.pheatmap[[4]], 
               nrow=3)
  dev.off()
}
