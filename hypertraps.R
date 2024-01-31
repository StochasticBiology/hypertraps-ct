require(Rcpp)
require(ggplot2)
require(ggpubr)
require(ggraph)
require(ggwordcloud)
require(igraph)
require(stringr)
require(stringdist)

DecToBin <- function(x, len) {
  s = c()
  for(j in (len-1):0)
  {
    if(x >= 2**j) { s=c(s,1); x = x-2**j } else { s=c(s,0)}
  }
  return(paste(s, collapse=""))
}

BinToDec <- function(state) {
  this.ref = 0
  for(j in 1:length(state)) {
    this.ref = this.ref + state[j]*(2**(length(state)-j))
  }
  return(this.ref)
}

plotHypercube.lik.trace = function(my.post) {
  ### likelihood traces
  return(ggplot(my.post$lik.traces) + geom_line(aes(x=Step, y=LogLikelihood1)) +
           geom_line(aes(x=Step, y=LogLikelihood2)) + theme_light() )
}

plotHypercube.bubbles = function(my.post, reorder=FALSE, transpose=FALSE) {
  toplot = my.post$bubbles
  if(reorder == TRUE) {
    toplot$Name = factor(toplot$Name, levels=unique(toplot$Name))
  }
  if(transpose == TRUE) {
    toplot$x = toplot$Name
    toplot$y = toplot$Time
  } else {
    toplot$x = toplot$Time
    toplot$y = toplot$Name
  }
  this.plot = ggplot(toplot, aes(x=x, y=y, size=Probability)) + geom_point() +
    theme_light() 
  if(transpose == TRUE){
    return(this.plot + theme(axis.text.x = element_text(angle=90)) +
             xlab("") + ylab("Ordinal time"))
  } else {
    return(this.plot + xlab("Ordinal time") + ylab(""))
  }
}

plotHypercube.graph = function(my.post, thresh = 0.05, node.labels = TRUE,
                               node.label.size = 2) {
  ### produce hypercube subgraph
  bigL = my.post$L
  trans.p = my.post$dynamics$trans[my.post$dynamics$trans$Flux > thresh,]
  trans.g = graph_from_data_frame(trans.p)
  bs = unlist(lapply(as.numeric(V(trans.g)$name), DecToBin, len=bigL))
  V(trans.g)$binname = bs
  layers = str_count(bs, "1")
  this.plot =  ggraph(trans.g, layout="sugiyama", layers=layers) + 
    geom_edge_link(aes(edge_width=Flux, edge_alpha=Flux)) + 
    scale_edge_width(limits=c(0,NA)) + scale_edge_alpha(limits=c(0,NA)) +
    theme_graph(base_family="sans") #aes(label=bs)) + theme_graph() 
  if(node.labels == TRUE) {
    this.plot = this.plot + geom_node_point() + geom_node_label(aes(label=binname),
                                                                size = node.label.size) 
  }
  return(this.plot)
}

plotHypercube.sampledgraph = function(my.post, max.samps = 1000, thresh = 0.05, node.labels = TRUE, node.label.size = 2) {
  edge.from = edge.to = c()
  bigL = my.post$L
  nsamps = min(max.samps, nrow(my.post$routes))
  for(i in 1:nsamps) {
    state = 0
    for(j in 1:ncol(my.post$routes)) {
      edge.from = c(edge.from, state)
      state = state + 2**(my.post$L-my.post$routes[i,j]-1)
      edge.to = c(edge.to, state)
    }
  }
  df = data.frame(From=edge.from, To=edge.to)
  dfu = unique(df)
  dfu$Flux = 0
  for(i in 1:nrow(dfu)) {
    dfu$Flux[i] = length(which(df$From==dfu$From[i] & df$To==dfu$To[i]))
  }
  dfu = dfu[dfu$Flux > thresh*nsamps,]
  trans.g = graph_from_data_frame(dfu)
  bs = unlist(lapply(as.numeric(V(trans.g)$name), DecToBin, len=bigL))
  #bs = unlist(lapply(as.numeric(as.vector(V(trans.g))), DecToBin, len=bigL))
  V(trans.g)$binname = bs
  layers = str_count(bs, "1")
  this.plot = ggraph(trans.g, layout="sugiyama", layers=layers) + geom_edge_link(aes(edge_width=Flux, edge_alpha=Flux)) + 
    scale_edge_width(limits=c(0,NA)) + scale_edge_alpha(limits=c(0,NA)) +
    theme_graph(base_family="sans") #aes(label=bs)) + theme_graph() 
  if(node.labels == TRUE) {
    this.plot = this.plot + geom_node_point() + geom_node_label(aes(label=binname),size=node.label.size) 
  }
  return(this.plot)
}

plotHypercube.sampledgraph2 = function(my.post, max.samps = 1000, thresh = 0.05, 
                                       node.labels = TRUE, use.arc = TRUE, no.times = FALSE, 
                                       edge.label.size = 2, edge.label.angle = "across",
                                       edge.label.colour = "#000000",
                                       featurenames = c(""), truncate = -1,
                                       node.label.size = 2) {
  edge.from = edge.to = edge.time = edge.change = c()
  bigL = my.post$L
  if(truncate == -1 | truncate > bigL) { truncate = bigL }
  nsamps = min(max.samps, nrow(my.post$routes))
  for(i in 1:nsamps) {
    state = paste0(rep("0", bigL), collapse = "")
    for(j in 1:truncate) {
      edge.from = c(edge.from, state)
      locus = my.post$routes[i,j]+1
      substr(state, locus, locus) <- "1"
      edge.to = c(edge.to, state)
      edge.change = c(edge.change, my.post$routes[i,j])
      edge.time = c(edge.time, my.post$times[i,j])
    }
  }
  df = data.frame(From=edge.from, To=edge.to, Change=edge.change, Time=edge.time)
  dfu = unique(df)
  if(length(featurenames) > 1) {
    dfu$Change = featurenames[dfu$Change+1]
  }
  dfu$Flux = dfu$MeanT = dfu$SDT = NA
  for(i in 1:nrow(dfu)) {
    this.set = which(df$From==dfu$From[i] & df$To==dfu$To[i])
    dfu$Flux[i] = length(this.set)
    dfu$MeanT[i] = mean(df$Time[this.set])
    if(length(this.set) > 1) {
      dfu$SDT[i] = sd(df$Time[this.set])
    }
    if(no.times == TRUE) {
      dfu$label[i] = paste(c("+", dfu$Change[i]), collapse="")
    } else {
      dfu$label[i] = paste(c("+", dfu$Change[i], ": ", signif(dfu$MeanT[i], digits=2), " +- ", signif(dfu$SDT[i], digits=2)), collapse="") 
    }
    
  }
  dfu$Flux = dfu$Flux / nsamps
  dfu = dfu[dfu$Flux > thresh,]
  trans.g = graph_from_data_frame(dfu)
  #bs = unlist(lapply(as.numeric(V(trans.g)$name), DecToBin, len=bigL))
  #bs = unlist(lapply(as.numeric(as.vector(V(trans.g))), DecToBin, len=bigL))
  bs = V(trans.g)$name
  V(trans.g)$binname = bs
  layers = str_count(bs, "1")
  
  if(use.arc == TRUE) {
    this.plot=  ggraph(trans.g, layout="sugiyama", layers=layers) + 
      geom_edge_arc(aes(edge_width=Flux, edge_alpha=Flux, label=label), 
                    label_size = edge.label.size, label_colour=edge.label.colour, color="#AAAAFF",
                    label_parse = TRUE, angle_calc = edge.label.angle) + 
      scale_edge_width(limits=c(0,NA)) + scale_edge_alpha(limits=c(0,NA)) +
      theme_graph(base_family="sans")
  } else {
    this.plot=  ggraph(trans.g, layout="sugiyama", layers=layers) + 
      geom_edge_link(aes(edge_width=Flux, edge_alpha=Flux, label=label), 
                     label_size = edge.label.size, label_colour=edge.label.colour, color="#AAAAFF",
                     label_parse = TRUE, angle_calc = edge.label.angle) + 
      scale_edge_width(limits=c(0,NA)) + scale_edge_alpha(limits=c(0,NA)) +
      theme_graph(base_family="sans")
  }
  if(node.labels == TRUE) {
    this.plot = this.plot + geom_node_point() + geom_node_label(aes(label=binname),size=node.label.size) 
  }
  return(this.plot)
}


plotHypercube.timehists = function(my.post, t.thresh = 20, featurenames = c(""), log.time = TRUE) {
  thdfp = data.frame()
  if(length(featurenames) > 1) {
    my.post$timehists$feature.label = featurenames[my.post$timehists$OriginalIndex+1]
  } else {
    my.post$timehists$feature.label = my.post$timehists$OriginalIndex
  }
  for(i in c(3,4)) {
    for(j in unique(my.post$timehists$feature.label)) {
      sub = my.post$timehists[my.post$timehists$feature.label == j & my.post$timehists$Time < t.thresh,]
      sub1 = my.post$timehists[my.post$timehists$feature.label == j & my.post$timehists$Time >= t.thresh,]
      thdfp = rbind(thdfp, sub)
      thdfp = rbind(thdfp, data.frame(OriginalIndex = j, feature.label=j, Time=t.thresh, Probability=sum(sub1$Probability)))
    }
  }

  if (log.time) {
    g.thist = ggplot(thdfp[thdfp$Time < t.thresh,],
                     aes(x=log(Time+1), y=Probability)) +
      labs(x = "log(t+1)", y="Probability")
  } else {
    g.thist = ggplot(thdfp[thdfp$Time < t.thresh,],
                     aes(x=Time, y=Probability))
    + labs(x = "t", y="Probability")
  }
  g.thist = g.thist + 
    #geom_col(position="dodge") + xlim(-0.1,thresh+0.5) + facet_wrap(~OriginalIndex, ncol=2, scales="free") +
    geom_line() + xlim(-0.1,log(t.thresh+1)) + facet_wrap(~feature.label, nrow=2) +
    theme_light() #+ scale_x_continuous(trans="log10")

  g.thist2 = ggplot(thdfp[thdfp$Time == t.thresh,], aes(x=feature.label, y=Probability)) + 
    #geom_col(position="dodge") + xlim(-0.1,thresh+0.5) + facet_wrap(~OriginalIndex, ncol=2, scales="free") +
    geom_col(position="dodge") + 
    theme_light() + labs(x = "Feature", y="Probability") #+ scale_x_continuous(trans="log10")
  
  return(ggarrange(g.thist, g.thist2))
}

plotHypercube.regularisation = function(my.post) {
  return(ggplot(my.post$regularisation$reg.process, 
                aes(x=nparam, y=AIC)) + geom_point() + 
           labs(x = "Number of non-zero parameters", y="AIC") + theme_light() )
}

plotHypercube.motifs = function(my.post, featurenames = c("")) {
  # motif plot
  if(length(featurenames) > 1) {
    labels = featurenames
  } else {
    labels = 1:my.post$L
  }
  rdf = data.frame()
  for(j in 1:ncol(my.post$routes)) {
    startprob = 0
    for(i in 0:max(my.post$routes)) {
      thisprob = length(which(my.post$routes[,j]==i))/nrow(my.post$routes)
      rdf = rbind(rdf, data.frame(Index=i, Label=labels[i+1], Time=j, Start=startprob, End=startprob+thisprob, Probability=thisprob))
      startprob = startprob+thisprob
    }
  }
  return(ggplot(rdf) + geom_rect(aes(xmin=Time-0.5,xmax=Time+0.5,ymin=Start,ymax=End,fill=factor(Label))) +
           geom_text(aes(x=Time,y=(Start+End)/2,label=Label), color="#FFFFFF") + ylab("Probability") + 
           scale_fill_brewer(palette = "PuRd") + theme_light())
}

plotHypercube.timeseries = function(my.post, log.axis = TRUE, featurenames=c("")) {
  # time series illustration
  if(length(featurenames) > 1) {
    labels = featurenames
  } else {
    labels = 1:my.post$L
  }
  rtdf = data.frame()
  for(i in 1:(min(nrow(my.post$routes),1000))) {
    prevtime = 0
    for(j in 1:ncol(my.post$routes)) {
      rtdf = rbind(rtdf, data.frame(Run=i, Step=j, Label=labels[my.post$routes[i,j]+1], Index=my.post$routes[i,j], PrevTime=prevtime, Time=my.post$times[i,j]))
      prevtime = my.post$times[i,j]
    }
  }
  if(log.axis == TRUE) {
    return( ggplot(rtdf) + geom_segment(aes(x=PrevTime,xend=Time,y=Step-1,yend=Step,color=factor(Label, levels=labels)), alpha=0.5) +
              scale_x_continuous(trans="log10") + scale_color_brewer(palette = "Spectral") +
              labs(x= "t", y="Number of features", color = "Feature") + theme_light())
  } else {
    return ( ggplot(rtdf) + geom_segment(aes(x=PrevTime,xend=Time,y=Step-1,yend=Step,color=factor(Label, levels=labels)), alpha=0.5) +
      scale_color_brewer(palette = "Spectral")+ 
        labs(x= "t", y="Number of features", color = "Feature") + theme_light() )
  }
}

plotHypercube.summary = function(my.post, f.thresh = 0.05, t.thresh = 20, continuous.time = TRUE) {
  if(continuous.time == TRUE) {
    return (ggarrange(plotHypercube.lik.trace(my.post),
                      plotHypercube.bubbles(my.post),
                      plotHypercube.sampledgraph2(my.post, thresh = f.thresh, use.arc=FALSE, edge.label.size=3) + 
                        theme(legend.position="none") + expand_limits(x = c(-1, 4)),
                      plotHypercube.timehists(my.post, t.thresh), nrow=2, ncol=2) )
  } else {
    return (ggarrange(plotHypercube.lik.trace(my.post),
                      plotHypercube.bubbles(my.post),
                      plotHypercube.sampledgraph2(my.post, thresh = f.thresh, use.arc=FALSE, edge.label.size=3, no.times = TRUE) + 
                        theme(legend.position="none") + expand_limits(x = c(-1, 4)) ) )
  }
}

# plot pairwise influences between features in the L^2 picture
plotHypercube.influences = function(my.post, 
                                    feature.names=c(""), 
                                    use.regularised = FALSE, 
                                    use.final = FALSE, 
                                    reorder = FALSE, 
                                    upper.right = FALSE) {
  if(my.post$model != 2) {
    stop("Influence plot currently only supported for model type 2 (pairwise influences)")
  }
  plot.df = data.frame()
  if(length(feature.names) > 1) {
    labels = feature.names
  } else {
    labels = 1:my.post$L
  }
  for(i in 1:my.post$L) {
    for(j in 1:my.post$L) {
      ref = (i-1)*my.post$L + j
      if(use.regularised == TRUE) {
        ref.mean = as.numeric(my.post$regularisation$best[ref])
        ref.sd = 0
      } else if(use.final == TRUE) {
        ref.mean = mean(my.post$posterior.samples[nrow(my.post$posterior.samples),ref])
        ref.sd = 0
      } else {
        ref.mean = mean(my.post$posterior.samples[,ref])
        ref.sd = sd(my.post$posterior.samples[,ref])
      }
      plot.df = rbind(plot.df, data.frame(x=i, y=j, mean=ref.mean, cv=abs(ref.sd/ref.mean)))
    }
  }  
  plot.df$xlab = labels[plot.df$x]
  plot.df$ylab = labels[plot.df$y]
  if(reorder == TRUE) {
    diag.means = plot.df[plot.df$x == plot.df$y,]
    my.order = order(diag.means$mean, decreasing=TRUE)
    plot.df$xlab = factor(plot.df$xlab, levels=labels[my.order])
    plot.df$ylab = factor(plot.df$ylab, levels=labels[my.order])
  } else {
    plot.df$xlab = factor(plot.df$xlab)
    plot.df$ylab = factor(plot.df$ylab)
  }
  if(upper.right == TRUE) {
    return(ggplot(plot.df, aes(x=xlab,y=factor(ylab, levels=rev(levels(plot.df$ylab))),fill=mean,alpha=cv)) + geom_tile() + 
             scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
             scale_alpha_continuous(range=c(1,0)) +
             theme_light() + xlab("Acquired trait") + ylab("(Influenced) rate") +
             theme(axis.text.x = element_text(angle=90)) ) #+
    #scale_x_continuous(breaks=1:my.post$L, labels=labels) +
    #scale_y_continuous(breaks=1:my.post$L, labels=labels))
  } else {
    return(ggplot(plot.df, aes(x=xlab,y=ylab,fill=mean,alpha=cv)) + geom_tile() + 
             scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
             scale_alpha_continuous(range=c(1,0)) +
             theme_light() + xlab("Acquired trait") + ylab("(Influenced) rate") +
             theme(axis.text.x = element_text(angle=90)) ) #+
    # scale_x_continuous(breaks=1:my.post$L, labels=labels) +
    # scale_y_continuous(breaks=1:my.post$L, labels=labels))
  }
  
}


mylabel = function(label, suffix) {
  return(paste(c(label, suffix), collapse=""))
}

readHyperinf = function(label, postlabel = "", fulloutput=FALSE, regularised = FALSE) {
  rL = list()
  rL$label = label
  rL$lik.traces = read.csv(mylabel(label, "-lik.csv"))
  rL$L = rL$lik.traces$L[1]
  rL$model = rL$lik.traces$model[1]
  rL$best = read.table(mylabel(label, "-best.txt"))
  rL$posterior.samples = read.table(mylabel(label, "-posterior.txt"))
  
  if(fulloutput == TRUE) {
    tmpL = list()
    tmpL$states = read.csv(mylabel(label, "-states.csv"))
    tmpL$trans = read.csv(mylabel(label, "-trans.csv"))
    rL$dynamics = tmpL
  }
  
  if(regularised == TRUE) {
    tmpL = list()
    tmpL$best = read.table(mylabel(label, "-regularised.txt"))
    tmpL$reg.process = read.csv(mylabel(label, "-regularising.csv"))
    rL$regularisation = tmpL
  }
  
  if(postlabel != "") {
    rL$bubbles = read.csv(mylabel(postlabel, "-bubbles.csv"))
    rL$timehists = read.csv(mylabel(postlabel, "-timehists.csv"))
    rL$routes = read.table(mylabel(postlabel, "-routes.txt"), sep=" ")
    rL$betas = read.table(mylabel(postlabel, "-betas.txt"), sep=" ")
    rL$times = read.table(mylabel(postlabel, "-times.txt"), sep=" ") 
  }
  
  return(rL)
}

writeHyperinf = function(wL, label, postlabel = "", fulloutput=FALSE, regularised=FALSE) {
  write.table(t(wL$best), mylabel(label, "-best.txt"), row.names=FALSE, col.names=FALSE)
  write.table(wL$posterior.samples, mylabel(label, "-posterior.txt"), row.names=FALSE, col.names=FALSE)
  write.table(wL$lik.traces, mylabel(label, "-lik.csv"), row.names=FALSE, sep=",", quote=FALSE)
  
  if(fulloutput == TRUE) {
    write.table(wL$dynamics$states, mylabel(label, "-states.csv"), row.names=FALSE, sep=",", quote=FALSE)
    write.table(wL$dynamics$trans, mylabel(label, "-trans.csv"), row.names=FALSE, sep=",", quote=FALSE)
  }
  
  if(regularised == TRUE) {
    write.table(t(wL$regularisation$best), mylabel(label, "-regularised.txt"), row.names=FALSE, col.names=FALSE)
    write.table(wL$regularisation$reg.process, mylabel(label, "-regularising.csv"), row.names=FALSE, sep=",", quote=FALSE)
  }
  
  if(postlabel != "") {
    write.table(wL$bubbles, mylabel(postlabel, "-bubbles.csv"), row.names=FALSE, sep=",", quote=FALSE)
    write.table(wL$timehists, mylabel(postlabel, "-timehists.csv"), row.names=FALSE, sep=",", quote=FALSE)
    write.table(wL$routes, mylabel(postlabel, "-routes.txt"), row.names=FALSE, col.names = FALSE, sep=" ", quote=FALSE)
    write.table(wL$betas, mylabel(postlabel, "-betas.txt"), row.names=FALSE, col.names = FALSE, sep=" ", quote=FALSE)
    write.table(wL$times, mylabel(postlabel, "-times.txt"), row.names=FALSE, col.names = FALSE, sep=" ", quote=FALSE)
  }
}

pullFeatureLabels = function(my.post) {
  sub = my.post$bubbles[1:my.post$L,]
  return(as.vector(sub$Name[order(sub$OriginalIndex)]))
}

# construct the (probability-weighted) q-gram distance
qgramdist = function(my.post.1, my.post.2) {
  # pull routes and probabilities for first cube
  routes = table(apply(my.post.1$routes, 1, paste, collapse=""))
  L = ncol(my.post.1$routes)
  route.set = rownames(routes)
  route.probs = as.numeric(routes)/sum(as.numeric(routes))
  df.1 = data.frame()
  # loop through q-gram length
  for(i in 2:L) {
    # loop through routes found on the cube
    for(j in 1:length(route.set)) {
      # get q-grams from this route
      qgramset = qgrams(route.set[j], q=i)
      # sloppy. if this q-gram exists in our set, increase its score, otherwise add it
      for(k in 1:ncol(qgramset)) {
        ref = which(df.1$gram == colnames(qgramset)[k])
        if(length(ref) == 0) {
          df.1 = rbind(df.1, data.frame(gram = colnames(qgramset)[k], prob.2 = 0, prob.1 = qgramset[1,k]*route.probs[j]))
        } else {
          df.1$prob.1[ref] =  df.1$prob.1[ref] + qgramset[1,k]*route.probs[j]
        }
      }
    }
  }
  
  # now pull the second cube
  routes = table(apply(my.post.2$routes, 1, paste, collapse=""))
  route.set = rownames(routes)
  route.probs = as.numeric(routes)/sum(as.numeric(routes))
  # same loop as above
  for(i in 2:L) {
    for(j in 1:length(route.set)) {
      qgramset = qgrams(route.set[j], q=i)
      for(k in 1:ncol(qgramset)) {
        ref = which(df.1$gram == colnames(qgramset)[k])
        if(length(ref) == 0) {
          df.1 = rbind(df.1, data.frame(gram = colnames(qgramset)[k], prob.1 = 0, prob.2 = qgramset[1,k]*route.probs[j]))
        } else {
          df.1$prob.2[ref] = df.1$prob.2[ref] + qgramset[1,k]*route.probs[j]
        }
      }
    }
  }
  
  # build a named list of q-gram scores and final value (for debugging/exploration)
  returnlist = list()
  returnlist$df = df.1
  returnlist$val = sum(abs(df.1$prob.1-df.1$prob.2))
  return(returnlist)
}

# predict the next evolutionary step from a given state
# only works so far if the posterior structure has transition dynamics information
predictNextStep = function(my.post, state) {
  # get this state reference and look up exit routes
  this.ref = BinToDec(state)
  out.edges = my.post$dynamics$trans[my.post$dynamics$trans$From==this.ref,]
  out.probs = out.edges$Probability
  predictions = data.frame(states = unlist(lapply(out.edges$To, DecToBin, my.post$L)),
                           probs=out.probs)
  return(predictions)
}

# get the representation of different hypercube levels in a dataset
dataLevels = function(data.mat) {
  data.mat[data.mat==2] = 0
  counts = rowSums(data.mat)
  counts = as.numeric(table(counts))/length(counts)
  return(data.frame(level=1:length(counts), prop=counts))
}

# predict unobserved values in a given observation
# only works so far if the posterior structure has transition dynamics information
predictHiddenVals = function(my.post, state, level.weight=1) {
  # assign uniform weights to levels on the hypercube if not specified
  if(length(level.weight)==1) {
    level.weight = rep(1, my.post$L+1)
  }
  n1s = length(which(state==1))
  if(n1s > 0) {
    level.weight[1:n1s] = 0
  }
  n0s = length(which(state==0))
  if(n0s > 0) {
    level.weight[(my.post$L+2-n0s):(my.post$L+1)] = 0
  }
  # normalise level weights
  level.weight = level.weight/sum(level.weight)
  # get the unobserved indices and construct all binary combinations that could fill them
  hidden.set = which(state == 2)
  hidden.options = expand.grid(rep(list(0:1), length(hidden.set)))
  
  # initialise results
  res.df = data.frame()
  if(nrow(hidden.options) > 0) {
    # loop through each possible binary combination
    for(i in 1:nrow(hidden.options)) {
      # get reference to this particular state
      tmpstate = state
      tmpstate[hidden.set] = as.numeric(hidden.options[i,])
      ref = BinToDec(tmpstate)
      # pull this state's probability from learned hypercube
      raw.prob = my.post$dynamics$states$Probability[my.post$dynamics$states$State==ref]
      res.df = rbind(res.df, data.frame(state=paste(tmpstate,collapse=""), 
                                        level=sum(tmpstate),
                                        raw.prob=raw.prob))
    }
    # normalise probabilities across all options per level, and across all weighted levels
    res.df$prob = res.df$level.prob =  0
    for(i in n1s:(length(state)-n0s)) {
      res.df$level.prob[res.df$level==i] = res.df$raw.prob[res.df$level==i]/sum(res.df$raw.prob[res.df$level==i])
      res.df$prob[res.df$level==i] = res.df$raw.prob[res.df$level==i] * level.weight[i+1]
    }
    res.df$prob = res.df$prob/sum(res.df$prob)
    
    # produce parallel output describing aggregated 0/1 probabilities for each locus
    hidden.options.probs = cbind(hidden.options, as.numeric(res.df$prob))
    locus.probs = data.frame()
    for(i in 1:length(hidden.set)) {
      locus = hidden.set[i]
      prob = sum(hidden.options.probs[which(hidden.options.probs[,i]==1),ncol(hidden.options.probs)])
      locus.probs = rbind(locus.probs, data.frame(locus=locus, prob=prob))
    }
  } else {
    tmpstate = state
    ref = BinToDec(tmpstate)
    # pull this state's probability from learned hypercube
    res.df = rbind(res.df, data.frame(state=paste(tmpstate,collapse=""), 
                                      level=sum(tmpstate),
                                      raw.prob=my.post$dynamics$states$Probability[my.post$dynamics$states$State==ref],
                                      level.prob=1,
                                      prob=1))
    locus.probs = data.frame()
  }
  
  output.list = list()
  output.list$state.probs = res.df
  output.list$locus.probs = locus.probs
  
  return(output.list)
}

plotHypercube.prediction = function(prediction, max.size = 30) {
  if(length(prediction$states) > 0) {
    g.1 = ggplot(prediction, aes(label=states, size=probs), angle=0) + 
      geom_text_wordcloud() + scale_size_area(max_size = max.size) +
      theme_minimal()
    g.2 = ggplot(prediction, aes(x=states, y=probs)) +
      labs(x = "State", y = "Probability") +
      geom_col() + theme_light()  
  } else {
    g.1 = ggplot(prediction$state.probs, aes(label=state, size=prob), angle=0) + 
      geom_text_wordcloud() + scale_size_area(max_size = max.size) +
      theme_minimal()
    g.2 = ggplot(prediction$locus.probs, aes(x=factor(locus), y=prob)) + 
      labs(x = "State", y = "Probability") +
      geom_col() + theme_light() + xlab("Probability feature is a 1")
  }
  return(ggarrange(g.1, g.2))
}

# return sampled state probabilities at a given time tau (only well-posed for the continuous time case)
prob.by.time = function(my.post, tau) {
  tset = my.post$times
  fset = my.post$routes
  for(i in 2:ncol(tset)) {
    tset[,i] = tset[,i] + tset[,i-1]
  }
  sset = fset*ifelse(tset<tau, 1, NA)
  states = rowSums(2**sset, na.rm=TRUE)
  binstates = unlist(lapply(states, DecToBin, len=my.post$L))
  freqs = as.data.frame(table(binstates))
  df = data.frame(State = freqs$binstates, Probability = freqs$Freq/sum(freqs$Freq))
  return(df)
}

sourceCpp("hypertraps-r.cpp")

## Note: the warning
## "Using the `size` aesthetic in this geom was deprecated in ggplot2 3.4.0."
## is a known issue: https://github.com/thomasp85/ggraph/issues/333
