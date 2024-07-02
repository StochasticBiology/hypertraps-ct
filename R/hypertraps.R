#' Hypercubic inference in continuous time
#'
#' @docType _PACKAGE
#' @name hypertrapsct
#' @useDynLib hypertrapsct
#' @importFrom Rcpp evalCpp sourceCpp
NULL
#> NULL

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

#' Likelihood trace
#' 
#' Creates a plot
#' 
#' @param my.post A fitted hypercube returned from HyperTraPS
#' @return a ggplot
#' @export
#' @examples
#' observations <- matrix(c(0,0,0,
#'                          0,0,1,
#'                          0,1,1,
#'                          1,1,1), byrow=TRUE, ncol=3)
#' fitted.cube <- HyperTraPS(observations)
#' plotHypercube.lik.trace(fitted.cube)
plotHypercube.lik.trace = function(my.post) {
  ### likelihood traces
  this.plot = ggplot2::ggplot(my.post$lik.traces)  
  if("CurrentLogLikelihood" %in% colnames(my.post$lik.traces)) {
    this.plot = this.plot + ggplot2::geom_line(ggplot2::aes(x=Step, y=CurrentLogLikelihood), color="#FF0000") 
  }
  
  this.plot = this.plot + ggplot2::geom_line(ggplot2::aes(x=Step, y=LogLikelihood1), linetype="dotted") +
    ggplot2::geom_line(ggplot2::aes(x=Step, y=LogLikelihood2), linetype="dotted")
  return(this.plot + ggplot2::theme_light() )
}

#' Acquisition ordering
#' 
#' Creates a plot
#' 
#' @param my.post A fitted hypercube returned from HyperTraPS
#' @return a ggplot
#' @export
#' @examples
#' observations <- matrix(c(0,0,0,
#'                          0,0,1,
#'                          0,1,1,
#'                          1,1,1), byrow=TRUE, ncol=3)
#' fitted.cube <- HyperTraPS(observations)
#' plotHypercube.bubbles(fitted.cube)
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
  this.plot = ggplot2::ggplot(toplot, ggplot2::aes(x=toplot$x, y=toplot$y, size=toplot$Probability)) + ggplot2::geom_point() +
    ggplot2::theme_light() 
  if(transpose == TRUE){
    return(this.plot + ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90)) +
             ggplot2::xlab("") + ggplot2::ylab("Ordinal time"))
  } else {
    return(this.plot + ggplot2::xlab("Ordinal time") + ggplot2::ylab(""))
  }
}

#' Binned acquisition ordering
#' 
#' Creates a plot
#' 
#' @param my.post A fitted hypercube returned from HyperTraPS
#' @param reorder Logical. Defaults to FALSE.
#' @param transpose Logical. Defaults to FALSE.
#' @param bins Numeric. Defaults to 5.
#' @return a ggplot
#' @export
#' @examples
#' observations <- matrix(c(0,0,0,
#'                          0,0,1,
#'                          0,1,1,
#'                          1,1,1), byrow=TRUE, ncol=3)
#' fitted.cube <- HyperTraPS(observations)
#' plotHypercube.bubbles.coarse(fitted.cube)
plotHypercube.bubbles.coarse = function(my.post, reorder=FALSE, transpose=FALSE, bins=5) {
  nbub = data.frame()
  toplot = my.post$bubbles
  tmax = max(toplot$Time)
  for(this.name in unique(toplot$Name)) {
    for(j in 1:bins) {
      prob = sum(toplot$Probability[toplot$Name==this.name & floor(bins*toplot$Time/tmax)+1==j])
      nbub = rbind(nbub, data.frame(Name=this.name, Time=j, Probability=prob))
    }
  }
  toplot = nbub
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
  this.plot = ggplot2::ggplot(toplot, ggplot2::aes(x=toplot$x, y=toplot$y, size=toplot$Probability)) + ggplot2::geom_point() +
    ggplot2::theme_light() 
  if(transpose == TRUE){
    return(this.plot + ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90)) +
             ggplot2::xlab("") + ggplot2::ylab("Ordinal time"))
  } else {
    return(this.plot + ggplot2::xlab("Ordinal time") + ggplot2::ylab(""))
  }
}

#' Plot transition graph
#' 
#' Creates a plot
#' 
#' @param my.post A fitted hypercube returned from HyperTraPS
#' @param thresh A numeric from 0 to 1. Lower flux threshold to display. Defaults to 0.05.
#' @param node.labels Logical. Defaults to TUE
#' @param node.label.size A numeric. Defaults to 2
#' @param node.label.box A logical. Defaults to FALSE.
#' @return a ggplot
#' @export
#' @examples
#' observations <- matrix(c(0,0,0,
#'                          0,0,1,
#'                          0,1,1,
#'                          1,1,1), byrow=TRUE, ncol=3)
#' fitted.cube <- HyperTraPS(observations)
#' plotHypercube.graph(fitted.cube, node.labels = FALSE)
plotHypercube.graph = function(my.post, thresh = 0.05, 
                               node.labels = TRUE,
                               node.label.size = 2,
                               node.labels.box = FALSE) {
  ### produce hypercube subgraph
  bigL = my.post$L
  trans.p = my.post$dynamics$trans[my.post$dynamics$trans$Flux > thresh,]
  trans.g = igraph::graph_from_data_frame(trans.p)
  bs = unlist(lapply(as.numeric(igraph::V(trans.g)$name), DecToBin, len=bigL))
  igraph::V(trans.g)$binname = bs
  layers = stringr::str_count(bs, "1")
  this.plot =  ggraph::ggraph(trans.g, layout="sugiyama", layers=layers) + 
    ggraph::geom_edge_link(ggplot2::aes(edge_width=trans.g$Flux, edge_alpha=trans.g$Flux), color="#AAAAFF") + 
    ggraph::scale_edge_width(limits=c(0,NA)) + ggraph::scale_edge_alpha(limits=c(0,NA)) +
    ggraph::theme_graph(base_family="sans") #aes(label=bs)) + theme_graph() 
  if(node.labels == TRUE) {
    if(node.labels.box == TRUE) {
      this.plot = this.plot + ggraph::geom_node_label(ggplot2::aes(label=trans.g$binname),
                                                                size = node.label.size) 
    } else {
      this.plot = this.plot + ggraph::geom_node_text(ggplot2::aes(label=trans.g$binname),
                                              size = node.label.size) 
    } 
  }
  return(this.plot)
}

#' Transition graph plot
#' 
#' Creates a plot
#' 
#' @param my.post A fitted hypercube returned from HyperTraPS
#' @param max.samps A numeric. Number of trajectories to plot. Defaults to 1000.
#' @param thresh A numeric from 0 to 1. Lower flux threshold to display. Defaults to 0.05.
#' @param node.labels Logical. Whether to display node labels. Defaults to TRUE.
#' @param node.label.size A numeric. Defaults to 2
#' @return a ggplot
#' @export
#' @examples
#' observations <- matrix(c(0,0,0,
#'                          0,0,1,
#'                          0,1,1,
#'                          1,1,1), byrow=TRUE, ncol=3)
#' fitted.cube <- HyperTraPS(observations)
#' plotHypercube.sampledgraph(fitted.cube, node.labels = FALSE)
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
  trans.g = igraph::graph_from_data_frame(dfu)
  bs = unlist(lapply(as.numeric(igraph::V(trans.g)$name), DecToBin, len=bigL))
  #bs = unlist(lapply(as.numeric(as.vector(V(trans.g))), DecToBin, len=bigL))
  igraph::V(trans.g)$binname = bs
  layers = stringr::str_count(bs, "1")
  this.plot = ggraph::ggraph(trans.g, layout="sugiyama", layers=layers) + ggraph::geom_edge_link(ggplot2::aes(edge_width=Flux, edge_alpha=Flux)) + 
    ggraph::scale_edge_width(limits=c(0,NA)) + ggraph::scale_edge_alpha(limits=c(0,NA)) +
    ggraph::theme_graph(base_family="sans") #aes(label=bs)) + theme_graph() 
  if(node.labels == TRUE) {
    this.plot = this.plot + ggraph::geom_node_point() + ggraph::geom_node_label(ggplot2::aes(label=trans.g$binname),size=node.label.size) 
  }
  return(this.plot)
}

#' Transition graph plot
#' 
#' Creates a plot
#' 
#' @param my.post A fitted hypercube returned from HyperTraPS
#' @param max.samps A numeric. Number of trajectories to plot. Defaults to 1000.
#' @param thresh A numeric from 0 to 1. Lower flux threshold to display. Defaults to 0.05.
#' @param node.labels Logical. Whether to display node labels. Defaults to TRUE.
#' @param use.arc Logical. Whether to use curved edges. Defaults to FALSE.
#' @param no.times Logical. Whether to display timing on edges. Defaults to FALSE.
#' @param small.times Logical. Whether to make the timing information tiny.
#' @param times.offset A numeric vector with relative x and y coordinates. 
#'      Defaults to c(0.1, -0.1).
#' @param edge.label.size A numeric. Defaults to 2.
#' @param edge.label.angle A string. Unused
#' @param edge.label.colour A hex color code. Defaults to "#000000".
#' @param edge.check.overlap A logical. Passed to [ggraph::geom_edge_link()]. 
#'      If TRUE, text that overlaps previous text in the same layer will not be plotted. 
#'      Defaults to TRUE.
#' @param featurenames A character vector
#' @param truncate An optional integer. Limits the number of features plotted.
#' @param node.label.size A numeric. Defaults to 2
#' @param use.timediffs Logical. Defaults to TRUE.
#' @return a ggplot
#' @export
#' @examples
#' observations <- matrix(c(0,0,0,
#'                          0,0,1,
#'                          0,1,1,
#'                          1,1,1), byrow=TRUE, ncol=3)
#' fitted.cube <- HyperTraPS(observations)
#' plotHypercube.sampledgraph2(fitted.cube, 
#'                             no.times = TRUE, 
#'                             use.timediffs = FALSE, 
#'                             node.labels = FALSE)
plotHypercube.sampledgraph2 = function(my.post, max.samps = 1000, thresh = 0.05, 
                                       node.labels = TRUE, use.arc = FALSE, no.times = FALSE, 
                                       small.times = FALSE, times.offset = c(0.1,-0.1),
                                       edge.label.size = 2, edge.label.angle = "across",
                                       edge.label.colour = "#000000", edge.check.overlap = TRUE,
                                       featurenames = TRUE, truncate = -1,
                                       node.label.size = 2, use.timediffs = TRUE) {
  if(featurenames == TRUE) {
    featurenames = my.post$featurenames
  } else {
    featurenames = c("")
  }
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
      if(use.timediffs == TRUE) {
        edge.time = c(edge.time, my.post$timediffs[i,j])
      } else {
        edge.time = c(edge.time, my.post$times[i,j])
      }
    }
  }
  
  df = data.frame(From=edge.from, To=edge.to, Change=edge.change, Time=edge.time)
  dfu = unique(df[,1:3])
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
      dfu$tlabel[i] = paste(c(signif(dfu$MeanT[i], digits=2), " +- ", signif(dfu$SDT[i], digits=2)), collapse="") 
    } else {
      dfu$label[i] = paste(c("+", dfu$Change[i], ":\n", signif(dfu$MeanT[i], digits=2), " +- ", signif(dfu$SDT[i], digits=2)), collapse="") 
    }
    
  }
  dfu$Flux = dfu$Flux / nsamps
  dfu = dfu[dfu$Flux > thresh,]
  trans.g = igraph::graph_from_data_frame(dfu)
  #bs = unlist(lapply(as.numeric(V(trans.g)$name), DecToBin, len=bigL))
  #bs = unlist(lapply(as.numeric(as.vector(V(trans.g))), DecToBin, len=bigL))
  bs = igraph::V(trans.g)$name
  igraph::V(trans.g)$binname = bs
  layers = stringr::str_count(bs, "1")
  
  if(truncate > bigL/2) { 
    this.plot=  ggraph::ggraph(trans.g, layout="sugiyama", layers=layers) 
  } else {
    this.plot=  ggraph::ggraph(trans.g, layout="tree")
  }
  if(use.arc == TRUE) {
    this.plot= this.plot +
      ggraph::geom_edge_arc(aes(edge_width=Flux, edge_alpha=Flux, label=label, angle=45), 
                    label_size = edge.label.size, label_colour=edge.label.colour, color="#AAAAFF",
                    label_parse = TRUE, angle_calc = edge.label.angle, check_overlap = edge.check.overlap) + 
      ggraph::scale_edge_width(limits=c(0,NA)) + ggraph::scale_edge_alpha(limits=c(0,NA)) +
      ggraph::theme_graph(base_family="sans")
  } else {
    this.plot=  this.plot +
      ggraph::geom_edge_link(ggplot2::aes(edge_width=Flux, edge_alpha=Flux, label=label, angle=45), 
                     label_size = edge.label.size, label_colour=edge.label.colour, color="#AAAAFF",
                     label_parse = TRUE, angle_calc = edge.label.angle, check_overlap = edge.check.overlap) + 
      ggraph::scale_edge_width(limits=c(0,NA)) + ggraph::scale_edge_alpha(limits=c(0,NA)) +
      ggraph::theme_graph(base_family="sans")
  }
  if(small.times == TRUE) {
    this.plot = this.plot + ggraph::geom_edge_link(ggplot2::aes(edge_width=Flux, edge_alpha=Flux, label=tlabel, angle=45), 
                                           label_size = edge.label.size-1, label_colour=edge.label.colour, alpha = 0, color="#AAAAFF",
                                           label_parse = TRUE, angle_calc = edge.label.angle, check_overlap = TRUE,
                                           position = ggplot2::position_nudge(x=times.offset[1], y = times.offset[2])) 
  }
  if(node.labels == TRUE) {
    this.plot = this.plot + ggraph::geom_node_text(ggplot2::aes(label=igraph::V(trans.g)$binname),size=node.label.size) 
  }
  
  return(this.plot)
}


#' Timing histograms
#' 
#' Creates a plot
#' 
#' @param my.post A fitted hypercube returned from HyperTraPS
#' @param t.thresh Upper threshold of time to display. Defaults to 20.
#' @param featurenames A logical. Whether or not to use featurenames. Defaults to TRUE
#' @param featurenames A character vector.
#' @param log.time Whether or not to log-transform time. Defaults to TRUE.
#' @return a ggplot
#' @export
#' @examples
#' observations <- matrix(c(0,0,0,
#'                          0,0,1,
#'                          0,1,1,
#'                          1,1,1), byrow=TRUE, ncol=3)
#' fitted.cube <- HyperTraPS(observations)
#' plotHypercube.timehists(fitted.cube)
plotHypercube.timehists = function(my.post, t.thresh = 20, featurenames = TRUE, log.time = TRUE) {
  thdfp = data.frame()
  if(featurenames == TRUE) {
    featurenames = my.post$featurenames
  } else {
    featurenames = c("")
  }
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
    g.thist = ggplot2::ggplot(thdfp[thdfp$Time < t.thresh,],
                     ggplot2::aes(x=log(Time+1), y=Probability)) +
      ggplot2::labs(x = "log(t+1)", y="Probability")
  } else {
    g.thist = ggplot2::ggplot(thdfp[thdfp$Time < t.thresh,],
                     ggplot2::aes(x=Time, y=Probability)) +
      ggplot2::labs(x = "t", y="Probability")
  }
  g.thist = g.thist + 
    #geom_col(position="dodge") + xlim(-0.1,thresh+0.5) + facet_wrap(~OriginalIndex, ncol=2, scales="free") +
    ggplot2::geom_line() + ggplot2::xlim(-0.1,log(t.thresh+1)) + ggplot2::facet_wrap(~feature.label, nrow=2) +
    ggplot2::theme_light() #+ scale_x_continuous(trans="log10")
  
  g.thist2 = ggplot2::ggplot(thdfp[thdfp$Time == t.thresh,], ggplot2::aes(x=feature.label, y=Probability)) + 
    #geom_col(position="dodge") + xlim(-0.1,thresh+0.5) + facet_wrap(~OriginalIndex, ncol=2, scales="free") +
    ggplot2::geom_col(position="dodge") + 
    ggplot2::theme_light() + ggplot2::labs(x = "Feature", y="Probability") #+ scale_x_continuous(trans="log10")
  
  return(ggpubr::ggarrange(g.thist, g.thist2))
}

#' Regularisation info and AIC
#' 
#' Creates a plot
#' 
#' @param my.post A fitted hypercube returned from HyperTraPS
#' @return a ggplot
#' @export
#' @examples
#' observations <- matrix(c(0,0,0,
#'                          0,0,1,
#'                          0,1,1,
#'                          1,1,1), byrow=TRUE, ncol=3)
#' fitted.cube <- HyperTraPS(observations, regularise = TRUE)
#' plotHypercube.regularisation(fitted.cube)
plotHypercube.regularisation = function(my.post) {
  return(ggplot2::ggplot(my.post$regularisation$reg.process, 
                ggplot2::aes(x=nparam, y=AIC)) + ggplot2::geom_point() + 
           ggplot2::labs(x = "Number of non-zero parameters", y="AIC") + ggplot2::theme_light() )
}

#' Motifs
#' 
#' Summarize the inferred acquisition order. Analogous to [plotHypercube.bubbles()]
#' 
#' @param my.post A fitted hypercube returned from HyperTraPS
#' @param featurenames Optional character vector
#' @param label.size Numeric. Defaults to 3.
#' @param label.scheme String. Defaults to "full". Will remove the labels if any
#'  other value is passed.
#' @return a ggplot
#' @export
#' @examples
#' observations <- matrix(c(0,0,0,
#'                          0,0,1,
#'                          0,1,1,
#'                          1,1,1), byrow=TRUE, ncol=3)
#' fitted.cube <- HyperTraPS(observations)
#' plotHypercube.motifs(fitted.cube)
plotHypercube.motifs = function(my.post, 
                                featurenames = TRUE, 
                                label.size=3, 
                                label.scheme = "full") {
  # motif plot
  if(featurenames == TRUE) {
    featurenames = my.post$featurenames
  } else {
    featurenames = c("")
  }
  if(length(featurenames) > 1) {
    labels = featurenames
  } else {
    labels = 1:my.post$L
  }
  sample.set = c(2, my.post$L-1, round(my.post$L/2), round(my.post$L/4), round(3*my.post$L/4))
  rdf = data.frame()
  for(j in 1:ncol(my.post$routes)) {
    startprob = 0
    for(i in 0:max(my.post$routes)) {
      thisprob = length(which(my.post$routes[,j]==i))/nrow(my.post$routes)
      this.label = labels[i+1]
      if(label.scheme != "full" & !(j %in% sample.set)) { this.label = "" }
      rdf = rbind(rdf, data.frame(Index=i, TextLabel=this.label, Label=labels[i+1], Time=j, Start=startprob, End=startprob+thisprob, Probability=thisprob))
      
      startprob = startprob+thisprob
    }
  }
  return(ggplot2::ggplot(rdf) + ggplot2::geom_rect(ggplot2::aes(xmin=Time-0.5,xmax=Time+0.5,ymin=Start,ymax=End,fill=factor(Label))) +
           ggplot2::geom_text(ggplot2::aes(x=Time,y=(Start+End)/2,label=TextLabel), color="#FFFFFF", size=label.size) + 
           ggplot2::labs(x = "Ordering", y="Probability", fill="Feature") + 
           viridis::scale_fill_viridis(discrete = TRUE, option="inferno", begin=0.2, end=0.8) +
           ggplot2::theme_light())
}

#' Time series plot
#' 
#' Visualize the acquisition of features in continuous time.
#'
#' @param my.post A fitted hypercube returned from HyperTraPS.
#' @param log.time Whether or not to log-transform timings. Defaults to TRUE
#' @param featurenames A logical. Whether or not to use feature names. Defaults to TRUE
#' @return a ggplot
#' @export
#' @examples
#' observations <- matrix(c(0,0,0,
#'                          0,0,1,
#'                          0,1,1,
#'                          1,1,1), byrow=TRUE, ncol=3)
#' fitted.cube <- HyperTraPS(observations)
#' plotHypercube.timeseries(fitted.cube)
plotHypercube.timeseries = function(my.post, log.time = TRUE, featurenames=TRUE) {
  # time series illustration
  if(featurenames == TRUE) {
    featurenames = my.post$featurenames
  } else {
    featurenames = c("")
  }
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
  if(log.time == TRUE) {
    return( ggplot2::ggplot(rtdf) + ggplot2::geom_segment(ggplot2::aes(x=PrevTime,xend=Time,y=Step-1,yend=Step,color=factor(Label, levels=labels)), alpha=0.5) +
              ggplot2::scale_x_continuous(trans="log10") + ggplot2::scale_color_brewer(palette = "Spectral") +
              ggplot2::labs(x= "t", y="Number of features", color = "Feature") + ggplot2::theme_light())
  } else {
    return ( ggplot2::ggplot(rtdf) + ggplot2::geom_segment(ggplot2::aes(x=PrevTime,xend=Time,y=Step-1,yend=Step,color=factor(Label, levels=labels)), alpha=0.5) +
               ggplot2::scale_color_brewer(palette = "Spectral")+ 
               ggplot2::labs(x= "t", y="Number of features", color = "Feature") + ggplot2::theme_light() )
  }
}

#' Rapidly visualize key metrics from a HyperTraPS run
#' 
#' Creates a combined likelihood trace, bubble plot and sampledgraph.
#' 
#' @param my.post A fitted hypercube returned from [HyperTraPS()]
#' @param f.thresh Lower threshold passed to [plotHypercube.sampledgraph2()]
#' @param t.thresh A upper threshold for time passed to [plotHypercube.timehists()]
#' @param continuous.time Whether or not to include [plotHypercube.timehists()]
#' @return a ggplot
#' @export
#' @examples
#' ancestors   <- matrix(c(0,0,0,
#'                         0,0,1), ncol=3)
#' descendants <- matrix(c(0,1,1,
#'                         1,1,1), ncol=3)
#' start       <- c(0,1)
#' end         <- c(2,5)
#' featurenames = c("A", "B", "C")
#'
#' fitted.cube <- HyperTraPS(descendants, 
#'                           initialstates = ancestors,
#'                           starttimes = start,
#'                           endtimes = end,
#'                           featurenames = featurenames)
#' plotHypercube.summary(fitted.cube)
plotHypercube.summary = function(my.post, f.thresh = 0.05, t.thresh = 20, continuous.time = TRUE) {
  if(continuous.time == TRUE) {
    return (ggpubr::ggarrange(plotHypercube.lik.trace(my.post),
                      plotHypercube.bubbles(my.post),
                      plotHypercube.sampledgraph2(my.post, thresh = f.thresh, use.arc=FALSE, edge.label.size=3) + 
                        ggplot2::theme(legend.position="none") + ggplot2::expand_limits(x = c(-1, 4)),
                      plotHypercube.timehists(my.post, t.thresh), nrow=2, ncol=2) )
  } else {
    return (ggpubr::ggarrange(plotHypercube.lik.trace(my.post),
                      plotHypercube.bubbles(my.post),
                      plotHypercube.sampledgraph2(my.post, thresh = f.thresh, use.arc=FALSE, edge.label.size=3, no.times = TRUE) + 
                        ggplot2::theme(legend.position="none") + ggplot2::expand_limits(x = c(-1, 4)) ) )
  }
}

#' Visualize pairwise influences
#' 
#' Plot pairwise influences between features in the L^2 picture.
#' 
#' @param my.post A fitted hypercube returned from HyperTraPS
#' @param featurenames A character vector.
#' @param use.regularised Logical. Defaults to FALSE.
#' @param use.final Logical. Defaults to FALSE.
#' @param reorder Logical. Defaults to FALSE.
#' @param upper.right Logical. Defaults to FALSE.
#' @param cv.thresh Numeric. Upper threshold. Defaults to Inf.
#' @param red.green Logical. Adjust colors for red-green color-blindness. 
#'      Defaults to FALSE.
#' @return a ggplot
#' @export
#' @examples
#' observations <- matrix(c(0,0,0,
#'                          0,0,1,
#'                          0,1,1,
#'                          1,1,1), byrow=TRUE, ncol=3)
#' fitted.cube <- HyperTraPS(observations)
#' plotHypercube.bubbles(fitted.cube)
plotHypercube.influences = function(my.post, 
                                    featurenames=TRUE, 
                                    use.regularised = FALSE, 
                                    use.final = FALSE, 
                                    reorder = FALSE, 
                                    upper.right = FALSE,
                                    cv.thresh = Inf,
                                    red.green = FALSE) {
  if(my.post$model != 2) {
    stop("Influence plot currently only supported for model type 2 (pairwise influences)")
  }
  if(featurenames == TRUE) {
    featurenames = my.post$featurenames
  } else {
    featurenames = c("")
  }
  plot.df = data.frame()
  if(length(featurenames) > 1) {
    labels = featurenames
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
  plot.df$cv[is.na(plot.df$cv)] = 0
  plot.df$precision = 1-plot.df$cv
  plot.df$precision[plot.df$precision < 0] = 0
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
  plot.df = plot.df[plot.df$cv < cv.thresh | plot.df$x==plot.df$y,]
  
  if(red.green == TRUE) {
    low.col = "#FF8888"
    high.col = "#338833"
  } else {
    low.col = "#FF8888"
    high.col = "#8888FF"
  }
  if(upper.right == TRUE) {
    if(cv.thresh == Inf) {
      this.plot = ggplot2::ggplot(plot.df, ggplot2::aes(x=plot.df$xlab,y=factor(plot.df$ylab, levels=rev(levels(plot.df$ylab))),fill=mean,
                                      alpha=precision))
    } else {
      this.plot = ggplot2::ggplot(plot.df, ggplot2::aes(x=plot.df$xlab,y=factor(plot.df$ylab, levels=rev(levels(plot.df$ylab))),fill=mean))
    } } else {
      if(cv.thresh == Inf) {
        this.plot = ggplot2::ggplot(plot.df, ggplot2::aes(x=plot.df$xlab,y=plot.df$ylab,fill=plot.df$mean,alpha=plot.df$precision)) 
      } else {
      this.plot = ggplot2::ggplot(plot.df, ggplot2::aes(x=plot.df$xlab,y=plot.df$ylab,fill=plot.df$mean)) 
    } }

   return(this.plot + ggplot2::geom_tile() +
             ggplot2::scale_fill_gradient2(low = low.col, mid = "white", high = high.col, midpoint = 0) +
             ggplot2::scale_alpha_continuous(range=c(0,1)) +
             ggplot2::theme_light() + 
             ggplot2::labs(x="Acquired trait", y="(Influenced) rate", fill="Posterior\nmean", alpha="Posterior\nprecision") +
             ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90)) ) #+
}


mylabel = function(label, suffix) {
  return(paste(c(label, suffix), collapse=""))
}

#' Read data written to file
#' 
#' `readHyperinf` reads csv-files created by [writeHyperinf()].
#' 
#' @param label A string that describes the path and prefix that will be used.
#' @param postlabel A string that describes the path and prefix that will be 
#'      used to save posterior analysis information to file.
#' @param fulloutput whether states and transitions should be read.
#'      Defaults to FALSE.
#' @param regularised Whether regularisation data should be read.
#'      Defaults to FALSE.
#' @return A list of loaded data.
#' @export
#' @examples
#' observations <- matrix(c(0,0,0,
#'                          0,0,1,
#'                          0,1,1,
#'                          1,1,1), byrow=TRUE, ncol=3)
#' fitted.cube <- HyperTraPS(observations)
#' plotHypercube.bubbles(fitted.cube)
readHyperinf = function(label, postlabel = "", fulloutput=FALSE, regularised = FALSE) {
  rL = list()
  rL$label = label
  rL$lik.traces = utils::read.csv(mylabel(label, "-lik.csv"))
  rL$L = rL$lik.traces$L[1]
  rL$model = rL$lik.traces$model[1]
  rL$best = utils::read.table(mylabel(label, "-best.txt"))
  rL$posterior.samples = utils::read.table(mylabel(label, "-posterior.txt"))
  
  if(fulloutput == TRUE) {
    tmpL = list()
    tmpL$states = utils::read.csv(mylabel(label, "-states.csv"))
    tmpL$trans = utils::read.csv(mylabel(label, "-trans.csv"))
    rL$dynamics = tmpL
  }
  
  if(regularised == TRUE) {
    tmpL = list()
    tmpL$best = utils::read.table(mylabel(label, "-regularised.txt"))
    tmpL$reg.process = utils::read.csv(mylabel(label, "-regularising.csv"))
    rL$regularisation = tmpL
  }
  
  if(postlabel != "") {
    rL$bubbles = utils::read.csv(mylabel(postlabel, "-bubbles.csv"))
    rL$timehists = utils::read.csv(mylabel(postlabel, "-timehists.csv"))
    rL$routes = utils::read.table(mylabel(postlabel, "-routes.txt"), sep=" ")
    rL$betas = utils::read.table(mylabel(postlabel, "-betas.txt"), sep=" ")
    rL$times = utils::read.table(mylabel(postlabel, "-times.txt"), sep=" ") 
    # older versions of code didn't output this, so catch no-file errors
    tryCatch( { rL$timediffs = utils::read.table(mylabel(postlabel, "-timediffs.txt"), sep=" ") }, 
              error=function(e) {
                cat("Didn't find time diff data\n")
              } )
  }
  
  return(rL)
}

#' Save a fitted hypercube
#' 
#' `writeHyperinf` saves output from [HyperTraPS()] to file.
#' 
#' @param wL A fitted hypercube returned from [HyperTraPS()]
#' @param label A string that describes the path and prefix that will be used.
#' @param postlabel A string that describes the path and prefix that will be 
#'      used to save posterior analysis information to file.
#' @param fulloutput whether states and transitions should be written as
#'      csv-files. Defaults to FALSE.
#' @param regularised Whether regularisation data should be written to file.
#'      Defaults to FALSE.
#' @export
#' @examples
#' observations <- matrix(c(0,0,0,
#'                          0,0,1,
#'                          0,1,1,
#'                          1,1,1), byrow=TRUE, ncol=3)
#' fitted.cube <- HyperTraPS(observations)
#' writeHyperinf(fitted.cube, "my_cube")
writeHyperinf = function(wL, label, postlabel = "", fulloutput=FALSE, regularised=FALSE) {
    utils::write.table(t(wL$best), mylabel(label, "-best.txt"), row.names=FALSE, col.names=FALSE)
  utils::write.table(wL$posterior.samples, mylabel(label, "-posterior.txt"), row.names=FALSE, col.names=FALSE)
  utils::write.table(wL$lik.traces, mylabel(label, "-lik.csv"), row.names=FALSE, sep=",", quote=FALSE)
  
  if(fulloutput == TRUE) {
      utils::write.table(wL$dynamics$states, mylabel(label, "-states.csv"), row.names=FALSE, sep=",", quote=FALSE)
    utils::write.table(wL$dynamics$trans, mylabel(label, "-trans.csv"), row.names=FALSE, sep=",", quote=FALSE)
  }
  
  if(regularised == TRUE) {
      utils::write.table(t(wL$regularisation$best), mylabel(label, "-regularised.txt"), row.names=FALSE, col.names=FALSE)
    utils::write.table(wL$regularisation$reg.process, mylabel(label, "-regularising.csv"), row.names=FALSE, sep=",", quote=FALSE)
  }
  
  if(postlabel != "") {
      utils::write.table(wL$bubbles, mylabel(postlabel, "-bubbles.csv"), row.names=FALSE, sep=",", quote=FALSE)
    utils::write.table(wL$timehists, mylabel(postlabel, "-timehists.csv"), row.names=FALSE, sep=",", quote=FALSE)
    utils::write.table(wL$routes, mylabel(postlabel, "-routes.txt"), row.names=FALSE, col.names = FALSE, sep=" ", quote=FALSE)
    utils::write.table(wL$betas, mylabel(postlabel, "-betas.txt"), row.names=FALSE, col.names = FALSE, sep=" ", quote=FALSE)
    utils::write.table(wL$times, mylabel(postlabel, "-times.txt"), row.names=FALSE, col.names = FALSE, sep=" ", quote=FALSE)
    utils::write.table(wL$timediffs, mylabel(postlabel, "-timediffs.txt"), row.names=FALSE, col.names = FALSE, sep=" ", quote=FALSE)
  }
}

pullFeatureLabels = function(my.post) {
  sub = my.post$bubbles[1:my.post$L,]
  return(as.vector(sub$Name[order(sub$OriginalIndex)]))
}

#' Calculate q-gram distance
#' 
#' Construct the (probability-weighted) q-gram distance between two hypercubes.
#' 
#' @param my.post.1 A fitted hypercube returned from [HyperTraPS()].
#' @param my.post.2 A fitted hypercube returned from [HyperTraPS()].
#' @return a list containing the distance between the two cubes.
#' @export
#' @examples
#' observations.1 <- matrix(c(0,0,0,
#'                            0,0,1,
#'                            0,1,1,
#'                            1,1,1), byrow=TRUE, ncol=3)
#' fitted.cube.1 <- HyperTraPS(observations.1)
#' observations.2 <- matrix(c(0,0,0,
#'                            1,0,0,
#'                            1,1,0,
#'                            1,1,1), byrow=TRUE, ncol=3)
#' fitted.cube.2 <- HyperTraPS(observations.2)
#' qgramdist(fitted.cube.1, fitted.cube.2)
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
      qgramset = stringdist::qgrams(route.set[j], q=i)
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
      qgramset = stringdist::qgrams(route.set[j], q=i)
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

#' Predict the next likely step
#' 
#' predict the next evolutionary step from a given state
#' only works so far if the posterior structure has transition dynamics information
#' 
#' @param my.post A fitted hypercube returned from HyperTraPS
#' @param state A numeric vector of 1s and 0s. 
#' @return A data frame of states and probabilities for the next step
#' @export
#' @examples
#' observations <- matrix(c(0,0,0,
#'                          0,0,1,
#'                          0,1,1,
#'                          1,1,1), byrow=TRUE, ncol=3)
#' fitted.cube <- HyperTraPS(observations)
#' predictNextStep(fitted.cube, c(0,0,1))
predictNextStep = function(my.post, state) {
  # get this state reference and look up exit routes
  this.ref = BinToDec(state)
  out.edges = my.post$dynamics$trans[my.post$dynamics$trans$From==this.ref,]
  out.probs = out.edges$Probability
  predictions = data.frame(states = unlist(lapply(out.edges$To, DecToBin, my.post$L)),
                           probs=out.probs)
  return(predictions)
}

#' Levels of a dataset
#' 
#' get the representation of different hypercube levels in a dataset
#' 
#' @param data.mat A numeric matrix
#' @return A data frame containing counts for each level of the hypercube
#' @export
#' @examples
#' observations <- matrix(c(0,0,0,
#'                          0,0,1,
#'                          0,1,1,
#'                          1,1,1), byrow=TRUE, ncol=3)
#' dataLevels(observations)
dataLevels = function(data.mat) {
  data.mat[data.mat==2] = 0
  counts = rowSums(data.mat)
  counts = as.numeric(table(counts))/length(counts)
  return(data.frame(level=1:length(counts), prop=counts))
}

#' Impute data
#' 
#' predict unobserved values in a given observation
#' only works so far if the posterior structure has transition dynamics information
#' 
#' @param my.post A fitted hypercube returned from HyperTraPS
#' @param state An integer vector detailing the state to impute. Hidden values 
#'      are denoted with 2.
#' @return A list containing two data frames: state.probs and locus.probs.
#' @export
#' @examples
#' observations <- matrix(c(0,0,0,
#'                          0,0,1,
#'                          0,1,1,
#'                          1,1,1), byrow=TRUE, ncol=3)
#' fitted.cube <- HyperTraPS(observations)
#' predictHiddenVals(fitted.cube, c(0,1,2))
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

#' Visualize probabilities at a given set of times
#' 
#' Creates a plot
#' plot genotype probabilities at a given set of times as a "motif" plot
#' 
#' @param my.post A fitted hypercube returned from HyperTraPS.
#' @param t.set A numeric vector of times to visualize.
#' @param thresh A numeric of the lower bound threshold for probabilities 
#'      to display. Defaults to 0.05.
#' @param label.size Defaults to 3.
#' @return a ggplot
#' @export
#' @examples
#' observations <- matrix(c(0,0,0,
#'                          0,0,1,
#'                          0,1,1,
#'                          1,1,1), byrow=TRUE, ncol=3)
#' fitted.cube <- HyperTraPS(observations)
#' plotHypercube.motifseries(fitted.cube)
plotHypercube.motifseries = function(my.post, t.set=0, thresh = 0.05, label.size=3) {
  df = data.frame()
  t.index = 1
  # build up dataframe with rectangle co-ordinates and labels for high-probability states
  for(this.t in t.set) {
    tmp.df = prob.by.time(my.post, this.t)
    tmp.df$t = this.t
    tmp.df$t.index = t.index
    tmp.df$s1 = tmp.df$s2 = 0
    tmp.df$label = ""
    tmp.df$s2[1] = tmp.df$Probability[1]
    if(tmp.df$Probability[1] > thresh) { tmp.df$label[1] = as.character(tmp.df$State[1]) }
    if(nrow(tmp.df) > 1) {
      for(i in 2:nrow(tmp.df)) {
        tmp.df$s1[i] = tmp.df$s2[i-1]
        tmp.df$s2[i] = tmp.df$s1[i]+tmp.df$Probability[i]
        if(tmp.df$Probability[i] > thresh) { tmp.df$label[i] = as.character(tmp.df$State[i]) }
      }
    }
    df = rbind(df, tmp.df)
    t.index = t.index + 1
  }
  return(
    ggplot2::ggplot(df) + ggplot2::geom_rect(ggplot2::aes(xmin=t.index-0.5,xmax=t.index+0.5,ymin=s1,ymax=s2,fill=factor(State)), color="#FFFFFF22") +
      ggplot2::geom_text(ggplot2::aes(x=t.index,y=(s1+s2)/2,label=label), color="#FFFFFF", size=label.size) + 
      ggplot2::labs(x = "Time", y="Probability", fill="State") + 
      viridis::scale_fill_viridis(discrete = TRUE, option="inferno", begin=0.2, end=0.8) +
      ggplot2::theme_light() + ggplot2::theme(legend.position = "none") +
      ggplot2::scale_x_continuous(breaks = 1:length(t.set), labels=as.character(t.set))
  )
}

#' Visualize pairwise influences
#' 
#' plot pairwise influences between features in the L^2 or L^3 picture
#' 
#' @param my.post A fitted hypercube returned from HyperTraPS.
#' @param featurenames A character vector of length L. Optional 
#' @param use.regularised A logical. Defaults to FALSE
#' @param use.final A logical. Defaults to FALSE.
#' @param thresh A numeric describing the threshold. Defaults to 0.05
#' @param cv.thresh A numeric describing the upper bound. Defaults to Inf
#' @param label.size A numeric describing label size in mm. Defaults to 2
#' @param red.green A logical.
#' @return a ggplot
#' @export
#' @examples
#' observations <- matrix(c(0,0,0,
#'                          0,0,1,
#'                          0,1,1,
#'                          1,1,1), byrow=TRUE, ncol=3)
#' fitted.cube <- HyperTraPS(observations)
#' plotHypercube.influencegraph(fitted.cube, featurenames=FALSE)
plotHypercube.influencegraph = function(my.post, 
                                        featurenames=TRUE, 
                                        use.regularised = FALSE, 
                                        use.final = FALSE,
                                        thresh=0.05,
                                        cv.thresh = Inf,
                                        label.size = 2,
                                        red.green = FALSE) {
  plot.df = data.frame()
  if(featurenames == TRUE) {
    featurenames = my.post$featurenames
  } else {
    featurenames = c("")
  }
  if(length(featurenames) > 1) {
    labels = featurenames
  } else {
    labels = 1:my.post$L
  }
  
  if(my.post$model == 2) {
    for(i in 1:my.post$L) {
      for(j in 1:my.post$L) {
        ref = (i-1)*my.post$L + (j-1) + 1
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
        if(i != j) {
          plot.df = rbind(plot.df, data.frame(x=labels[i], y=labels[j], mean=ref.mean, cv=abs(ref.sd/ref.mean)))
        }
      }
    }  
  } else if(my.post$model == 3) {
    for(i in 1:my.post$L) {
      for(j in 1:my.post$L) {
        for(k in 1:my.post$L) {
          ref = (j-1)*my.post$L*my.post$L + (i-1)*my.post$L + k
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
          if(i == j) { this.xlab = labels[i] }
          if(i != j) { this.xlab = paste0(labels[i], "+", labels[j]) }
          if(i != k & j != k) {
            plot.df = rbind(plot.df, data.frame(x=this.xlab, y=labels[k], mean=ref.mean, cv=abs(ref.sd/ref.mean)))
          }
        }
      }  
    }
    
  } else {
    stop("Influence plot currently only supported for model type 2 or 3 (pairwise/tripletwise influences)")
  }
  plot.df = plot.df[plot.df$cv < cv.thresh,]
  to.g.df = plot.df[abs(plot.df$mean)>thresh,1:3]
  colnames(to.g.df) = c("From", "To", "Weight")
  to.g.df$Direction = factor(as.character(sign(to.g.df$Weight)), levels=c("-1", "1"))
  g = igraph::graph_from_data_frame(to.g.df)
  this.plot = ggraph::ggraph(g, layout="kk")
  if(red.green == TRUE) {
    low.col = "#FF8888"
    high.col = "#338833"
  } else {
    low.col = "#FF8888"
    high.col = "#8888FF"
  }
  return(  this.plot + ggraph::geom_edge_arc(ggplot2::aes(colour=Direction, alpha=abs(to.g.df$Weight), width=abs(to.g.df$Weight)),
                                     strength=0.1, arrow=grid::arrow(length=grid::unit(0.2, "inches"), type="closed")) +
             ggraph::geom_node_label(ggplot2::aes(label=name), size=label.size) + ggplot2::theme_void() +
             ggplot2::labs(edge_width="Magnitude", edge_alpha="Magnitude", colour="Direction") +
             ggraph::scale_edge_colour_manual(values = setNames(c(low.col, high.col), factor(c("-1", "1"))))
  )
  
}

#' Visualize predicted hidden values
#' 
#' `plotHypercube.prediction` visualizes output from `predictHiddenVals`
#' 
#' @param prediction Created with `predictHiddenVals`
#' @param max.size Defaults to 30. Size of largest points
#' @return a ggplot
#' @export
#' @examples
#' observations <- matrix(c(0,0,0,
#'                          0,0,1,
#'                          0,1,1,
#'                          1,1,1), byrow=TRUE, ncol=3)
#' fitted.cube <- HyperTraPS(observations)
#' prediction <- predictHiddenVals(fitted.cube, c(0,2,1))
#' plotHypercube.prediction(prediction)
plotHypercube.prediction = function(prediction, max.size = 30) {
  if(length(prediction$states) > 0) {
    g.1 = ggplot2::ggplot(prediction, ggplot2::aes(label=states, size=probs), angle=0) + 
      ggwordcloud::geom_text_wordcloud() + ggplot2::scale_size_area(max_size = max.size) +
      ggplot2::theme_minimal()
    g.2 = ggplot2::ggplot(prediction, ggplot2::aes(x=states, y=probs)) +
      ggplot2::labs(x = "State", y = "Probability") +
      ggplot2::geom_col() + ggplot2::theme_light()  
  } else {
    g.1 = ggplot2::ggplot(prediction$state.probs, ggplot2::aes(label=state, size=prob), angle=0) + 
      ggwordcloud::geom_text_wordcloud() + ggplot2::scale_size_area(max_size = max.size) +
      ggplot2::theme_minimal()
    g.2 = ggplot2::ggplot(prediction$locus.probs, ggplot2::aes(x=factor(locus), y=prob)) + 
      ggplot2::labs(x = "State", y = "Probability") +
      ggplot2::geom_col() + ggplot2::theme_light() + ggplot2::labs(x="Feature",y="Probability feature is a 1")
  }
  return(ggpubr::ggarrange(g.1, g.2))
}

# return sampled state probabilities at a given time tau (only well-posed for the continuous time case)
prob.by.time = function(my.post, tau) {
  tset = my.post$times
  fset = my.post$routes
  for(i in 2:ncol(tset)) {
    tset[,i] = tset[,i] + tset[,i-1]
  }
  sset = fset*ifelse(tset<tau, 1, NA)
  states = rowSums(2**(my.post$L-1-sset), na.rm=TRUE)
  binstates = unlist(lapply(states, DecToBin, len=my.post$L))
  freqs = as.data.frame(table(binstates))
  df = data.frame(State = freqs$binstates, Probability = freqs$Freq/sum(freqs$Freq))
  return(df)
}

#' Reconstruct ancestral states and create an annotated tree
#' 
#' `curate.tree` reads a phylogeny from file in newick tree format and a 
#' csv-file containing tip-labels and features. 
#' 
#' @param tree.filename a file containing a phylogeny in newick-tree format
#' @param data.filename a file with tip-labels in the first column and features
#'  in the remaining ones.
#' @param losses defaults to FALSE. Determines whether irreversible gains or 
#'  losses of features are considered.
#' @param data.header defaults to TRUE. Determines Whether data.filename is
#'  read with a header or not.
#' @return a list containing tree, data, transitions, srcs, dests and times
#' @export
curate.tree = function(tree.src, data.src, losses = FALSE, data.header=TRUE) {
  if(is.character(tree.src)) {
    # read in Newick tree and root
    my.tree = ape::read.tree(tree.src)
    my.rooted.tree = ape::root(my.tree, 1, resolve.root = TRUE)
  } else {
    my.rooted.tree = tree.src
  }
 
  if(is.character(data.src)) {
    # read in barcode data
    my.data = utils::read.csv(data.src, header=data.header)
    colnames(my.data)[1] = "label"
  } else {
    my.data = data.src
  }
 
  match.set = match(my.data$label, my.rooted.tree$tip.label)
  if(any(is.na(match.set))) {
    message("Found observations that didn't correspond to tips of the tree!")
    my.data = my.data[-which(is.na(match.set)),]
    match.set = match(my.data$label, my.rooted.tree$tip.label)
  }
  
  if(any(duplicated(my.data$label))) {
    message("Duplicates in observation set!")
  }
  
  # prune tree to include only those tips in the barcode dataset
  tree = ape::drop.tip(my.rooted.tree,
                  my.rooted.tree$tip.label[-match.set])
  
  tree$node.label = as.character(length(tree$tip.label) + 1:tree$Nnode)
  tree.labels = c(tree$tip.label, tree$node.label)
  
  if(any(duplicated(tree$tip.label))) {
    message("Duplicates in tree tips!")
  }
  
  cat("\n------- Painting ancestors...\n  ")
  
  # initialise "recursive" algorithm
  change = T
  new.row = my.data[1,]
  changes = data.frame()
  
  # while we're still making changes
  while(change == T) {
    change = F
    # loop through all nodes
    for(tree.ref in 1:length(tree.labels)) {
      this.label = tree.labels[tree.ref]
      # see if this node exists in our barcode dataset
      if(!(this.label %in% my.data$label)) {
        # if not, check to see if its children are all characterised
        descendant.refs = phangorn::Children(tree, tree.ref)
        if(all(tree.labels[descendant.refs] %in% my.data$label)) {
          
          ## ancestral state reconstruction
          # pull the rows in our barcode dataset corresponding to children of this node
          descendant.rows = which(my.data$label %in% tree.labels[descendant.refs])
          if(losses == FALSE) {
            # bitwise AND to construct the ancestral state
            new.barcode = apply(my.data[descendant.rows,2:ncol(my.data)], 2, prod)
          } else {
            # bitwise OR to construct the ancestral state
            new.barcode = apply(my.data[descendant.rows,2:ncol(my.data)], 2, max)
          }
          # add this ancestral state to the barcode dataset
          new.data = new.row
          new.data$label[1] = this.label
          new.data[1,2:ncol(new.data)] = new.barcode
          my.data = rbind(my.data, new.data)
          
          ## adding transitions to our observation set
          # loop through children
          for(d.ref in descendant.refs) {
            # pull the barcodes and branch lengths together and add to data frame
            d.row = which(my.data$label == tree.labels[d.ref])
            # get the reference for this edge
            e.ref = which(tree$edge[,2] == d.ref)
            changes = rbind(changes, 
                            data.frame(from=paste0(new.data[2:ncol(new.data)], collapse=""),
                                       to=paste0(my.data[d.row,2:ncol(new.data)], collapse=""),
                                       time=tree$edge.length[e.ref],
                                       from.node=this.label,
                                       to.node=tree.labels[d.ref]))
          }
          # we made a change, so keep the loop going
          change = T
        }
      }
    }
  }
  srcs = matrix(as.numeric(unlist(lapply(changes$from, strsplit, split=""))), byrow=TRUE, ncol=ncol(new.data)-1)
  dests = matrix(as.numeric(unlist(lapply(changes$to, strsplit, split=""))), byrow=TRUE, ncol=ncol(new.data)-1)
  
  rL = list("tree" = tree,
            "data" = my.data,
            "transitions" = changes,
            "srcs" = srcs,
            "dests" = dests,
            "times" = changes$time)
  return(rL)
}

#' Visualize phylogenetic tree with features
#' 
#' Creates a plot
#' 
#' @param curated.tree combined phylogeny and feature data created with [curate.tree()]
#' @return a ggplot
#' @export
plotHypercube.curated.tree = function(tree.set, scale.fn = geom_treescale(y=20, linesize=3, width =0.01)) {
  data.m = tree.set$data[,2:ncol(tree.set$data)]
  rownames(data.m) = tree.set$data[,1]
  data.m = tree.set$data[1:length(tree.set$tree$tip.label), 2:ncol(tree.set$data)]
  rownames(data.m) = tree.set$data$label[1:length(tree.set$tree$tip.label)]
  this.plot = ggtree::gheatmap(ggtree::ggtree(tree.set$tree) + scale.fn, data.m, low="white", high="#AAAAAA",
                               colnames_angle=90, hjust=0) +
    ggplot2::theme(legend.position="none")
  return(this.plot)
}
