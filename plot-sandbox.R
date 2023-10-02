#### x range, scale and log type to be edited

library(ggplot2)
library(gridExtra)
source("plot-trans.R")
library(stringdist)
library(pheatmap)
library(igraph)
library(ggraph)
library(ggplotify)

#### bubble and hypercube plots for XOR experiments

fname = c("test-xor-mod--1", "test-xor-mod-2", "test-xor-mod-3"); bigL=4
fname = c("test-cross-mod-1", "test-cross-mod-2", "test-cross-mod-3"); bigL=5
fname = 
bdf = thdf = data.frame()
g.bubbles = g.cube = g.gcube = g.motif = g.pheatmap = list()
for(i in 1:length(fname)) {
  bubble.name = paste(c("VerifyData/", fname[i], "-posterior.txt-bubbles.csv"), collapse="")
  tmpdf = read.csv(bubble.name)
  tmpdf$Expt=i
  bdf = rbind(bdf, tmpdf)
  thist.name = paste(c("VerifyData/", fname[i], "-posterior.txt-timehists.csv"), collapse="")
  tmpdf = read.csv(thist.name)
  tmpdf$Expt=i
  thdf = rbind(thdf, tmpdf)
g.bubbles[[i]] = ggplot(bdf, aes(x=Time, y=OriginalIndex, size=Probability)) +
  geom_point() +
  theme_light()

trans.name = paste(c("VerifyData/", fname[i], "-trans.txt"), collapse="")
states.name = paste(c("VerifyData/", fname[i], "-states.txt"), collapse="")
trans.1 = read.csv(trans.name, sep=" ")
trans.s.1 = read.csv(states.name, sep=" ")
g.cube[[i]] = plot.hypercube3(trans.1, statesdf=trans.s.1, 
                            node.labels = FALSE, seg.labels = TRUE, threshold=5e-2)

trans.1$Flux = trans.1$Probability*trans.s.1$Probability[trans.1$From+1]

trans.p = trans.1[trans.1$Flux > 1e-4,]
trans.g = graph_from_data_frame(trans.p)
bs = unlist(lapply(as.numeric(V(trans.g)$name), DecToBin, len=bigL))
V(trans.g)$binname = bs
layers = str_count(bs, "1")
g.gcube[[i]] = ggraph(trans.g, layout="sugiyama", layers=layers) + geom_edge_link(aes(edge_width=Flux, edge_alpha=Flux)) + 
  geom_node_point() + geom_node_label(aes(label=name)) + theme_graph() #aes(label=bs)) + theme_graph() 


routes.name = paste(c("VerifyData/", fname[i], "-posterior.txt-routes.txt"), collapse="")
routes = read.table(routes.name)

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

sf = 2
png("plot-sandbox.png", width=1600*sf, height=800*sf, res=72*sf)
grid.arrange(g.bubbles[[1]], g.gcube[[1]], g.motif[[1]], g.pheatmap[[1]], 
             g.bubbles[[2]], g.gcube[[2]], g.motif[[2]], g.pheatmap[[2]],
             g.bubbles[[3]], g.gcube[[3]], g.motif[[3]], g.pheatmap[[3]], nrow=3)
dev.off()

