library(phytools)
library(ggtree)
library(ggtreeExtra)

# this function converts a species name string from the Newick format which Common Taxonomy Tree gives us into a simpler lower-case, no quotes version comparable to Kostas' dataset
convname = function(str) {
  return(tolower(gsub("\'", "", gsub("_", " ", str))))
}

for(expt in c(1,2)) {
  
if(expt == 1) {
my.tree = read.tree("../Data/mro-ncbi-tree.phy")
my.data = read.csv("../Data/mro-barcodes.csv")
my.tree$tip.label = convname(my.tree$tip.label)
my.data$Organism = convname(my.data$Organism)
}
  if(expt == 2) {
my.tree = read.tree("../Data/ng.2878-S2.txt")
my.data = read.csv("../Data/tuberculosis-v5-header-19-29.csv")
my.data$Organism = my.data$Isolate
  }
  
p.df = data.frame()
for(i in 1:nrow(my.data)) {
  for(j in 2:ncol(my.data)) {
    if(my.data[i,j] == 1) {
      p.df = rbind(p.df, data.frame(label=my.data$Organism[i], feature=j))
    }
  }
}

#p.tree = ggtree(my.tree, layout="circular") + geom_tiplab2(offset=2,hjust=1)
if(expt == 1) {
p.tree = ggtree(my.tree, layout="circular") + 
  geom_fruit(geom=geom_text, offset=1, 
             mapping = aes(x=label, label=label), 
             size=3, color="#888888") 
  data.m = as.matrix(my.data[2:ncol(my.data)])
}
if(expt == 2){
  p.tree = ggtree(my.tree, branch.length="none",layout="circular")
  data.m = as.matrix(my.data[2:(ncol(my.data)-1)])
}

rownames(data.m) = my.data$Organism
fname = paste(c("tree-", expt, ".png"), collapse="")
sf = 2
png(fname, width=800*sf, height=800*sf, res=72*sf)
print(gheatmap(p.tree, data.m, low="white", high="black", colnames_angle=90) + theme(legend.position = "none"))
dev.off()

fname = paste(c("tree-", expt, "a.png"), collapse="")
png(fname, width=300*sf, height=500*sf, res=72*sf)
if(expt == 1) {
  p.tree = ggtree(my.tree) + geom_tiplab(aes(label=label), size=3)
  hoff = 13
  voff = 1
} else {
  p.tree = ggtree(my.tree) 
  hoff = 0
  voff = 50
}
print(gheatmap(p.tree, data.m, low="white", high="#AAAAAA", offset=hoff, colnames_angle=90, colnames_offset_y=voff, font.size=3, width=0.5) + theme(legend.position = "none"))
dev.off()
}
