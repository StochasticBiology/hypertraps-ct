library(ggplot2)
library(ggpubr)

#### bubble and hypercube plots for test cases

fname = c("cross-0", "cross-1", "cross-2")
n.name = c("10", "40", "160")
bdf = data.frame()
for(i in 1:length(fname)) {
bubble.name = paste(c("../VerifyData/", fname[i], "-bubbles.csv"), collapse="")
tmpdf = read.csv(bubble.name)
tmpdf$Expt=i
tmpdf$Expt.name = n.name[i]
bdf = rbind(bdf, tmpdf)
}
g.bubbles = ggplot(bdf, aes(x=Time+Expt/10, y=OriginalIndex, size=Probability, color=factor(Expt.name, levels=n.name))) +
  geom_point() + labs(x="Ordering", y="Feature", size="Probability", color="Observations") +
  theme_light()

#### let's try to reproduce the previous paper figures

### easy cube with gains not losses with new code

df = read.csv("../VerifyData/easycube-posterior.txt", header=FALSE, sep=" ")

# old parameterisation
#colnames(df) = c("t0on0", "t1on1", "t2on2", "na0", "t0on1", "t0on2", "t1on0", "na1", "t1on2", "t2on0", "t2on1", "na2")
# new parameterisation
colnames(df) = c("t0on0", "t0on1", "t0on2", "t1on0", "t1on1", "t1on2", "t2on0", "t2on1", "t2on2")

df$t0.1 = exp(df$t2on2)
df$t0.2 = exp(df$t1on1)
df$t0.4 = exp(df$t0on0)
df$t1.3 = exp(df$t1on1 + df$t2on1)
df$t1.5 = exp(df$t0on0 + df$t2on0)
df$t2.3 = exp(df$t2on2 + df$t1on2)
df$t2.6 = exp(df$t0on0 + df$t1on0)
df$t4.5 = exp(df$t2on2 + df$t0on2)
df$t4.6 = exp(df$t1on1 + df$t0on1)
df$t3.7 = exp(df$t0on0 + df$t2on0 + df$t1on0)
df$t5.7 = exp(df$t1on1 + df$t2on1 + df$t0on1)
df$t6.7 = exp(df$t2on2 + df$t1on2 + df$t0on2)

nodes.df = data.frame(x = c(0, 1, 1, 2, 1, 2, 2, 3),
                      y = c(0, 1, 0, 1,-1, 0,-1, 0),
                      label = 0:7)
edge.src = c(0,0,0,1,1,2,2,4,4,3,5,6)
edge.dst = c(1,2,4,3,5,3,6,5,6,7,7,7)
edges.df = data.frame()
for(i in 1:length(edge.src)) {
  edges.df = rbind(edges.df, data.frame(x=nodes.df$x[edge.src[i]+1], y=nodes.df$y[edge.src[i]+1],
                                        xend=nodes.df$x[edge.dst[i]+1], yend=nodes.df$y[edge.dst[i]+1]))
}
means = c(mean(df$t0.1), mean(df$t0.2), mean(df$t0.4),
          mean(df$t1.3), mean(df$t1.5), mean(df$t2.3),
          mean(df$t2.6), mean(df$t4.5), mean(df$t4.6),
          mean(df$t3.7), mean(df$t5.7), mean(df$t6.7))
sds = c(sd(df$t0.1), sd(df$t0.2), sd(df$t0.4),
        sd(df$t1.3), sd(df$t1.5), sd(df$t2.3),
        sd(df$t2.5), sd(df$t4.5), sd(df$t4.6),
        sd(df$t3.7), sd(df$t5.7), sd(df$t6.7))
edges.df$label = ""
for(i in 1:nrow(edges.df)) {
  edges.df$label[i] = paste(c(signif(means[i], digits=2), "+-", signif(sds[i], digits=2)), collapse="")
}
true.df = read.csv("../VerifyData/synth-easycube.txt", header=FALSE, sep=" ")
edges.df$truelabel = ""
for(i in 1:nrow(edges.df)) {
  edges.df$truelabel[i] = true.df$V3[true.df$V1==edge.src[i] & true.df$V2==edge.dst[i]]
}

g.easy = ggplot() + geom_segment(data=edges.df, aes(x=x,y=y,xend=xend,yend=yend), color="#AAAAAA") +
  geom_point(data=nodes.df, aes(x=x,y=y), size=6) + 
  geom_text(data=nodes.df, aes(x=x,y=y,label=label), color="white") + 
  geom_text(data=edges.df, aes(x=(x+2*xend)/3, y=(y+2*yend)/3+0.05, label=label), color="red") +
  geom_text(data=edges.df, aes(x=(x+2*xend)/3, y=(y+2*yend)/3-0.05, label=truelabel), color="blue") +
  theme_void()


### hard cube with gains not losses with new code

df = read.csv("../VerifyData/hardcube-posterior.txt", header=FALSE, sep=" ")

# old parameterisation
#colnames(df) = c("t0on0", "t1on1", "t2on2", "na0", "t0on1", "t0on2", "t1on0", "na1", "t1on2", "t2on0", "t2on1", "na2")
# new parameterisation
colnames(df) = c("t0on0", "t0on1", "t0on2", "t1on0", "t1on1", "t1on2", "t2on0", "t2on1", "t2on2")

df$t0.1 = exp(df$t2on2)
df$t0.2 = exp(df$t1on1)
df$t0.4 = exp(df$t0on0)
df$t1.3 = exp(df$t1on1 + df$t2on1)
df$t1.5 = exp(df$t0on0 + df$t2on0)
df$t2.3 = exp(df$t2on2 + df$t1on2)
df$t2.6 = exp(df$t0on0 + df$t1on0)
df$t4.5 = exp(df$t2on2 + df$t0on2)
df$t4.6 = exp(df$t1on1 + df$t0on1)
df$t3.7 = exp(df$t0on0 + df$t2on0 + df$t1on0)
df$t5.7 = exp(df$t1on1 + df$t2on1 + df$t0on1)
df$t6.7 = exp(df$t2on2 + df$t1on2 + df$t0on2)

nodes.df = data.frame(x = c(0, 1, 1, 2, 1, 2, 2, 3),
                      y = c(0, 1, 0, 1,-1, 0,-1, 0),
                      label = 0:7)
edge.src = c(0,0,0,1,1,2,2,4,4,3,5,6)
edge.dst = c(1,2,4,3,5,3,6,5,6,7,7,7)
edges.df = data.frame()
for(i in 1:length(edge.src)) {
  edges.df = rbind(edges.df, data.frame(x=nodes.df$x[edge.src[i]+1], y=nodes.df$y[edge.src[i]+1],
                                        xend=nodes.df$x[edge.dst[i]+1], yend=nodes.df$y[edge.dst[i]+1]))
}
means = c(mean(df$t0.1), mean(df$t0.2), mean(df$t0.4),
          mean(df$t1.3), mean(df$t1.5), mean(df$t2.3),
          mean(df$t2.6), mean(df$t4.5), mean(df$t4.6),
          mean(df$t3.7), mean(df$t5.7), mean(df$t6.7))
sds = c(sd(df$t0.1), sd(df$t0.2), sd(df$t0.4),
        sd(df$t1.3), sd(df$t1.5), sd(df$t2.3),
        sd(df$t2.5), sd(df$t4.5), sd(df$t4.6),
        sd(df$t3.7), sd(df$t5.7), sd(df$t6.7))
edges.df$label = ""
for(i in 1:nrow(edges.df)) {
  edges.df$label[i] = paste(c(signif(means[i], digits=2), "+-", signif(sds[i], digits=2)), collapse="")
}
true.df = read.csv("../VerifyData/synth-hardcube.txt", header=FALSE, sep=" ")
edges.df$truelabel = ""
for(i in 1:nrow(edges.df)) {
  edges.df$truelabel[i] = signif(true.df$V3[true.df$V1==edge.src[i] & true.df$V2==edge.dst[i]], digits=3)
}

g.hard = ggplot() + geom_segment(data=edges.df, aes(x=x,y=y,xend=xend,yend=yend), color="#AAAAAA") +
  geom_point(data=nodes.df, aes(x=x,y=y), size=6) + 
  geom_text(data=nodes.df, aes(x=x,y=y,label=label), color="white") + 
  geom_text(data=edges.df, aes(x=(x+2*xend)/3, y=(y+2*yend)/3+0.05, label=label), color="red") +
  geom_text(data=edges.df, aes(x=(x+2*xend)/3, y=(y+2*yend)/3-0.05, label=truelabel), color="blue") +
  theme_void()

hist.df = data.frame()
for(i in 1:12) {
  colref = 10+i
  trueval = as.numeric(edges.df$truelabel[i])
  hist.df = rbind(hist.df, data.frame(edge=i, probscale=df[,colref]/trueval))
}
g.hard.hist = ggplot(hist.df, aes(x=log10(probscale),fill=factor(edge))) + 
  geom_histogram(position="identity", alpha=0.2) + 
  labs(x="log10 Phat / P", y="Count", fill="Edge label") +
  theme_light()

#### analytic vs sampling simulation

rcdf = read.csv("../Verify/randomcubes.txt", header=FALSE, sep=" ")
g.timehist = ggplot(rcdf[rcdf$V1!=0,]) + 
  geom_line(aes(x=V2,y=V3, color=factor(V1))) +
  geom_point(aes(x=V2,y=V4, color=factor(V1)), size=5, alpha=0.2) +
  geom_point(aes(x=V2,y=V5, color=factor(V1))) +
  xlim(0,12) + labs(x="t", y="P(110,t|000,0)", color="Parameter set") +
  theme_light()

#### time histograms for inferred cross case
h1.df = read.csv("../VerifyData/cross-0-timehists.csv")
h1.df$expt=1
h2.df = read.csv("../VerifyData/cross-1-timehists.csv")
h2.df$expt=2
h3.df = read.csv("../VerifyData/cross-2-timehists.csv")
h3.df$expt=3
h.df = rbind(h1.df, h2.df, h3.df)

n.name = c("10", "40", "160")
h.df$expt.name = factor(n.name[h.df$expt], levels=n.name)
g.thist = ggplot(h.df, aes(x=Time, y=Probability, fill=factor(OriginalIndex+1))) + 
  geom_col(position="dodge") + xlim(-0.1,1.05) + facet_wrap(~expt.name, nrow=3) +
  labs(x="Time", y="Probability", fill="Feature") +
  theme_light() #+ scale_x_continuous(trans="log10") 

#### overall plot
sf = 2
png("plot-verify.png", width=800*sf, height=800*sf, res=72*sf)
ggarrange(g.timehist, g.easy, g.hard, g.hard.hist, g.bubbles, g.thist, 
          labels = c("A", "B", "C", "D", "E", "F"), nrow=3, ncol=2)
dev.off()
