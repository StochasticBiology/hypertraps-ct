library(ggplot2)
library(ggpubr)

gset = list()
for(L in c(10, 30, 50, 70, 90)) {
  fname = paste(c("VerifyData/test-bigcross-hard-", L, "-posterior.txt-bubbles.csv"), collapse="")
  df = read.csv(fname)
  gset[[length(gset)+1]] = ggplot(df, aes(x=Time, y=OriginalIndex, 
                                          size=sqrt(Probability)/sqrt(max(Probability)), 
                                          alpha=sqrt(Probability)/sqrt(max(Probability)))) + 
    geom_point() + theme_light() + theme(legend.position="none") +
    scale_alpha_continuous(range=c(0,ifelse(L<50,1,0.3)))
}
sf = 2
png("plot-tests-L.png", width=800*sf, height=500*sf, res=72*sf)
ggarrange(plotlist=gset, nrow=2, ncol=3)
dev.off()

gset2 = list()
for(L in c(10, 30, 50, 70, 90)) {
  fname = paste(c("VerifyData/test-bigcross-hard-", L, "-posterior.txt-bubbles.csv"), collapse="")
  df = read.csv(fname)
  new.df = data.frame(sTime=rep(1:10, each=10), sIndex=rep(1:10, rep=10), sProbability=0)
  for(i in 1:nrow(df)) {
     ref = which(new.df$sTime == round(df$Time[i]/L*10) & 
                   new.df$sIndex == round(df$OriginalIndex[i]/L*10))
     new.df$sProbability[ref] = new.df$sProbability[ref]+df$Probability[i]
  }  
  gset2[[length(gset2)+1]] = ggplot(new.df, aes(x=sTime, y=sIndex, 
                                          size=sProbability/max(sProbability), 
                                          alpha=sProbability/max(sProbability))) + 
    geom_point() + theme_light() + theme(legend.position="none") 
}
ggarrange(plotlist=gset2, nrow=2, ncol=3)
