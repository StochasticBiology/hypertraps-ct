library(ggplot2)
library(gridExtra)

# compare bubble plot summaries with walker numbers and DT/CT
df.1 = read.csv("../VerifyData/test-cross-1-bubbles.csv") 
df.1$expt = 1; df.1$t = 1; df.1$label = "DT n=200"
df.2 = read.csv("../VerifyData/test-cross-2-bubbles.csv") 
df.2$expt = 2; df.2$t = 1; df.2$label = "DT n=2000"
df.3 = read.csv("../VerifyData/test-cross-3-bubbles.csv") 
df.3$expt = 3; df.3$t = 1; df.3$label = "DT n=20"

df.4 = read.csv("../VerifyData/test-cross-ct-1-bubbles.csv") 
df.4$expt = 1; df.4$t = 2; df.4$label = "CT n=200"
df.5 = read.csv("../VerifyData/test-cross-ct-2-bubbles.csv") 
df.5$expt = 2; df.5$t = 2; df.5$label = "CT n=2000"
df.6 = read.csv("../VerifyData/test-cross-ct-3-bubbles.csv") 
df.6$expt = 3; df.6$t = 2; df.6$label = "CT n=20"
df = rbind(df.1, df.2, df.3, df.4, df.5, df.6)

g.bubbles = ggplot(df, aes(x=Time+expt/10,y=OriginalIndex+t/10,size=Probability,color=label)) + 
  geom_point() + xlab("Time") + ylab("Feature") + theme_light()

# compare likelihood values by walker number for direct time experiments
lik.1 = read.csv("../VerifyData/test-cross-1-lik.csv")
lik.1$expt = "MCMC n=200"
lik.2 = read.csv("../VerifyData/test-cross-2-lik.csv")
lik.2$expt = "MCMC n=2000"
lik.3 = read.csv("../VerifyData/test-cross-3-lik.csv")
lik.3$expt = "MCMC n=20"
lik.4 = read.csv("../VerifyData/test-cross-sa-lik.csv")
lik.4$expt = "SA"
lik.5 = read.csv("../VerifyData/test-cross-sgd-lik.csv")
lik.5$expt = "SGD"

lik.df = rbind(lik.1, lik.2, lik.3)

ggplot(lik.df, aes(x=LogLikelihood1,y=LogLikelihood2, color=factor(expt))) + 
  geom_point() + 
  facet_wrap(~expt)

lik.ts.df = rbind(lik.1, lik.2, lik.3, lik.4, lik.5)

g.dt.lik = ggplot(lik.ts.df, aes(x=Step, y=LogLikelihood1, color=expt)) +
  geom_line() + ylab("log L") + theme_light()

# compare likelihood values by walker number for continuous time experiments
lik.1 = read.csv("../VerifyData/test-cross-ct-1-lik.csv")
lik.1$expt = "MCMC n=200"
lik.2 = read.csv("../VerifyData/test-cross-ct-2-lik.csv")
lik.2$expt = "MCMC n=2000"
lik.3 = read.csv("../VerifyData/test-cross-ct-3-lik.csv")
lik.3$expt = "MCMC n=20"
lik.4 = read.csv("../VerifyData/test-cross-ct-sa-lik.csv")
lik.4$expt = "SA"
lik.5 = read.csv("../VerifyData/test-cross-ct-sgd-lik.csv")
lik.5$expt = "SGD"

lik.df = rbind(lik.1, lik.2, lik.3)

ggplot(lik.df, aes(x=LogLikelihood1,y=LogLikelihood2, color=factor(expt))) + 
  geom_point() + 
  facet_wrap(~expt)

lik.ts.df = rbind(lik.1, lik.2, lik.3, lik.4, lik.5)

g.ct.lik = ggplot(lik.ts.df, aes(x=Step, y=LogLikelihood1, color=expt)) +
  geom_line() + ylab("log L") + theme_light()

g.liks = grid.arrange(g.dt.lik, g.ct.lik, nrow=1)

sf = 2
png("plot-tests.png", width=800*sf, height=600*sf, res=72*sf)
grid.arrange(g.bubbles, g.liks, nrow=2)
dev.off()
