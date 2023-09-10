library(ggplot2)

# compare bubble plot summaries with walker numbers and DT/CT
df.1 = read.csv("testbed-cross-1-posterior.txt-bubbles.csv") 
df.1$expt = 1; df.1$t = 1
df.2 = read.csv("testbed-cross-2-posterior.txt-bubbles.csv") 
df.2$expt = 2; df.2$t = 1
df.3 = read.csv("testbed-cross-3-posterior.txt-bubbles.csv") 
df.3$expt = 3; df.3$t = 1

df.4 = read.csv("testbed-cross-ct-1-posterior.txt-bubbles.csv") 
df.4$expt = 1; df.4$t = 2
df.5 = read.csv("testbed-cross-ct-2-posterior.txt-bubbles.csv") 
df.5$expt = 2; df.5$t = 2
df.6 = read.csv("testbed-cross-ct-3-posterior.txt-bubbles.csv") 
df.6$expt = 3; df.6$t = 2
df = rbind(df.1, df.2, df.3, df.4, df.5, df.6)

ggplot(df, aes(x=Time+expt/10,y=OriginalIndex+t/10,size=Probability,color=factor(expt+10*t))) + geom_point()

# compare likelihood values by walker number
lik.1 = read.csv("testbed-cross-1-lik.txt")
lik.1$expt = 1
lik.2 = read.csv("testbed-cross-2-lik.txt")
lik.2$expt = 2
lik.3 = read.csv("testbed-cross-3-lik.txt")
lik.3$expt = 3
lik.4 = read.csv("testbed-cross-sa-lik.txt")
lik.4$expt = 4
lik.5 = read.csv("testbed-cross-sgd-lik.txt")
lik.5$expt = 5

lik.df = rbind(lik.1, lik.2, lik.3)

ggplot(lik.df, aes(x=LogLikelihood1,y=LogLikelihood2, color=factor(expt))) + 
  geom_point() + 
  facet_wrap(~expt)

lik.ts.df = rbind(lik.1, lik.2, lik.3, lik.4, lik.5)

ggplot(lik.ts.df, aes(x=Step, y=LogLikelihood1, color=factor(expt))) +
  geom_line()

######

# compare bubble plot summaries with walker numbers and DT/CT
df.1 = read.csv("VerifyData/test-cross-1-posterior.txt-bubbles.csv") 
df.1$expt = 1; df.1$t = 1
df.2 = read.csv("VerifyData/test-cross-2-posterior.txt-bubbles.csv") 
df.2$expt = 2; df.2$t = 1
df.3 = read.csv("VerifyData/test-cross-3-posterior.txt-bubbles.csv") 
df.3$expt = 3; df.3$t = 1

df.4 = read.csv("VerifyData/test-cross-ct-1-posterior.txt-bubbles.csv") 
df.4$expt = 1; df.4$t = 2
df.5 = read.csv("VerifyData/test-cross-ct-2-posterior.txt-bubbles.csv") 
df.5$expt = 2; df.5$t = 2
df.6 = read.csv("VerifyData/test-cross-ct-3-posterior.txt-bubbles.csv") 
df.6$expt = 3; df.6$t = 2
df = rbind(df.1, df.2, df.3, df.4, df.5, df.6)

ggplot(df, aes(x=Time+expt/10,y=OriginalIndex+t/10,size=Probability,color=factor(expt+10*t))) + geom_point()

# compare likelihood values by walker number for direct time experiments
lik.1 = read.csv("VerifyData/test-cross-1-lik.txt")
lik.1$expt = 1
lik.2 = read.csv("VerifyData/test-cross-2-lik.txt")
lik.2$expt = 2
lik.3 = read.csv("VerifyData/test-cross-3-lik.txt")
lik.3$expt = 3
lik.4 = read.csv("VerifyData/test-cross-sa-lik.txt")
lik.4$expt = 4
lik.5 = read.csv("VerifyData/test-cross-sgd-lik.txt")
lik.5$expt = 5

lik.df = rbind(lik.1, lik.2, lik.3)

ggplot(lik.df, aes(x=LogLikelihood1,y=LogLikelihood2, color=factor(expt))) + 
  geom_point() + 
  facet_wrap(~expt)

lik.ts.df = rbind(lik.1, lik.2, lik.3, lik.4, lik.5)

ggplot(lik.ts.df, aes(x=Step, y=LogLikelihood1, color=factor(expt))) +
  geom_line()

# compare likelihood values by walker number for continuous time experiments
lik.1 = read.csv("VerifyData/test-cross-ct-1-lik.txt")
lik.1$expt = 1
lik.2 = read.csv("VerifyData/test-cross-ct-2-lik.txt")
lik.2$expt = 2
lik.3 = read.csv("VerifyData/test-cross-ct-3-lik.txt")
lik.3$expt = 3
lik.4 = read.csv("VerifyData/test-cross-ct-sa-lik.txt")
lik.4$expt = 4
lik.5 = read.csv("VerifyData/test-cross-ct-sgd-lik.txt")
lik.5$expt = 5

lik.df = rbind(lik.1, lik.2, lik.3)

ggplot(lik.df, aes(x=LogLikelihood1,y=LogLikelihood2, color=factor(expt))) + 
  geom_point() + 
  facet_wrap(~expt)

lik.ts.df = rbind(lik.1, lik.2, lik.3, lik.4, lik.5)

ggplot(lik.ts.df, aes(x=Step, y=LogLikelihood1, color=factor(expt))) +
  geom_line()

