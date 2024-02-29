#################

# if we just want the HyperTraPS code without other dependencies, we only need Rcpp
# library(Rcpp)
# sourceCpp("hypertraps-r.cpp")

# if we want various helper functions and ggplot functions in addition, use this
source("hypertraps.R")

m.1 = matrix(rep(c(0,0,0,0,0,
               1,0,0,0,0,
               1,1,0,0,0,
               1,1,1,0,0,
               1,1,1,1,0,
               0,0,0,0,0,
               0,0,0,0,1,
               0,0,0,1,1,
               0,0,1,1,1,
               0,1,1,1,1),5), byrow=TRUE, ncol=5)
m.2 = matrix(rep(c(1,0,0,0,0,
               1,1,0,0,0,
               1,1,1,0,0,
               1,1,1,1,0,
               1,1,1,1,1,
               0,0,0,0,1,
               0,0,0,1,1,
               0,0,1,1,1,
               0,1,1,1,1,
               1,1,1,1,1),5), byrow=TRUE, ncol=5)
times.cs = rep(c(0.1, 0.2, 0.3, 0.4, 0.5), 10)
times = rep(0.1, 50)

### simple demo
my.post = HyperTraPS(m.2, initialstates = m.1, 
                     starttimes = times, endtimes = times, 
                     length = 4,
                     samplegap = 10,
                     output_transitions = 1,
                     featurenames = c("A", "B", "C", "D", "E")); 
plotHypercube.summary(my.post)
plotHypercube.sampledgraph2(my.post, thresh=0.1, use.arc=FALSE, edge.label.size=3) + 
  theme(legend.position="none") + expand_limits(x = c(-0.1, 1.1))
plotHypercube.influences(my.post, cv.thresh = Inf)
plotHypercube.influencegraph(my.post, cv.thresh = 1)
plotHypercube.motifseries(my.post, c(0.001, 0.01, 0.5, 1, 2, 5))

sf = 2
png("plot-new-fig-3.png", width=600*sf, height=500*sf, res=72*sf)
ggarrange(plotHypercube.sampledgraph2(my.post, thresh=0.1, use.arc=FALSE, edge.label.size=3, 
                                      edge.label.angle = "none", node.labels=TRUE,
                                      no.times=TRUE, small.times=TRUE,
                                      times.offset = c(0.1,0.2)) + 
            theme(legend.position="none") + expand_limits(x = c(-0.05,1.1)) ,
          plotHypercube.influences(my.post, cv.thresh = Inf),
          plotHypercube.influencegraph(my.post, cv.thresh = 1),
          plotHypercube.motifseries(my.post, c(0.001, 0.01, 0.5, 1, 2, 5)),
          labels = c("A", "B", "C", "D")
          )
dev.off()

# demonstrate predictions of future behaviour
prediction.step = predictNextStep(my.post, c(1,1,0,0,0))
plotHypercube.prediction(prediction.step)

# demonstrate predictions of hidden values...
# ... with no assumptions about progress on the hypercube
prediction.hidden = predictHiddenVals(my.post, c(1,2,2,2,2))
plotHypercube.prediction(prediction.hidden)
# ... enforcing given belief about progress on the hypercube
prediction.hidden = predictHiddenVals(my.post, c(1,2,2,2,2), level.weight=c(0,0,1,0,0,0))
plotHypercube.prediction(prediction.hidden)

# impose priors -- here disallowing every pairwise effect
# prior format: n_param rows, 2 columns. [i,1] = min i; [i,2] = max i
# here we impose priors on the base rates: feature 1 > feature 5 >> all others
priors = matrix(0, ncol=2, nrow=5*5)
priors[,1] = -10
priors[,2] = 10
for(i in 0:4) {
  priors[i*5+i+1,1] = -10
  priors[i*5+i+1,2] = -10
}
priors[0*5+0+1,1] = 1
priors[0*5+0+1,2] = 1
priors[4*5+4+1,1] = 0
priors[4*5+4+1,2] = 0

my.post.priors = HyperTraPS(m.2, initialstates = m.1, 
                     starttimes = times, endtimes = times, 
                     priors = priors,
                     featurenames = c("A", "B", "C", "D", "E")); 
plotHypercube.summary(my.post.priors)
# compare to results from L model (independent features)  
my.post.model1 = HyperTraPS(m.2, initialstates = m.1, 
                            starttimes = times, endtimes = times, 
                            model = 1, kernel=3,
                            featurenames = c("A", "B", "C", "D", "E")); 
plotHypercube.summary(my.post.model1)


### different levels of uncertainty in timings
# precisely specified timings, as above
my.post.time.precise = HyperTraPS(m.2, initialstates = m.1, 
                               starttimes = times, endtimes = times, 
                               length = 3, outputinput = 1,
                               featurenames = c("A", "B", "C", "D", "E")); 
# infinite width time window for transitions (just inferring ordering)
my.post.time.inf = HyperTraPS(m.2, initialstates = m.1, 
                                starttimes = times*0, endtimes = times*Inf, 
                                length = 3, outputinput = 1,
                                featurenames = c("A", "B", "C", "D", "E"));
# finite time window for each uncertain transition time
my.post.time.uncertain = HyperTraPS(m.2, initialstates = m.1, 
                     starttimes = times*0.25, endtimes = times*4, 
                     length = 3, outputinput = 1,
                     featurenames = c("A", "B", "C", "D", "E")); 
ggarrange(plotHypercube.timehists(my.post.time.precise, t.thresh=3), 
          plotHypercube.timehists(my.post.time.uncertain, t.thresh=3),
          plotHypercube.timehists(my.post.time.inf, t.thresh=3),
          nrow=3)

# output selected demos to file
sf = 2
png("plot-demos-si.png", width=800*sf, height=800*sf, res=72*sf)
ggarrange(plotHypercube.prediction(prediction.step, max.size=15), plotHypercube.prediction(prediction.hidden, max.size=15),
          plotHypercube.timehists(my.post.time.uncertain, t.thresh=3),
          plotHypercube.timehists(my.post.time.inf, t.thresh=3),
          plotHypercube.sampledgraph2(my.post.priors, use.arc = FALSE) + theme(legend.position="none"), 
          ncol=2, nrow=3,
          labels=c("A","B","C","D","E"))
dev.off()

sf = 2
png("plot-demos-timings-si.png", width=800*sf, height=600*sf, res=72*sf)
ggarrange(plotHypercube.timehists(my.post.time.precise, t.thresh=3),
          plotHypercube.timehists(my.post.time.uncertain, t.thresh=3),
          plotHypercube.timehists(my.post.time.inf, t.thresh=3),
nrow=2, ncol=2, labels=c("A","B","C"))
dev.off()

fig.1a = plotHypercube.sampledgraph2(my.post, use.arc = FALSE, edge.label.size = 3) + 
  theme(legend.position="none") + 
  expand_limits(x = c(-0.5, 4))
fig.1b = plotHypercube.timehists(my.post)
fig.1c = plotHypercube.influences(my.post)
fig.1d = plotHypercube.timeseries(my.post)

png("plot-demos.png", width=800*sf, height=600*sf, res=72*sf)
ggarrange(fig.1a, fig.1b, fig.1c, fig.1d, labels=c("A","B","C","D"))
dev.off()

# write output to files
writeHyperinf(my.post, "simpledemo", my.post$L, postlabel = "simpledemo", fulloutput=TRUE)

# retrieve output from files
my.post.r = readHyperinf("simpledemo", postlabel = "simpledemo", fulloutput=TRUE)
plotHypercube.summary(my.post.r)

# run an example with fewer walkers
my.post.sparse = HyperTraPS(m.2, initialstates = m.1, 
                            starttimes = times, endtimes = times,
                            featurenames = c("A", "B", "C", "D", "E"), walkers = 2); 
plotHypercube.summary(my.post.sparse, t.thresh = 2)

# direct time run (no time window specified)
my.post.dt = HyperTraPS(m.2, initialstates = m.1, featurenames = c("A", "B", "C", "D", "E")); 
plotHypercube.summary(my.post.dt, continuous.time = FALSE)
ggarrange(plotHypercube.timehists(my.post.dt, t.thresh=3),
          plotHypercube.timehists(my.post.time.inf, t.thresh=3),
          nrow=2)

### various other demos
# other plots
plotHypercube.motifs(my.post)
plotHypercube.timeseries(my.post)
plotHypercube.graph(my.post)

# regularisation
my.post.regularise = HyperTraPS(m.2, initialstates = m.1, regularise = 1, walkers = 20)
plotHypercube.regularisation(my.post.regularise)
plotHypercube.summary(my.post.regularise, continuous.time = FALSE)

# simulated annealing output
my.post.sa = HyperTraPS(m.2, initialstates = m.1, sa = 1)
plotHypercube.summary(my.post.sa, continuous.time = FALSE)

# phenotypic landscape inference
my.post.pli = HyperTraPS(m.2, initialstates = m.1, pli = 1)
plotHypercube.summary(my.post.pli, continuous.time = FALSE)

# start with every edge parameterised, then regularise
my.post.bigmodel.regularise = HyperTraPS(m.2, initialstates = m.1, model = -1, regularise = 1, walkers = 20)
plotHypercube.regularisation(my.post.bigmodel.regularise)

# this example demonstrates different model choices -- the data is generated using a process where pairs of features influence other features
# the (inappropriate) L^1 and L^2 parameterisations cannot capture this, but the "all edge" (model -1) and L^3 parameterisations can
logic.mat = readLines("RawData/old-hi-order.txt")
logic.mat = do.call(rbind, lapply(strsplit(logic.mat, ""), as.numeric))
logic.mat = rbind(logic.mat, logic.mat)
logic.mat.i = readLines("RawData/old-hi-order-init.txt")
logic.mat.i = do.call(rbind, lapply(strsplit(logic.mat.i, ""), as.numeric))
logic.mat.i = rbind(logic.mat.i, logic.mat.i)
logic.starts = logic.mat.i
logic.ends = logic.mat

logic.post.m1 = HyperTraPS(logic.ends, initialstates = logic.starts, length = 4, model = -1, walkers = 20)
logic.post.m1r = HyperTraPS(logic.ends, initialstates = logic.starts, length = 4, model = -1, walkers = 20, regularise = 1)
logic.post.1 = HyperTraPS(logic.ends, initialstates = logic.starts, length = 4, model = 1, walkers = 20)
logic.post.2 = HyperTraPS(logic.ends, initialstates = logic.starts, length = 4, model = 2, walkers = 20)
logic.post.3 = HyperTraPS(logic.ends, initialstates = logic.starts, length = 4, model = 3, walkers = 20)

#plotHypercube.influencegraph(logic.post.3, cv.thresh = 0.4)

png("plot-demo-logic.png", width=800*sf, height=400*sf, res=72*sf)
ggarrange(plotHypercube.graph(logic.post.m1) + ggtitle("All edges") + theme(legend.position="none"),
        #  plotHypercube.graph(logic.post.1) + ggtitle("L") + theme(legend.position="none"),
          plotHypercube.graph(logic.post.2)+ ggtitle("L^2") + theme(legend.position="none"),
          plotHypercube.graph(logic.post.3)+ ggtitle("L^3") + theme(legend.position="none"),
          plotHypercube.graph(logic.post.m1r)+ ggtitle("All edges, regularised") + theme(legend.position="none"),
          plotHypercube.influencegraph(logic.post.3, cv.thresh = 0.4),
          plotHypercube.regularisation(logic.post.m1r))
dev.off()

g.logic = ggarrange(plotHypercube.graph(logic.post.m1) + theme(legend.position="none"),
      #    plotHypercube.graph(logic.post.1) +  theme(legend.position="none"),
          plotHypercube.graph(logic.post.2)+ theme(legend.position="none"),
          plotHypercube.graph(logic.post.3)+ theme(legend.position="none"),
          plotHypercube.graph(logic.post.m1r)+ theme(legend.position="none"),
          plotHypercube.influencegraph(logic.post.3, cv.thresh = 0.4) + theme(legend.position="none"),
          plotHypercube.regularisation(logic.post.m1r),
          labels=c("i, all edges", "ii, L^2", "iii, L^3", "iv, all + reg", "v", "vi"))

prediction.hidden = predictHiddenVals(my.post, c(1,2,2,2,2))
g.hidden = plotHypercube.prediction(prediction.hidden, max.size = 10)
g.step = plotHypercube.prediction(prediction.step, max.size=10)
g.priors = plotHypercube.sampledgraph2(my.post.priors, use.arc = FALSE) + theme(legend.position="none")
my.post.70 = readHyperinf("VerifyData/test-bigcross-hard-70", postlabel = "VerifyData/test-bigcross-hard-70")
g.big.70 = ggplot(my.post.70$bubbles, aes(x=Time, y=OriginalIndex, 
               size=Probability, alpha=Probability)) + 
  labs(x = "Ordering", y = "Feature") +
  geom_point() + theme_light() + theme(legend.position="none") +
  scale_alpha_continuous(range=c(0,1))

png("plot-demo-features.png", width=800*sf, height=800*sf, res=72*sf)
ggarrange(g.logic,
          ggarrange(
            g.big.70, 
            g.priors, 
            ggarrange(g.hidden, g.step, labels=c("D i", "D ii"), nrow=2),
            widths = c(1,0.5,1), labels=c("B", "C", ""), nrow=1), 
          labels=c("A", ""), nrow=2, heights=c(1,0.75))
dev.off()

#refs = which(logic.post.m1r$regularisation$best != -20)-1
#data.frame(state=floor(refs/5), locus=refs%%5)

# compare likelihoods across model structures
c(max(logic.post.m1$lik.traces$LogLikelihood1),
  max(logic.post.1$lik.traces$LogLikelihood1),
  max(logic.post.2$lik.traces$LogLikelihood1),
  max(logic.post.3$lik.traces$LogLikelihood1))

#### short-form examples from past studies -- should run in a few minutes and give approximations to the original results

# ovarian cancer case study reproduction
# traits are chromosomal aberrations, observations are independent patient samples
cgh.mat = readLines("RawData/ovarian.txt")
cgh.mat = do.call(rbind, lapply(strsplit(cgh.mat, ""), as.numeric))
cgh.names = as.vector(read.table("RawData/ovarian-names.txt", sep=","))[[1]]

my.post.cgh = HyperTraPS(cgh.mat, 
                        length = 4, outputinput = 1, 
                        featurenames = cgh.names) 
ggarrange(plotHypercube.lik.trace(my.post.cgh), 
          plotHypercube.bubbles(my.post.cgh, reorder=TRUE), 
          plotHypercube.sampledgraph2(my.post.cgh, no.times=TRUE, node.labels=FALSE, use.arc=FALSE), ncol=3)
plotHypercube.sampledgraph2(my.post.cgh, no.times=TRUE)

# C4 paper reproduction
# traits are physical/genetic features associated with C4, observations are (incomplete) phylogenetically independent intermediate species
c4.mat = as.matrix(read.table("RawData/c4-curated.csv", sep=","))
c4.names = as.vector(read.table("RawData/c4-trait-names.txt", sep=","))[[1]]

my.post.c4 = HyperTraPS(c4.mat, 
                        length = 4, outputinput= 1, 
                        featurenames = c4.names) 
ggarrange(plotHypercube.lik.trace(my.post.c4), plotHypercube.bubbles(my.post.c4, reorder=TRUE), nrow=2)

# malaria paper reproduction
# traits are clinical features, observations are (incomplete) independent patient presentations
malaria.df = read.csv("RawData/jallow_dataset_binary_with2s.csv")
malaria.mat = as.matrix(malaria.df[,2:ncol(malaria.df)])
malaria.names = as.vector(read.table("RawData/malaria-names.txt", sep=","))[[1]]

my.post.malaria = HyperTraPS(malaria.mat, 
                        length = 3, outputinput= 1, 
                        kernel = 2,
                        walkers = 20,
                        featurenames = malaria.names) 
ggarrange(plotHypercube.lik.trace(my.post.malaria), plotHypercube.bubbles(my.post.malaria, reorder=TRUE, transpose=TRUE), nrow=2)

# tool use paper reproduction
# traits are modes of tool use, observations are phylogenetically coupled species observations (phylogeny has been accounted for, giving transition pairs)
tools.mat = as.matrix(read.table("RawData/total-observations.txt-trans.txt"))
tools.names = as.vector(read.table("RawData/tools-names.txt"))[[1]]
tools.starts = tools.mat[seq(from=1, to=nrow(tools.mat), by=2),]
tools.ends = tools.mat[seq(from=2, to=nrow(tools.mat), by=2),]

my.post.tools = HyperTraPS(tools.ends, initialstates = tools.starts, 
                           length = 4, outputinput= 1, 
                           featurenames = tools.names) 
ggarrange(plotHypercube.lik.trace(my.post.tools), plotHypercube.bubbles(my.post.tools, reorder=TRUE, transpose=TRUE), nrow=2)
plotHypercube.sampledgraph(my.post.tools, max=100)

sf = 2
png("demo-science-plots.png", width=1200*sf, height=600*sf, res=72*sf)
ggarrange( plotHypercube.bubbles(my.post.cgh, reorder=TRUE) + xlab("Order") + ylab("Aberration") + ggtitle("Ovarian cancer progression"),
           plotHypercube.bubbles(my.post.c4, reorder=TRUE) + xlab("Order") + ylab("C4 feature") + ggtitle("C4 photosynthesis"),
           plotHypercube.bubbles(my.post.malaria, reorder=TRUE, transpose=TRUE) + xlab("Symptom") + ylab("Order") + ggtitle("Severe malaria progression"),
           plotHypercube.bubbles(my.post.tools, reorder=TRUE, transpose=TRUE) + xlab("Tool use mode") + ylab("Order") + ggtitle("Tool use emergence"))
dev.off()
