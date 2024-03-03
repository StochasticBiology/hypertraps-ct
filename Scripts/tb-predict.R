setwd("..")
source("hypertraps.R")
setwd("Scripts")

require(parallel)

# read and curate data on transitions for TB case study
tb.set = curate.tree("../Data/ng.2878-S2.txt", "../Data/tuberculosis-v5-header-19-29.csv")
tb.srcs = tb.set$srcs
tb.dests = tb.set$dests
tb.times = tb.set$times*1000
tb.names = colnames(tb.set$data[2:ncol(tb.set$data)])

nsplits = 10
L = ncol(tb.set$srcs)

training.refs = list()
nobs = nrow(tb.set$srcs)
set.seed(1)
for(i in 1:nsplits) {
  training.refs[[i]] = sample(1:nobs, 0.75*nobs)
}

# for parallelisation -- with different control parameter "fork" this function will return different HyperTraPS experiments
parallel.fn = function(fork, srcs, dests, times, training.refs) {
  train.srcs = srcs[training.refs[[fork]],]
  train.dests = dests[training.refs[[fork]],]
  train.times = times[training.refs[[fork]]]
  return(HyperTraPS(train.dests, initialstates = train.srcs, endtimes = train.times, penalty = 1, seed = 1, length = 5, kernel = 3))
}

# run these experiments in parallel
parallelised.runs <- mcmapply(parallel.fn, fork=1:nsplits,
                               MoreArgs = list(src = tb.srcs,
                                               dests = tb.dests,
                                               times = tb.times,
                                               training.refs = training.refs),
                               SIMPLIFY = FALSE,
                               mc.cores = min(detectCores(), nsplits))

# labels for experiments and features
res.df = res.state.df = data.frame()
for(i in 1:nsplits) {
  test.refs = setdiff(1:nobs, training.refs[[i]])
  for(j in 1:length(test.refs)) {
    state = tb.dests[test.refs[j],]
    
    prev.state = tb.srcs[test.refs[j],]
    if(sum(state) - sum(prev.state) == 1) {
      prediction = predictNextStep(parallelised.runs[[i]], prev.state)
      true.dest = paste0(state, collapse="")
      predicted.ref = which(prediction$states == true.dest)
      predicted.ranking = nrow(prediction)+1-rank(prediction$probs)[predicted.ref]
      predicted.prob = prediction$probs[predicted.ref]
      predicted.mode = prediction$states[prediction$probs==max(prediction$probs)]
      tdf = data.frame(splitref = i,
                       true.dest = true.dest, 
                       predicted.mode = predicted.mode,
                       predicted.ranking=predicted.ranking, 
                       predicted.prob=predicted.prob)
      res.state.df = rbind(res.state.df, tdf)
    }
    tdf = as.data.frame(t(state))
    hidden.locus = sample(1:L, 1)
    true.val = state[hidden.locus]
    state[hidden.locus] = 2
    prediction = predictHiddenVals(parallelised.runs[[i]], state) 
    prob.1 = prediction$locus.probs$prob
    tdf$splitref = i
    tdf$sample = j
    tdf$hidden.locus = hidden.locus
    tdf$true.val = true.val
    tdf$prob.1 = prob.1
    tdf$hit = ifelse( (prob.1 > 0.5 & true.val == 1) | (prob.1 < 0.5 & true.val == 0), TRUE, FALSE)
    res.df = rbind(res.df, tdf)
    
    
  }
}
#ggplot(res.state.df, aes(x=predicted.ranking, fill=factor(splitref))) + 
g.future = ggplot(res.state.df[res.state.df$splitref==1,], aes(x=factor(predicted.ranking))) + #, fill=factor(splitref))) + 
  geom_histogram(stat="count") +
  labs(x="Prediction ranking of next step", y="Number of test observations", fill="Sample") +
  theme_minimal()

this.split = 1
test.refs = setdiff(1:nobs, training.refs[[this.split]])
x = as.data.frame(tb.dests[test.refs,])
colnames(x) = tb.names
x$sample = 1:nrow(x)
xm = melt(x, id.vars = c("sample"))
preds = data.frame(variable=tb.names[res.df$hidden.locus], splitref=res.df$splitref, sample=res.df$sample, hit=res.df$hit)
preds = preds[preds$splitref==this.split,]
g.hidden = ggplot() + geom_tile(data=xm, aes(x=variable, y=sample, fill=factor(value))) + 
  scale_fill_manual(values = c("white", "#CCCCCC")) +
  geom_point(data=preds, aes(x=variable, y=sample, color=hit)) + theme_minimal() +
  labs(x="Drug", y="Test observation", color="Prediction correct", fill="True values")

roc.df = data.frame()
threshes = c((0:200)/400, 0.4999, 0.49999)
for(thresh in threshes) {
  unclass = length(which( (res.df$prob.1 < 0.5+thresh) &
                          (res.df$prob.1 > 0.5-thresh) ))
  trues = length(which( (res.df$true.val==1 & res.df$prob.1 > 0.5+thresh) |
                            (res.df$true.val==0 & res.df$prob.1 < 0.5-thresh) ))
  falses = (nrow(res.df)-unclass) - trues
  roc.df = rbind(roc.df, data.frame(thresh=thresh,
                                    unclass=unclass/nrow(res.df),
                                    trues=trues/(trues+falses)))
}
g.roc = ggplot(roc.df, aes(x=unclass, y=trues)) + geom_point() +
  theme_light() +
  labs(x="Proportion of observations\nwithout predictions", y="Proportion of correct predictions")

sf = 2
png("plot-tb-predictions.png", width=800*sf, height=300*sf, res=72*sf)
ggarrange(g.future, g.hidden, g.roc, nrow=1, labels=c("A", "B", "C"), widths=c(1,2,1))
dev.off()
