setwd("..")
source("hypertraps.R")
setwd("Scripts")

#### various test bed experiments

for(expt in 1:5) {
  # cross
  if(expt == 1) {
    fname = c("test-cross-mod-1", "test-cross-mod-2", "test-cross-mod-3"); 
    oneshot="no"; bigL=5
  }
  # hi-order logic (normal inference)
  if(expt == 2) {
    fname = c("test-ho-mod--1", "test-ho-mod-2", "test-ho-mod-3"); 
    oneshot="no"; bigL=5
  }
  # hi-order logic (longer MCMC chain)
  if(expt == 3) {
    fname = c("test-ho-mod--1-l", "test-ho-mod-2-l", "test-ho-mod-3-l"); 
    oneshot="no"; bigL=5
  }
  # hi-order logic (simulated annealing)
  if(expt == 4) {
    fname = c("test-ho-mod--1-sa", "test-ho-mod-2-sa", "test-ho-mod-3-sa"); 
    oneshot="yes"; bigL=5
  }
  # cross (simulated annealing, regularised)
  if(expt == 5) {
    fname = c("test-cross-mod--1-sa", "test-cross-mod-2-sa", "test-cross-mod-3-sa"); 
    oneshot="regularised"; bigL=5
  }
  
  g.lik = g.bubbles = g.cube = g.reg = list()  
  # loop over files in output
  for(i in 1:length(fname)) {
    this.fname = paste(c("../VerifyData/", fname[i]), collapse="")
    if(oneshot=="regularised") {
      this.post = readHyperinf(this.fname, postlabel = this.fname, regularised = TRUE)
      # get trajectories from regularisation, if appropriate
      g.reg[[i]] = plotHypercube.regularisation(this.post)
    } else {
      this.post = readHyperinf(this.fname, postlabel = this.fname)
    }
    g.lik[[i]] = plotHypercube.lik.trace(this.post)
    g.bubbles[[i]] = plotHypercube.bubbles(this.post)
    g.cube[[i]] = plotHypercube.sampledgraph(this.post)
  }    
  
  # produce summary output
  out.name = paste(c("plot-sandbox-", expt, ".png"), collapse="")
  sf = 2
  png(out.name, width=1100*sf, height=800*sf, res=72*sf)
  if(oneshot == "regularised") {
    grid.arrange(g.reg[[1]], g.bubbles[[1]], g.cube[[1]],
                 g.reg[[2]], g.bubbles[[2]], g.cube[[2]],
                 g.reg[[3]], g.bubbles[[3]], g.cube[[3]],
                 nrow=3)
  } else {
    grid.arrange(g.lik[[1]], g.bubbles[[1]], g.cube[[1]],
                 g.lik[[2]], g.bubbles[[2]], g.cube[[2]],
                 g.lik[[3]], g.bubbles[[3]], g.cube[[3]],
                 nrow=3)
  }
  dev.off()
}



