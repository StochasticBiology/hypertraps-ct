# The simulated data can be generated from evamtools,
# and this requires installing the package,
# or it can be generated from the web-based tool 
# https://www.iib.uam.es/evamtools/
# and then uploaded. Alternatively, and the easiest here,
# is to load the pregenerated data.

simulate_from_evamtools <- FALSE

if (simulate_from_evamtools) {
  require(evamtools)

  # Simplified from generate_sample_from_dag
  sim_data_evam <- function(edges,
                            noise = 0.01,
                            N = 500) {
    parent_set <- evamtools:::parent_set_from_edges(edges)
    gene_names <- names(parent_set)
    n_genes <- length(parent_set)
    dag_trm <- evamtools:::HESBCN_model_2_output(edges, parent_set)$HESBCN_trans_rate_mat
    dag_probs <- evamtools:::probs_from_trm(dag_trm)
    tmp_samples_as_vector <-
      evamtools:::genot_probs_2_pD_ordered_sample(x = dag_probs,
                                                  ngenes = n_genes,
                                                  gene_names = gene_names,
                                                  N = N,
                                                  out = "vector"
                                                  )
    data_with_noise <- evamtools:::genotypeCounts_to_data(tmp_samples_as_vector,
                                                          e = noise)
    csd_counts <- evamtools:::data_to_counts(data_with_noise, out="data.frame")
    return(list(csd_counts = csd_counts,
                data = data_with_noise))
  }

  # XOR, AND, OR
  # Same structure, but different lambdas, compared to default in the web GUI
  dag_XOR_AND_OR <- data.frame(From = c("Root", "Root", rep(c("A", "B"), 3)),
                               To = c("A", "B", rep(c("C", "D", "E"), c(2, 2, 2))),
                              Relation = rep(c("Single", "AND", "OR", "XOR"), rep(2, 4)),
                              Lambdas = c(3.0, 2.5, rep(c(1.5, 1.2, 1.3), c(2, 2, 2))))
  sim_XOR_AND_OR <- sim_data_evam(dag_XOR_AND_OR)

  # D depends with XOR on A, B, C
  dag_X3 <- data.frame(From = c(rep("Root", 3), "A", "B", "C"),
                       To = c("A", "B", "C", rep("D", 3)),
                       Relation = c(rep("Single", 3), rep("XOR", 3)),
                       Lambdas = c(1.6, 1.8, 2.2, rep(2.0, 3)))
  sim_X3 <- sim_data_evam(dag_X3)

  save(file = "sim_XOR_AND_OR-X3.RData", sim_X3, sim_XOR_AND_OR)
  
  # Simplified from plot_evam
  plot_model <- function(x, main = "") {
    method_info <- igraph::graph_from_data_frame(x[, c("From", "To")])
    parent_set <- evamtools:::parent_set_from_edges(x)
    evamtools:::plot_method(method_info = method_info,
                            parent_set = parent_set,
                            edges = x, 
                            method = main)
  }

  plot_hypercube_transitions <- function(x, main = "") {
    obj <- list(edges = x,
                parent_set = evamtools:::parent_set_from_edges(x))
    tmat <- as.matrix(evamtools:::cpm2tm(obj)$trans_mat_genots)
    evamtools:::plot_genot_fg(tmat, plot_type = "trans_mat")
  }
  
  pdf(file = "true_models.pdf", width = 11, height = 11)
  par(mfrow = c(2, 2))
  par(oma = c(2, 2, 2, 2))
  plot_model(dag_X3, "X3")
  plot_model(dag_XOR_AND_OR, "XOR_AND_OR")
  plot_hypercube_transitions(dag_X3)
  title("X3: \nTrue transitions between genotypes")
  plot_hypercube_transitions(dag_XOR_AND_OR)
  title("XOR_AND_OR: \nTrue transitions between genotypes")
  dev.off()
  
} else {
  # Load RData. Otherwise, provide the RDS files downloaded from evamtools
  load("sim_XOR_AND_OR-X3.RData")
  ## d_x3  <- readRDS(DAG_X3_data.rds)$data
  ## d_xao <- readRDS("DAG_A_O_X_data.rds")$data
  d_x3  <- sim_X3$data
  d_xao <- sim_XOR_AND_OR$data
  # Genotype frequencies
  sim_X3$csd_counts
  sim_XOR_AND_OR$csd_counts
}

pwd <- getwd()
setwd("../.")
source("hypertraps.R")
setwd(pwd)

## Analyses take about 50 and 80 minutes for x3 and xao
## datasets. The pre-run analyses are available as RData files
run_analyses <- FALSE
if (run_analyses) {
  library(parallel)
  ## (We could parallelise over data set too. Not done here
  ## to keep the examples separate)
  x3.runs <- mcmapply(HyperTraPS,
                      model = c(2, 3, -1),
                      MoreArgs = list(
                        obs = d_x3,
                        featurenames = c("A", "B", "C", "D"),
                        length = 5,
                        kernel = 3,
                        penalty = 1,
                        limited_output = 1),
                      SIMPLIFY = FALSE,
                      mc.cores = detectCores()
                      )
  save(file = "x3.runs.RData", x3.runs)

  xao.runs <- mcmapply(HyperTraPS,
                       model = c(2, 3, -1),
                       MoreArgs = list(
                         obs = d_xao,
                         featurenames = c("A", "B", "C", "D", "E"),
                         length = 5,
                         kernel = 3,
                         penalty = 1,
                         limited_output = 1),
                       SIMPLIFY = FALSE,
                       mc.cores = detectCores()
                       )
  save(file = "xao.runs.RData", xao.runs)

  ## Stepwise regularise. These take about 50 and 80 minutes
  x3_step.runs <- mcmapply(HyperTraPS,
                           model = c(2, 3, -1),
                          MoreArgs = list(
                            obs = d_x3,
                            featurenames = c("A", "B", "C", "D"),
                            length = 5,
                            kernel = 3,
                            penalty = 0,
                            regularise = 1,
                            limited_output = 1),
                          SIMPLIFY = FALSE,
                          mc.cores = detectCores())
  save(file = "x3_step.runs.RData", x3_step.runs)

  xao_step.runs <- mcmapply(HyperTraPS,
                            model = c(2, 3, -1),
                           MoreArgs = list(
                             obs = d_xao,
                             featurenames = c("A", "B", "C", "D", "E"),
                             length = 5,
                             kernel = 3,
                             penalty = 0,
                             regularise = 1,
                             limited_output = 1),
                           SIMPLIFY = FALSE,
                           mc.cores = detectCores())
  save(file = "xao_step.runs.RData", xao_step.runs)
} else {
  load("x3.runs.RData")
  load("xao.runs.RData")
  load("x3_step.runs.RData")
  load("xao_step.runs.RData")
}


## Check traces before continuing
## FIXME: I am not sure that what the traces are displaying is detailed in
## the README or anywhere else. Looking at the plots, they do not seem to be
## two separate, independent, chains (the two traces always trace very similar
## paths). And, as something for further work, would it make sense
## to have multiple chains, possibly parallelized, to better assess convergence
## and have more samples from the posterior?
do.call(ggarrange, lapply(x3.runs, plotHypercube.lik.trace))
do.call(ggarrange, lapply(x3_step.runs, plotHypercube.lik.trace))
do.call(ggarrange, lapply(xao.runs, plotHypercube.lik.trace))
do.call(ggarrange, lapply(xao_step.runs, plotHypercube.lik.trace))

## Easier names
x3.l2 <- x3.runs[[1]]
x3.l3 <- x3.runs[[2]]
x3.m1 <- x3.runs[[3]]

x3_step.l2 <- x3_step.runs[[1]]
x3_step.l3 <- x3_step.runs[[2]]
x3_step.m1 <- x3_step.runs[[3]]

xao.l2 <- xao.runs[[1]]
xao.l3 <- xao.runs[[2]]
xao.m1 <- xao.runs[[3]]

xao_step.l2 <- xao_step.runs[[1]]
xao_step.l3 <- xao_step.runs[[2]]
xao_step.m1 <- xao_step.runs[[3]]



########################################
##
##        x3 dataset
##
########################################


## Regularisation
do.call(ggarrange, lapply(x3_step.runs, plotHypercube.regularisation))



## As in cancer-examples.R, get more samples from the regularised parameterisation
## To make results comparable between stepwise regularisation and penalty =1,
## we will represent results using the output from PosteriorAnalysis.

## FIXME: It is unclear to me if all plots are using the regularised results:
## see this comment in cancer-examples.R: "add use.regularised to all plot functions"
## FIXME Calling PosteriorAnalysis with "use_regularised = 1" and then plotting
## those samples, is this the way to see the paths from the regularised results?
## FIXME "use_regularised" fails if the model has been fit using
## penalty. Makes sense, but should probably be documented?

x3.l2.more.samples <- PosteriorAnalysis(x3.l2,
                                        featurenames = c("A", "B", "C", "D"),
                                        samples_per_row = 1000)

x3.l3.more.samples <- PosteriorAnalysis(x3.l3,
                                        featurenames = c("A", "B", "C", "D"),
                                        samples_per_row = 1000)

x3.m1.more.samples <- PosteriorAnalysis(x3.m1,
                                        featurenames = c("A", "B", "C", "D"),
                                        samples_per_row = 1000)

## Same as above. Note "use_regularised = TRUE"
x3_step.l2.more.samples <- PosteriorAnalysis(x3_step.l2,
                                             featurenames = c("A", "B", "C", "D"),
                                             samples_per_row = 1000,
                                             use_regularised = TRUE)

x3_step.l3.more.samples <- PosteriorAnalysis(x3_step.l3,
                                             featurenames = c("A", "B", "C", "D"),
                                             samples_per_row = 1000,
                                             use_regularised = TRUE)

x3_step.m1.more.samples <- PosteriorAnalysis(x3_step.m1,
                                             featurenames = c("A", "B", "C", "D"),
                                             samples_per_row = 1000,
                                             use_regularised = TRUE)


## FIXME: I don't understand the number of parameters for L = -1
dim(x3.runs[[3]]$posterior.samples) ## 4^3 = 64.
## But the number of possible transitions is 32 = 2^(4-1) * 4
## Similarly for the example in the Rmd vignette.
## I am missing something obvious here.
lapply(x3.runs, function(x) dim(x$posterior.samples))
lapply(x3_step.runs, function(x) dim(x$posterior.samples))


####### Comparing stepwise regularisation with penalty = 1

## Influences. The results seem very different from those of MHN
## (compare with p. 1 of mhn_hesbcn_plots.pdf) and between
## types of regularisation
## FIXME: is this difference between regularisations expected?
ggarrange(plotHypercube.influences(x3.l2.more.samples,
                                   feature.names = c("A", "B", "C", "D"),
                                   upper.right = TRUE),
          plotHypercube.influences(x3_step.l2.more.samples,
                                   feature.names = c("A", "B", "C", "D"),
                                   upper.right = TRUE,
                                   use.regularised = TRUE)
          )

ggarrange(
  plotHypercube.sampledgraph2(x3.l2.more.samples, no.times = TRUE,
                              use.arc = FALSE, node.label.size = 4,
                              edge.label.size = 0),
  plotHypercube.sampledgraph2(x3_step.l2.more.samples, no.times = TRUE,
                              use.arc = FALSE, node.label.size = 4,
                              edge.label.size = 0)
)

ggarrange(
  plotHypercube.sampledgraph2(x3.l3.more.samples, no.times = TRUE,
                              use.arc = FALSE, node.label.size = 4,
                              edge.label.size = 0),
  plotHypercube.sampledgraph2(x3_step.l3.more.samples, no.times = TRUE,
                              use.arc = FALSE, node.label.size = 4,
                              edge.label.size = 0)
)

ggarrange(
  plotHypercube.sampledgraph2(x3.m1.more.samples, no.times = TRUE,
                              use.arc = FALSE, node.label.size = 4,
                              edge.label.size = 0),
  plotHypercube.sampledgraph2(x3_step.m1.more.samples, no.times = TRUE,
                              use.arc = FALSE, node.label.size = 4,
                              edge.label.size = 0)
) 



## The three models, penalty = 1.
## The full model captures the true dependency patterns for D (fourth locus)
## on second acquisition, but then it is incorrect on third acquisition
## (non-negligible transitions to 1011, 0111, 1101) and
## on the final one (but see below)

ggarrange(
  plotHypercube.sampledgraph2(x3.l2.more.samples, no.times = TRUE,
                              use.arc = FALSE, node.label.size = 4,
                              edge.label.size = 0) + ggtitle("L^2"),
  plotHypercube.sampledgraph2(x3.l3.more.samples, no.times = TRUE,
                              use.arc = FALSE, node.label.size = 4,
                              edge.label.size = 0) + ggtitle("L^3"),
  plotHypercube.sampledgraph2(x3.m1.more.samples, no.times = TRUE,
                              use.arc = FALSE, node.label.size = 4,
                              edge.label.size = 0) + ggtitle("Full")
)

## The three models, stepwise regularisation. This improves upon the penalty=1
## model in not showing transitions to 1011 (though 0111, 1101, and 1111 are
## shown).
ggarrange(
  plotHypercube.sampledgraph2(x3_step.l2.more.samples, no.times = TRUE,
                              use.arc = FALSE, node.label.size = 4,
                              edge.label.size = 0) + ggtitle("L^2"),
  plotHypercube.sampledgraph2(x3_step.l3.more.samples, no.times = TRUE,
                              use.arc = FALSE, node.label.size = 4,
                              edge.label.size = 0) + ggtitle("L^3"),
  plotHypercube.sampledgraph2(x3_step.m1.more.samples, no.times = TRUE,
                              use.arc = FALSE, node.label.size = 4,
                              edge.label.size = 0) + ggtitle("Full")
)



## All of them think D can be acquired last, but that is not possible.
## All of them, but mostly L^2 and L^3, also incorrectly think D can be acquired
## first. (These mistakes are not surprising, given previous plot)


## FIXME: Is it correct to call it on the output from
## PosteriorAnalysis? No major changes, though
ggarrange(
  plotHypercube.bubbles(x3.l2.more.samples) + ggtitle("L^2"),
  plotHypercube.bubbles(x3.l3.more.samples) + ggtitle("L^3"),
  plotHypercube.bubbles(x3.m1.more.samples) + ggtitle("Full")
)

## FIXME: As above: can we use PosteriorAnalysis output here?  The full model
## seems to be getting wrong the acquisition of the fourth character, in contrast
## to L^2 and L^3, where that acquisition takes place a very long time after the
## others.
ggarrange(
  plotHypercube.timeseries(x3.l2.more.samples,
                           featurenames = LETTERS[1:4]) + ggtitle("L^2"),
  plotHypercube.timeseries(x3.l3.more.samples,
                           featurenames = LETTERS[1:4]) + ggtitle("L^3"),
  plotHypercube.timeseries(x3.m1.more.samples,
                           featurenames = LETTERS[1:4]) + ggtitle("Full")
)


## L^3 is doing a better job than Full, even when it cannot really
## represent the true dependency. Acquiring D last, which is not possible
## under the true model, are happening a long time after the rest.
ggarrange(
  plotHypercube.bubbles(x3.l3.more.samples) + ggtitle("L^3"),
  plotHypercube.timeseries(x3.l3.more.samples,
                           featurenames = LETTERS[1:4]) + ggtitle("L^3"),
  plotHypercube.sampledgraph2(x3.l3.more.samples, no.times = FALSE,
                              featurenames = LETTERS[1:4],
                              use.arc = FALSE, node.label.size = 4,
                              edge.label.size = 4)
)
## But the full model is missing this.
ggarrange(
  plotHypercube.bubbles(x3.m1.more.samples) + ggtitle("Full"),
  plotHypercube.timeseries(x3.m1.more.samples,
                           featurenames = LETTERS[1:4]) + ggtitle("Full"),
  plotHypercube.sampledgraph2(x3.m1.more.samples, no.times = FALSE,
                              featurenames = LETTERS[1:4],
                              use.arc = FALSE, node.label.size = 4,
                              edge.label.size = 4)
)

## Repeat the above with stepwise regularisation
ggarrange(
  plotHypercube.bubbles(x3_step.l2.more.samples) + ggtitle("L^2"),
  plotHypercube.bubbles(x3_step.l3.more.samples) + ggtitle("L^3"),
  plotHypercube.bubbles(x3_step.m1.more.samples) + ggtitle("Full")
)

## Here, the full model is doing a much better job. The true forbidden
## acquisitions occur after very long times.
ggarrange(
  plotHypercube.timeseries(x3_step.l2.more.samples,
                           featurenames = LETTERS[1:4]) + ggtitle("L^2"),
  plotHypercube.timeseries(x3_step.l3.more.samples,
                           featurenames = LETTERS[1:4]) + ggtitle("L^3"),
  plotHypercube.timeseries(x3_step.m1.more.samples,
                           featurenames = LETTERS[1:4]) + ggtitle("Full")
)

## It can be confirmed here: all gains of D that we know are not possible
## have times of occurrence several orders of magnitude larger than the others.
ggarrange(
  plotHypercube.timeseries(x3_step.m1.more.samples,
                           featurenames = LETTERS[1:4]) + ggtitle("Full"),
  plotHypercube.sampledgraph2(x3_step.m1.more.samples, no.times = FALSE,
                              featurenames = LETTERS[1:4],
                              use.arc = FALSE, node.label.size = 4,
                              edge.label.size = 4)
)


## FIXME: Would it be possible to filter in sampledgraph* plots by time, so that
## we can filter out those transitions that take a long time?

## FIXME: I am now mentally filtering the sampledgraph one prompted by the
## timeseries one: the first and second events are gained quickly. Then two
## paths: gaining a third (any of the remaining, except D)

## FIXME: Can we get the expected times of gaining 1, 2, 3, ... events? This
## would help set the threshold for the filtering of long transitions.

## FIXME: I am implicitly thinking about the "transition to END", without
## formalizing it for now.



########################################
##
##        xao dataset
##
########################################

## Regularisation
do.call(ggarrange, lapply(xao_step.runs, plotHypercube.regularisation))

## Get samples from the posterior for all models.
xao.l2.more.samples <- PosteriorAnalysis(xao.l2,
                                         featurenames = c("A", "B", "C", "D", "E"),
                                         samples_per_row = 1000)

xao.l3.more.samples <- PosteriorAnalysis(xao.l3,
                                         featurenames = c("A", "B", "C", "D", "E"),
                                         samples_per_row = 1000)

xao.m1.more.samples <- PosteriorAnalysis(xao.m1,
                                         featurenames = c("A", "B", "C", "D", "E"),
                                         samples_per_row = 1000)

xao_step.l2.more.samples <- PosteriorAnalysis(xao_step.l2,
                                              featurenames = c("A", "B", "C", "D", "E"),
                                              samples_per_row = 1000,
                                              use_regularised = 1)

xao_step.l3.more.samples <- PosteriorAnalysis(xao_step.l3,
                                              featurenames = c("A", "B", "C", "D", "E"),
                                              samples_per_row = 1000,
                                              use_regularised = 1)

xao_step.m1.more.samples <- PosteriorAnalysis(xao_step.m1,
                                              featurenames = c("A", "B", "C", "D", "E"),
                                              samples_per_row = 1000,
                                              use_regularised = 1)

## Influences. These look more similar to MHN, and more similar
## between themselves, than in the x3 example.
ggarrange(plotHypercube.influences(xao.l2.more.samples,
                                   feature.names = c("A", "B", "C", "D", "E"),
                                   upper.right = TRUE),
          plotHypercube.influences(xao_step.l2.more.samples,
                                   feature.names = c("A", "B", "C", "D", "E"),
                                   upper.right = TRUE,
                                   use.regularised = TRUE)
          )

## In the rest, we focus on the L^3 and Full model, since we know the L^2
## cannot adequately recover the true model

## L^3
## Both types of regularisation do a fairly good job of showing mostly only genotypes that can exist
## for genotypes with 1, 2, and 3 mutations, and for 4 mutations most prob.
## flux goes through 11110, which is the only four mutations genotype that can
## exist.
ggarrange(
  plotHypercube.sampledgraph2(xao.l3.more.samples, no.times = TRUE,
                              use.arc = FALSE, node.label.size = 4,
                              edge.label.size = 0),
  plotHypercube.sampledgraph2(xao_step.l3.more.samples, no.times = TRUE,
                              use.arc = FALSE, node.label.size = 4,
                              edge.label.size = 0)
)
## Full: same as above
ggarrange(
  plotHypercube.sampledgraph2(xao.m1.more.samples, no.times = TRUE,
                              use.arc = FALSE, node.label.size = 4,
                              edge.label.size = 0),
  plotHypercube.sampledgraph2(xao_step.m1.more.samples, no.times = TRUE,
                              use.arc = FALSE, node.label.size = 4,
                              edge.label.size = 0)
) 

## Bubbles.
## All models are similar. 
ggarrange(
  plotHypercube.bubbles(xao.l2.more.samples) + ggtitle("L^2"),
  plotHypercube.bubbles(xao.l3.more.samples) + ggtitle("L^3"),
  plotHypercube.bubbles(xao.m1.more.samples) + ggtitle("Full")
)


## Generally OK but many forbidden transitions do not have huge times.
ggarrange(
  plotHypercube.bubbles(xao.l3.more.samples) + ggtitle("L^3"),
  plotHypercube.timeseries(xao.l3.more.samples,
                           featurenames = LETTERS[1:5]) + ggtitle("L^3"),
  plotHypercube.sampledgraph2(xao.l3.more.samples, no.times = FALSE,
                              featurenames = LETTERS[1:5],
                              use.arc = FALSE, node.label.size = 4,
                              edge.label.size = 4)
)

## Generally OK but many forbidden transitions do not have huge times.
ggarrange(
  plotHypercube.bubbles(xao.m1.more.samples) + ggtitle("Full"),
  plotHypercube.timeseries(xao.m1.more.samples,
                           featurenames = LETTERS[1:5]) + ggtitle("Full"),
  plotHypercube.sampledgraph2(xao.m1.more.samples, no.times = FALSE,
                              featurenames = LETTERS[1:5],
                              use.arc = FALSE, node.label.size = 4,
                              edge.label.size = 4)
)
