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
    tmp_samples_as_vector <- evamtools:::genot_probs_2_pD_ordered_sample(x = dag_probs,
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
  ## d_xao <- readRDS("DAG_A_O_X_data.rds")$data
  ## d_x3  <- readRDS(DAG_X3_data.rds)$data
  d_xao <- sim_XOR_AND_OR$data
  d_x3  <- sim_X3$data
  # Genotype frequencies
  sim_XOR_AND_OR$csd_counts
  sim_X3$csd_counts
}

pwd <- getwd()
setwd("../.")
source("hypertraps.R")
setwd(pwd)

## Analyses take about 2 and 5 hours for x3 and xao
## datasets. The pre-run analyses are available as an RData
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
                        length = 4,
                        kernel = 3,
                        penalty = 1,
                        limited_output = 1),
                      SIMPLIFY = FALSE,
                      mc.cores = detectCores()
                      )
  save(file = "x3.runs.RData", x3.runs)

  
  ## L = 4 is taking > 70 h, and is probably not sensible
  ## The runs here take about 4.5 hours (2, 4.5 and 3.5)
  xao.runs <- mcmapply(HyperTraPS,
                       model = c(2, 3, -1),
                       MoreArgs = list(
                         obs = d_xao,
                         featurenames = c("A", "B", "C", "D", "E"),
                         length = 4,
                         kernel = 3,
                         penalty = 1,
                         limited_output = 1),
                       SIMPLIFY = FALSE,
                       mc.cores = detectCores()
                       )
  save(file = "xao.runs.RData", xao.runs)
} else {
  load("xao.runs.RData")
  load("x3.runs.RData")
}


########################################
##
##        x3 dataset
##
########################################

## Check traces
## The second one might benefit from a longer run?
do.call(ggarrange, lapply(x3.runs, plotHypercube.lik.trace))

## Regularisation
do.call(ggarrange, lapply(x3.runs, plotHypercube.regularisation))

## Influences. The results seem very different from those of MHN
## (compare with p.2 of mhn_hesbcn_plots.pdf)
plotHypercube.influences(x3.runs[[1]], feature.names = c("A", "B", "C", "D"),
                         upper.right = TRUE)

## FIXME: I don't understand the number of parameters for L = -1
dim(x3.runs[[3]]$posterior.samples) ## 4^3 = 64.
## But the number of possible transitions is 32 = 2^(4-1) * 4
## Similarly for the example in the Rmd vignette.
## I am missing something obvious here.

lapply(x3.runs, function(x) dim(x$posterior.samples))

## Easier names
x3.l2 <- x3.runs[[1]]
x3.l3 <- x3.runs[[2]]
x3.m1 <- x3.runs[[3]]


## As in cancer-examples.R, get more samples from the regularised parameterisation

## FIXME: It is unclear to me if all plots are using the regularised results:
## see this comment in cancer-examples.R: "add use.regularised to all plot functions"
## Calling PosteriorAnalysis with "use_regularised = 1" and then plotting
## those samples, is this the way to see the paths from the regularised results?
x3.l2.more.samples <- PosteriorAnalysis(x3.l2,
                                        featurenames = c("A", "B", "C", "D"),
                                        samples_per_row = 1000,
                                        use_regularised = 1)

x3.l3.more.samples <- PosteriorAnalysis(x3.l3,
                                        featurenames = c("A", "B", "C", "D"),
                                        samples_per_row = 1000,
                                        use_regularised = 1)

x3.m1.more.samples <- PosteriorAnalysis(x3.m1,
                                        featurenames = c("A", "B", "C", "D"),
                                        samples_per_row = 1000,
                                        use_regularised = 1)



## Regularisation does not seem to make much of a difference
ggarrange(
  plotHypercube.sampledgraph2(x3.l2,  no.times = TRUE, use.arc = FALSE,
                              node.label.size = 4, edge.label.size = 0),
  plotHypercube.sampledgraph3(x3.l2.more.samples, no.times = TRUE,
                              use.arc = FALSE, node.label.size = 4,
                              edge.label.size = 0)) 

ggarrange(
  plotHypercube.sampledgraph2(x3.l3,  no.times = TRUE, use.arc = FALSE,
                              node.label.size = 4, edge.label.size = 0),
  plotHypercube.sampledgraph3(x3.l3.more.samples, no.times = TRUE,
                              use.arc = FALSE, node.label.size = 4,
                              edge.label.size = 0))

ggarrange(
  plotHypercube.sampledgraph2(x3.m1,  no.times = TRUE, use.arc = FALSE,
                              node.label.size = 4, edge.label.size = 0),
  plotHypercube.sampledgraph3(x3.m1.more.samples, no.times = TRUE,
                              use.arc = FALSE, node.label.size = 4,
                              edge.label.size = 0))


## The three models, regularised.
## L2 is different from MHN
## The full model captures the true dependency patterns for D (fourth locus)
## on second acquisition, but then it is incorrect on third acquisition
## (non-negligible transitions to 1011, 0111, 1101) and
## on the final one (but see below)

ggarrange(
  plotHypercube.sampledgraph3(x3.l2.more.samples, no.times = TRUE,
                              use.arc = FALSE, node.label.size = 4,
                              edge.label.size = 0) + ggtitle("L^2"),
  plotHypercube.sampledgraph3(x3.l3.more.samples, no.times = TRUE,
                              use.arc = FALSE, node.label.size = 4,
                              edge.label.size = 0) + ggtitle("L^3"),
  plotHypercube.sampledgraph3(x3.m1.more.samples, no.times = TRUE,
                              use.arc = FALSE, node.label.size = 4,
                              edge.label.size = 0) + ggtitle("Full")
)


## Bubbles.
## All of them think D can be acquired last, but that is not possible.
## All of them also incorrectly think D can be acquired first.
## (These mistakes not surprising, given previous plot)

ggarrange(
  plotHypercube.bubbles(x3.l2) + ggtitle("L^2"),
  plotHypercube.bubbles(x3.l3) + ggtitle("L^3"),
  plotHypercube.bubbles(x3.m1) + ggtitle("Full")
)

## FIXME: Is it correct to call it on the output from
## PosteriorAnalysis? No major changes, though
ggarrange(
  plotHypercube.bubbles(x3.l2.more.samples) + ggtitle("L^2"),
  plotHypercube.bubbles(x3.l3.more.samples) + ggtitle("L^3"),
  plotHypercube.bubbles(x3.m1.more.samples) + ggtitle("Full")
)

## As above: can we use PosteriorAnalysis output here?
ggarrange(
  plotHypercube.timeseries(x3.l2.more.samples,
                         feature.names = LETTERS[1:4]) + ggtitle("L^2"),
  plotHypercube.timeseries(x3.l3.more.samples,
                           feature.names = LETTERS[1:4]) + ggtitle("L^3"),
  plotHypercube.timeseries(x3.m1.more.samples,
                           feature.names = LETTERS[1:4]) + ggtitle("Full")
)


## Trying to see how the three plots complement each other
## So D is gained last, but a very long time after all others
## (bottom plot: 60, compared to value less than 0.5)
## which turns this last transition suspicious?
ggarrange(
    plotHypercube.bubbles(x3.l2.more.samples) + ggtitle("L^2"),
    plotHypercube.timeseries(x3.l2.more.samples,
                             feature.names = LETTERS[1:4]) + ggtitle("L^2"),
    plotHypercube.sampledgraph3(x3.l2.more.samples, no.times = FALSE,
                                feature.names = LETTERS[1:4],
                                use.arc = FALSE, node.label.size = 4,
                                edge.label.size = 4)
)

## Same thing
ggarrange(
  plotHypercube.bubbles(x3.l3.more.samples) + ggtitle("L^3"),
  plotHypercube.timeseries(x3.l3.more.samples,
                           feature.names = LETTERS[1:4]) + ggtitle("L^3"),
  plotHypercube.sampledgraph3(x3.l3.more.samples, no.times = FALSE,
                              feature.names = LETTERS[1:4],
                              use.arc = FALSE, node.label.size = 4,
                              edge.label.size = 4)
)

## Even more clear: all gains of D that we know are not possible
## have a time 10^8.
ggarrange(
  plotHypercube.bubbles(x3.m1.more.samples) + ggtitle("Full"),
  plotHypercube.timeseries(x3.m1.more.samples,
                           feature.names = LETTERS[1:4]) + ggtitle("Full"),
  plotHypercube.sampledgraph3(x3.m1.more.samples, no.times = FALSE,
                              feature.names = LETTERS[1:4],
                              use.arc = FALSE, node.label.size = 4,
                              edge.label.size = 4)
)

## Easier to see times if single plot
plotHypercube.sampledgraph3(x3.m1.more.samples, no.times = FALSE,
                            feature.names = LETTERS[1:4],
                            use.arc = FALSE, node.label.size = 4,
                            edge.label.size = 4)


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

## Follows the same logic as for x3.

## Check traces
do.call(ggarrange, lapply(xao.runs, plotHypercube.lik.trace))

## Regularisation
do.call(ggarrange, lapply(xao.runs, plotHypercube.regularisation))

## Influences. The results seem very more similar to from those of MHN
## (compare with p. 1 of mhn_hesbcn_plots.pdf)
plotHypercube.influences(xao.runs[[1]],
                         feature.names = c("A", "B", "C", "D", "E"),
                         upper.right = TRUE)

## Easier names
xao.l2 <- xao.runs[[1]]
xao.l3 <- xao.runs[[2]]
xao.m1 <- xao.runs[[3]]



xao.l2.more.samples <- PosteriorAnalysis(xao.l2,
                                         featurenames = c("A", "B", "C", "D", "E"),
                                         samples_per_row = 1000,
                                        use_regularised = 1)

xao.l3.more.samples <- PosteriorAnalysis(xao.l3,
                                         featurenames = c("A", "B", "C", "D", "E"),
                                         samples_per_row = 1000,
                                        use_regularised = 1)

xao.m1.more.samples <- PosteriorAnalysis(xao.m1,
                                         featurenames = c("A", "B", "C", "D", "E"),
                                         samples_per_row = 1000,
                                        use_regularised = 1)

## Regularisation does not seem to make much of a difference
ggarrange(
  plotHypercube.sampledgraph2(xao.l2,  no.times = TRUE, use.arc = FALSE,
                              node.label.size = 4, edge.label.size = 0),
  plotHypercube.sampledgraph3(xao.l2.more.samples, no.times = TRUE,
                              use.arc = FALSE, node.label.size = 4,
                              edge.label.size = 0)) 

ggarrange(
  plotHypercube.sampledgraph2(xao.l3,  no.times = TRUE, use.arc = FALSE,
                              node.label.size = 4, edge.label.size = 0),
  plotHypercube.sampledgraph3(xao.l3.more.samples, no.times = TRUE,
                              use.arc = FALSE, node.label.size = 4,
                              edge.label.size = 0))

ggarrange(
  plotHypercube.sampledgraph2(xao.m1,  no.times = TRUE, use.arc = FALSE,
                              node.label.size = 4, edge.label.size = 0),
  plotHypercube.sampledgraph3(xao.m1.more.samples, no.times = TRUE,
                              use.arc = FALSE, node.label.size = 4,
                              edge.label.size = 0))

## The three models, regularised.
## All models capture the key patterns, but L^3 does a better job than Full

ggarrange(
  plotHypercube.sampledgraph3(xao.l2.more.samples, no.times = TRUE,
                              use.arc = FALSE, node.label.size = 4,
                              edge.label.size = 0) + ggtitle("L^2"),
  plotHypercube.sampledgraph3(xao.l3.more.samples, no.times = TRUE,
                              use.arc = FALSE, node.label.size = 4,
                              edge.label.size = 0) + ggtitle("L^3"),
  plotHypercube.sampledgraph3(xao.m1.more.samples, no.times = TRUE,
                              use.arc = FALSE, node.label.size = 4,
                              edge.label.size = 0) + ggtitle("Full")
)

## Bubbles.
## All models are similar. In contrast to x3, there are no ordinal
## times precluded for any event in the true model.
ggarrange(
  plotHypercube.bubbles(xao.l2.more.samples) + ggtitle("L^2"),
  plotHypercube.bubbles(xao.l3.more.samples) + ggtitle("L^3"),
  plotHypercube.bubbles(xao.m1.more.samples) + ggtitle("Full")
)


## Using three plots together for each model. Again, the time to the forbidden
## transition has a much larger time.
ggarrange(
  plotHypercube.bubbles(xao.l2.more.samples) + ggtitle("L^2"),
  plotHypercube.timeseries(xao.l2.more.samples,
                           feature.names = LETTERS[1:5]) + ggtitle("L^2"),
  plotHypercube.sampledgraph3(xao.l2.more.samples, no.times = FALSE,
                              feature.names = LETTERS[1:5],
                              use.arc = FALSE, node.label.size = 4,
                              edge.label.size = 4)
)

## Same thing
ggarrange(
  plotHypercube.bubbles(xao.l3.more.samples) + ggtitle("L^3"),
  plotHypercube.timeseries(xao.l3.more.samples,
                           feature.names = LETTERS[1:5]) + ggtitle("L^3"),
  plotHypercube.sampledgraph3(xao.l3.more.samples, no.times = FALSE,
                              feature.names = LETTERS[1:5],
                              use.arc = FALSE, node.label.size = 4,
                              edge.label.size = 4)
)

## Even more clear: all gains of D that we know are not possible
## have a time 10^8.
ggarrange(
  plotHypercube.bubbles(xao.m1.more.samples) + ggtitle("Full"),
  plotHypercube.timeseries(xao.m1.more.samples,
                           feature.names = LETTERS[1:5]) + ggtitle("Full"),
  plotHypercube.sampledgraph3(xao.m1.more.samples, no.times = FALSE,
                              feature.names = LETTERS[1:5],
                              use.arc = FALSE, node.label.size = 4,
                              edge.label.size = 4)
)





## FIXME: evamTools should give bubble plots for all models and for the true model 
