# HyperTraPS(-CT)

Hypercubic transition path sampling: Flexible inference of accumulation pathways, in discrete or continuous time, under different model structures, using combinations of longitudinal, cross-sectional, and phylogenetically-linked observations.

General content
=========

Requirements
------
For the command-line-only version, you'll just need the ability to compile and run C code.

For the R version without helper and plotting functions, you'll need R and the `Rcpp` library.

For those helper and plotting functions, you'll need R with these libraries for full functionality: `Rcpp`, `ggplot2`, `ggpubr`, `ggraph`, `ggwordcloud`, `igraph`, `stringr`, `stringdist`.

For the specific scientific case studies involving mitochondrial and tuberculosis evolution, the data wrangling also needs Python with `ETE3`.

Setting up and running
-----------

This repo contains code to run HyperTraPS and its variants either from the command line or via R. First, pull the repo down and set the working directory to its root. For command-line work, compile the `hypertraps.c` source code, with (for example):

`gcc -o3 hypertraps.c -lm -o hypertraps.ce`

This will produce the executable `hypertraps.ce` (you can of course name it whatever you like), which can then be run from the command line with various arguments below.

For R work, you can load just the inference code using `Rcpp` alone with

`library(Rcpp)`  
`sourceCpp("hypertraps-r.cpp")`

or, to attach some useful helper and plotting functions (which depend on more libraries), you can use

`source("hypertraps.R")`

The core function in R is then `HyperTraPS` which can then be run from the R console with various arguments below.

Input
------

HyperTraPS deals fundamentally with *transitions between states labelled by binary strings*. This allows us to work with cross-sectional, longitudinally, and phylogenetically-coupled samples.

The fundamental data element that goes into HyperTraPS is an observed transition from a "before" state to an "after" state. In the case of cross-sectional data, the "before" state is assumed to be the state of all zeroes (0000...), corresponding to a state which has not acquired any features. For longitudinal and/or phylogenetic observations, "before" and "after" states must be specified.

HyperTraPS requires at least a matrix describing "after" states -- this is a required argument. In R, it is passed as a matrix; at the command line, it is passed as a file.

If a matrix supplying "before" states is absent, the data are assumed to be cross-sectional. For example, the matrix

`0 0 1`  
`0 1 1`

would reflect cross-sectional observations of states 001 and 011, implicitly corresponding to transitions 000->001 and 000->011.

A matrix of "before" states may be specified as a matrix in R or a file at the command line, using `initialstates`. For example, including the initial states

`0 0 1`  
`0 0 1`

with the above observations would now reflect the transitions 001->001 (i.e. remaining in state 001) and 001->011.

*If you're not interested in continuous time, uncertain data, priors, or old data formats, skip to the next section now.* 

For continuous-time inference, HyperTraPS works with a time window for each observed transition, specified via a start time and an end time. If the start time and end time are equal, the transition is specified as taking exactly that time. If start time = 0 and end time = Inf, the transition can take any amount of time, which is mathematically equivalent to the case without continuous time. In general, the start and end times specify an allowed set of durations for the given transition, allowing uncertain timings to be accounted for.

In R, these start and end times are vectors specified by `starttimes` and `endtimes`. At the command line, they are stored in files accessed by the `--starttimes` and `--endtimes` flags. In both cases, absent start times means that all start times are assumed to be zero; absent end times means that all end times are assumed to be Inf.

The digit 2 can be used to reflect uncertainty in a state. For example,

`0 2 1`

corresponds to an observation where the first feature is absent, the second feature may be present or absent, and the third feature is present.

HyperTraPS also accepts descriptions of prior distributions on parameters. For the moment these are assumed to be uniform in log parameter space, and are specified by the min and max for each distribution. At the command line these values should be passed as a single file, given by `--priors`, with two columns and N rows, where the ith row gives the minimum and maximum value for parameter i. In R this should be provided as a matrix `priors` with two columns and N rows. Remember that N, the number of parameters, will depend on the model structure chosen and the number of features in the system. For example, the 3-feature system above and the default model structure (2, referring to L^2 parameters) would have N = 9.

*If you're not interested in old data formats, skip to the next section now.* 

For compatibility with older studies, a different input format is possible at the command-line. When a single file containing a matrix of observations is provided, it is by default assumed that this gives cross-sectional observations. But the `--transitionformat` flag can be used to interpret this matrix as paired "before" and "after" observations. Now, the matrix should contain the "before" states as *odd* rows (starting from row 1) and the "after" states as *even* rows.

For example,

`0 0 1`  
`0 1 1`

with the `--transitionformat` flag would be interpreted as a transition 001->011 (odd row -> even row). Without the `--transitionformat` flag this would be interpreted as two independent observations, corresponding (as in the R case) to transitions 000->001 and 000->011.

Output
------

At the command line, HyperTraPS will output information about the progress of a run to the screen, and produce a set of files containing outputs from and descriptions of the inference process.

In R, HyperTraPS will output information about the progress of a run to the console, and return a named list containing outputs from and descriptions of the inference process.

*If you are just interested in plotting summary outputs from the inference process, skip to the next section now.*

The files from the command line version can be read into R (for plotting, analysis, etc) using `readHyperinf` from `hypertraps.R`. This reads a collection of output files and returns a named list. Similar, the named list structure in R can be written to files (for storage) using `writeHyperinf` from `hypertraps.R`. This takes a named list and produces the corresponding file set.

The output structures are

| Information | R | File output | Notes |
|-------------|---|-------------|-------|
| Number of features, model type, number of parameters, and likelihood traces over the run | *list*$lik.traces | *label*-lik.csv | |
| Best parameterisation found during run | *list*$best | *label*-best.txt | |
| (Posterior) samples of parameterisation | *list*$posterior.samples | *label*-posterior.txt | |
| Probabilities for individual states | *list*$dynamics$states | *label*-states.csv | Only produced for L < 16 |
| Probabilities for individual transitions | *list*$dynamics$trans | *label*-trans.csv | Only produced for L < 16 |
| Best parameterisation after regularisation | *list*$regularisation$best | *label*-regularised.txt | Optional |
| Number of parameters, likelihood, and information criteria during regularisation | *list*$regularisation$reg.process | *label*-regularising.csv | Optional |
| "Bubble" probabilities of trait *i* acquisition at ordinal time *j* | *list*$bubbles | *label*-bubbles.csv | |
| Histograms of times of trait *i* acquisition | *list*$timehists | *label*-timehists.csv | |
| Individual sampled routes of accumulation | *list*$routes | *label*-routes.txt | Matrix with L columns; jth element of a row gives the feature acquired at step j. Each row is a sampled trajectory; there are *samples_per_row* samples per output parameterisation in the (posterior) sample set |
| Transition times for individual sampled routes of accumulation | *list*$times | *label*-times.txt | Matrix with L columns, with elements corresponding to the timings of each of the steps in the routes matrix above |
| Dwelling statistics for individual sampled routes of accumulation | *list*$betas | *label*-betas.txt | Matrix with L columns, with elements corresponding to the "beta" (characteristic rate of departure) for each step in the routes matrix above |

The named list in R can be passed to plotting and prediction functions for analysis -- see "Visualising and using output" below.

Demonstration
-----------
A good place to start is `hypertraps-demos.R`, where the basic form of R commands for HyperTraPS, most of the more interesting arguments that can be provided, and several scientific case studies are demonstrated. The expected behaviour for most of this demonstration is in `docs/hypertraps-demos.html`, compiled in R Markdown format with `hypertraps-demos.Rmd`.

Arguments to HyperTraPS
------

HyperTraPS needs at least a set of observations. *This is the only essential input.* In R this should take the form of a matrix; on the command line it should be stored in a file. The table below gives other inputs that can be provided, with more technical points appearing in lower sections.

| Argument | R | Command-line | Default |
|----------|---|--------------|---------|
| Input data | obs=*matrix* | --obs *filename* | None (required) |
| Timings and initial states:||||
| Precursor states | initialstates=*matrix* | --initialstates *filename* | None; on command line can also be specified as odd-element rows in "Input data" |
| Time window start | starttimes=*vector* | --times *filename* | 0 |
| Time window end | endtimes=*vector* | --endtimes *filename* | starttimes if present (i.e. precisely specified times); otherwise Inf |
| More technical content: ||||
| Prior mins and maxs | priors=*matrix* | --priors *filename* | -10 to 10 in log space for each parameter (i.e. very broad range over orders of magnitude)
| Model structure | model=*N* | --model *N* | 2 |
| Number of walkers | walkers=*N* | --walkers *N* | 200 |
| Inference chain length | length=*N* | --length *N* | 3 |
| Perturbation kernel | kernel=*N* | --kernel *N*| 5 |
| Random seed | seed=*N* | --seed *N* | 1 |
| Gains (0) or losses (1) | losses=*N* | --losses *N* | 0 |
| Use APM (0/1) | apm=*N* | --apm | 0 |
| Use SA (0/1) | sa=*N* | --sa | 0 |
| Use SGD (0/1) | sgd=*N* | --sgd | 0 |
| Use PLI (0/1) | pli=*N* | --pli | 0 |
| Gap between samples | samplegap=*N* | (not available yet) | 1e3 (>1e4 steps) or 1e2 (>1e2 steps)
| Regularisation by penalised likelihood | penalty=*X* | --penalty *X* | 0 (controls penalty per non-zero parameter) 
| Regularisation by LASSO | lasso=*X* | (not available yet) | 0 (controls lambda, LASSO penalty for absolute parameter values) 
| Number of simulations per parameter sample | samples_per_row=*N* | (not available yet) | 10 (number of samples to use for each parameter set when simulating routes and times for output)
| Stepwise regularise model after best parameterisation (0/1) | regularise=*N* | --regularise | 0 |
| Transition format observations | (not available) | --transitionformat | (off) |
| Output exact transitions (0/1) | output_transitions=*N* | --outputtransitions *N* | Switched off for L > 15 to avoid large output; consider switching off for CT 

So some example calls are (see the various demo scripts for more):

| Task | R | Command-line |
|------|---|--------------|
| Run HyperTraPS with default settings | HyperTraPS(*matrix*) | ./hypertraps.ce --obs *filename* |
| Run HyperTraPS-CT with default settings | HyperTraPS(*matrix*, starttimes=*vector*, endtimes=*vector*) | ./hypertraps.ce --obs *filename* --times *filename* --endtimes *filename* |
| Run HyperTraPS with all-edges model and regularise by penalised likelihood | HyperTraPS(*matrix*, model=-1, penalty=1) | ./hypertraps.ce --obs *filename* --model -1 --penalty 1 |
| Run HyperTraPS with all-edges model, then stepwise regularise | HyperTraPS(*matrix*, model=-1, regularise=1) | ./hypertraps.ce --obs *filename* --model -1 --regularise |

Visualising and using output
--------

The various outputs of HyperTraPS can be used in the R plotting functions below, which summarise (amongst other things) the numerical behaviour of the inference processes, the ordering and timings (where appropriate) of feature acquisitions, the structure of the learned hypercubic transition network, and any outputs from regularisation. *To start, `plotHypercube.summary` gives an overview of what we've learned. All of these except `plotHypercube.prediction` take a fitted model -- the output of `HyperTraPS` -- as a required argument, and may take other options as the table describes. `plotHypercube.prediction` takes a required argument that is the output of prediction functions described below.

| Plot function | Description | Options and defaults |
|---------------|-------------|---------|
| `plotHypercube.summary` | Summary plot combining several of the above | *f.thresh*=0.05 (flux threshold for graph plot), *t.thresh*=20 (time threshold for time histograms), *continuous.time*=TRUE (plot continuous time summary information) |
| More specific plots: |||
| `plotHypercube.lik.trace` | Trace of likelihood over inference run, re-calculated twice with different samples (to show consistency or lack thereof), along with current "in use" likelihood | |
| `plotHypercube.bubbles` | "Bubble plot" of probability of acquiring trait *i* at ordinal step *j* | *transpose*=FALSE (horizontal and vertical axis), *reorder*=FALSE (order traits by mean acquisition ordering) |
| `plotHypercube.motifs` | Motif-style plot of probability of acquiring trait *i* at ordinal step *j* |  |
| `plotHypercube.motifseries` | Motif-style plot of probability of specific states at a set of given snapshot times | *t.set*=0 (a set of snapshot times); *thresh*=0.05 (minimum probability for a state to be labelled) |
| `plotHypercube.graph` | Transition graph with edge weights showing probability flux (from full output) | *thresh*=0.05 (minimum threshold of flux for drawing an edge) |
| `plotHypercube.sampledgraph` | Transition graph with edge weights showing probability flux (from sampled paths) | *thresh*=0.05 (minimum threshold of flux for drawing an edge), max=1000 (maximum number of sampled routes to consider) |
| `plotHypercube.sampledgraph2` | As above, with mean and s.d. of absolute timings for each step | *thresh*=0.05 (minimum threshold of flux for drawing an edge), *max*=1000 (maximum number of sampled routes to consider), *no.times*=FALSE (avoid annotating edges with time information), *use.arc*=TRUE (arc edge format -- looks messier but less prone to overlapping edge labels), *node.labels*=TRUE (binary labels for nodes), *edge.label.size*=2 (font size for edge labels) |
| `plotHypercube.timehists` | Histograms of absolute timings for each trait's acquisition | |
| `plotHypercube.regularisation` | Information criterion vs number of nonzero parameters during regularisation | |
| `plotHypercube.timeseries` | Time series of acquisitions across sampled routes | |
| `plotHypercube.prediction` | Visualise predictions of unobserved features or future behaviour, given a model fit | *prediction* (required, the output of `predictHiddenVals` or `predictNextStep` (see below) |
| `plotHypercube.influences` | For the L^2 model, visualise how each feature acquisition influences the rate of acquisition of other features as a matrix |  *featurenames*=c("") (set of names for features); *use.regularised*=FALSE (use stepwise-regularised param set); *reorder*=FALSE (order features by base rate); *upper.right*=FALSE (control orientation of diagonal); *cv.thresh*=Inf (threshold posterior coefficient of variation, only plot interactions below this)|
| `plotHypercube.influencegraph` | For the L^2 or L^3 model, visualise how each feature acquisition influences the rate of acquisition of other features as a network |  as `plotHypercube.influences` |

Some useful ones are demonstrated here:
![image](https://github.com/StochasticBiology/hypertraps-ct/assets/50171196/153ed0d7-88ea-4dc2-a3bc-0c24b25923db)

In addition, we can ask HyperTraPS to make predictions about (a) any unobserved features for a given observation (for example, what value the ?s might take in 01??), and (b) what future evolutionary behaviour is likely given that we are currently in a certain state. These are achieved with functions `predictHiddenVals` and `predictNextStep` respectively. Both of these require a fitted hypercubic model (the output of `HyperTraPS`), which we'll call `fit`.

To predict unobserved values in a given observation, you can invoke

`predictHiddenVals(fit, state)`

where `fit` is the fitted model and `state` is a vector describing the incomplete observations with 0s, 1s, and 2s, the latter of which mark unobserved features. So 0122 has uncertain features at positions 3 and 4. The function will output two dataframes, one giving the probability of observing each specific state that could correspond to the incomplete observation, and one giving the aggregate probability of each uncertain feature being 1.

You can optionally provide an argument `level.weights` which provides probability weightings for different "levels" of the hypercube, that is, the different total number of 1s in the observation. This allows you to specify how likely it is that the true observation has acquired a certain number of features. By default this is uniform between the number of certain 1s and the maximum number of possible 1s.

To predict future evolutionary behaviour, you can invoke

`predictNextStep(fit, state)`

where `fit` is the fitted model and `state` is a given state. This will return a dataframe describing the possible future dynamics from `state` and their probabilities.

You can pass the output of both these prediction functions to `plotHypercube.prediction` for a visualisation of the corresponding prediction.

HyperTraPS will also estimate the probability of a particular state for either the discrete or the continuous time cases. In the discrete case, this is returned in `fit$dynamics`, where the probabilities of different states and different transitions are reported. The state probabilities are reported in normalised form (assuming uniform sampling over all the L+1 states in a complete trajectory) and unnormalised (the probability of a state given that we have acquired that particular number of features).

In the continuous case, the function `prob.by.time` will use the sampled dynamics from the inferred parameterisation(s) to output a set of state probabilities at a given sampling time. For example

`prob.by.time(fit, 0.1)` 

will return a set of states and their associated probabilities at time 0.1.

Specific content for introduction paper
=======

In the `Scripts/`, `Verify/`, `RawData/`, and `Process/` directories are various scripts, datafiles, and data-generating code to demonstrate HyperTraPS functionality and explore two scientific case studies (anti-microbial resistance in tuberculosis and cancer progression). There's also data from previous studies using HyperTraPS for comparison and demonstration. The demonstration files in the root directory make use of these datafiles.

Raw data
----
In `RawData/` there is externally-derived raw data:
  * Phylogeny and feature list for TB
  * Feature lists and name sets for previously published case studies: C4 photosynthesis evolution, ovarian cancer progression, severe malaria clinical progression, tool use evolution

Producing and curating data
----

The code to generate/curate these datasets is in `Scripts/`. Run at least these from the command line:
  * `infer-verify.sh` -- generates verification datasets using code in `Verify/` (and runs HyperTraPS on the command-line)
  * `infer-tests.sh` -- generates tests of various parameters and experiments using code in `Verify/` (and runs HyperTraPS on the command-line)
  * `prepare-all.sh` -- Bash script using the code below in `Process/` to set up TB dataset. TB is processed using `cook-data.sh`.

You can also run HyperTraPS from the command-line for the TB and other scientific case studies (but you can also do this within R)
  * `infer-tb.sh` -- runs HyperTraPS on the command-line for TB case study (this can also be done in R, in `tb-case-study.R` below)
  * `infer-others.sh` -- runs HyperTraPS on the command-line for other scientific case studies (this can also be done in R in the demo script)

in `Verify` there is some code to generate synthetic datasets
  * `generate-easycube.c` -- C code to generate an easy L=3 case
  * `generate-hardcube.c` -- C code to generate a more difficult L=3 case
  * `generate-cross.c` -- C code to generate several L=5 cases supporting competing pathways

in `Process/` there is some code to curate the phylogenetic TB data:
  * `prune-tree.py` -- Python code to reduce a raw Newick tree to just those nodes contained a barcode datafile
  * `internal-labels.c` -- C code to introduce dummy internal node labels (for use in followup preparation)
  * `parse-new.py` -- Python code to infer internal node barcodes and produce transition datafiles ready for HyperTraPS (and summary graphic for checking)
  * `cook-data.sh` -- Bash script applying these steps to given data

Analysing and plotting data
----

After setting up the data with the scripts above
  * `tb-case-study.R` -- runs HyperTraPS and plots TB case study. 
  * `cancer-examples.R` -- R script running examples of cancer progression analysis  
  * `plot-verify.R` -- plot verification study (from `infer-verify.sh`)
  * `plot-tests.R` -- plot test studies (from `infer-tests.sh`)
  * `plot-trees.R` -- plot the TB phylogeny

