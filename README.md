# HyperTraPS(-CT)

Hypercubic transition path sampling: Flexible inference of accumulation pathways, in discrete or continuous time, under different model structures, using combinations of longitudinal, cross-sectional, and phylogenetically-linked observations.

General content
=========

Requirements
------
For the command-line-only version, you'll just need the ability to compile and run C code.

For the R version without helper and plotting functions, you'll need R and the `Rcpp` library.

For those helper and plotting functions, you'll need R with these libraries: `Rcpp`, `ggplot2`, `ggpubr`, `ggraph`, `igraph`, `stringr`, `stringdist`.

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

Demonstration
-----------
A good place to start is `hypertraps-demos.R`, where the basic form of R commands for HyperTraPS, most of the more interesting arguments that can be provided, and several scientific case studies are demonstrated. 

Arguments to HyperTraPS
------

HyperTraPS needs at least a set of observations. In R this should take the form of a matrix; on the command line it should be stored in a file.

| Argument | R | Command-line | Default |
|----------|---|--------------|---------|
| Input data | matrix_arg=*matrix* | --obs *filename* | None (required) |
| Precursor states | initialstates_arg=*matrix* | (odd-element rows in "Input data") | None |
| Cross-sectional observations | (assumed if "Precursor states" absent) | --crosssectional | 0 |
| Time window start | starttimes_arg=*vector* | --times *filename* | 0 |
| Time window end | endtimes_arg=*vector* | --endtimes *filename* | Inf |
| Model structure | model_arg=*N* | --model *N* | 2 |
| Number of walkers | walkers_arg=*N* | --walkers *N* | 200 |
| Inference chain length | length_index_arg=*N* | --length *N* | 3 |
| Perturbation kernel | kernel_index_arg=*N* | --kernel *N*| 5 |
| Random seed | seed_arg=*N* | --seed *N* | 1 |
| Gains (0) or losses (1) | losses_arg=*N* | --losses *N* | 0 |
| Use APM (0/1) | apm_type_arg=*N* | --apm | 0 |
| Use SA (0/1) | sa_arg=*N* | --sa | 0 |
| Use SGD (0/1) | sgd_arg=*N* | --sgd | 0 |
| Use PLI (0/1) | PLI_arg=*N* | --PLI | 0 |
| Regularise model (0/1) | regularise_arg=*N* | --regularise | 0 |

So some example calls are (see the various demo scripts for more):

| Task | R | Command-line |
|------|---|--------------|
| Run HyperTraPS with default settings | HyperTraPS(*matrix*) | ./hypertraps.ce --obs *filename* |
| Run HyperTraPS-CT with default settings | HyperTraPS(*matrix*, starttimes_arg=*vector*, endtimes_arg=*vector*) | ./hypertraps.ce --obs *filename* --times *filename* --endtimes *filename* |
| Run HyperTraPS with all-edges model, then regularise | HyperTraPS(*matrix*, model_arg=-1, regularise_arg=1) | ./hypertraps.ce --obs *filename* --model -1 --regularise |

Plots in R
--------

| Plot function | Description | Options and defaults |
|---------------|-------------|---------|
| `plotHypercube.lik.trace` | Trace of likelihood over inference run, calculated twice (to show consistency or lack thereof) | |
| `plotHypercube.bubbles` | "Bubble plot" of probability of acquiring trait *i* at ordinal step *j* | transpose=FALSE (horizontal and vertical axis), reorder=FALSE (order traits by mean acquisition ordering) |
| `plotHypercube.motifs` | Motif-style plot of probability of acquiring trait *i* at ordinal step *j* | |
| `plotHypercube.graph` | Transition graph with edge weights showing probability flux (from full output) | thresh=0.05 (minimum threshold of flux for drawing an edge) |
| `plotHypercube.sampledgraph` | Transition graph with edge weights showing probability flux (from sampled paths) | thresh=0.05 (minimum threshold of flux for drawing an edge), max=1000 (maximum number of sampled routes to consider) |
| `plotHypercube.sampledgraph2` | As above, with mean and s.d. of absolute timings for each step | thresh=0.05 (minimum threshold of flux for drawing an edge), max=1000 (maximum number of sampled routes to consider), no.times=FALSE (avoid annotating edges with time information) |
| `plotHypercube.timehists` | Histograms of absolute timings for each trait's acquisition | |
| `plotHypercube.regularisation` | Information criterion vs number of nonzero parameters during regularisation | |
| `plotHypercube.timeseries` | Time series of acquisitions across sampled routes | |
| `plotHypercube.summary` | Summary plot combining several of the above | |

All but the last are demonstrated here:
![image](https://github.com/StochasticBiology/hypertraps-ct/assets/50171196/153ed0d7-88ea-4dc2-a3bc-0c24b25923db)

Specific content for introduction paper
=======

`./`
----
Scripts wrapping the curation, inference, and analysis process for different case. In each case, datasets are produced and set in their own directory. The analysis code is then called to process these data, and plotting code follows.
  * `infer-verify.sh` -- generates verification datasets and runs HyperTraPS
  * `infer-tests.sh` -- generates tests of various parameters and experiments and runs HyperTraPS
  * `prepare-all.sh` -- Bash script using the code below in `Process/` to set up TB and MRO datasets. TB, and MRO with NCBI phylogeny, are straightforwardly processed using `cook-data.sh`. MRO with TimeTree phylogeny is a bit more involved and has its own script `mro-timetree-parse.sh`.
  * `infer-tb.sh` -- runs HyperTraPS for TB data
  * `infer-mro.sh` -- runs HyperTraPS for (different versions of) MRO data
  * `infer-others.sh` -- runs HyperTraPS for various other scientific cases
    
`Verify/`
---------
Code to generate synthetic test datasets
  * `generate-easycube.c` -- C code to generate an easy L=3 case
  * `generate-hardcube.c` -- C code to generate a more difficult L=3 case
  * `generate-cross.c` -- C code to generate several L=5 cases supporting competing pathways
    
`RawData/`
----------
  * Phylogenies and feature lists for MRO and TB
  * Feature lists and name sets for previously published case studies: C4 photosynthesis evolution, ovarian cancer progression, severe malaria clinical progression, tool use evolution

`Process/`
----------
Code for (1), distilling transition data suitable for HyperTraPS, for raw MRO and TB data. 
  * `prune-tree.py` -- Python code to reduce a raw Newick tree to just those nodes contained a barcode datafile
  * `internal-labels.c` -- C code to introduce dummy internal node labels (for use in followup preparation)
  * `parse-new.py` -- Python code to infer internal node barcodes and produce transition datafiles ready for HyperTraPS (and summary graphic for checking)
  * `cook-data.sh` -- Bash script applying these steps to given data
  * `mro-timetree-parse.sh` -- Bash script taking care of some subtleties when combining MRO data with TimeTree phylogeny
