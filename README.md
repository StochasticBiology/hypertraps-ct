# hctdump
HyperTraPS(-CT)
===============

Fall 2023 work
-----
In progress:
* updating workhorse code to output sensible L2 parameterisations and best state/transition values -- DONE
* retesting with old data -- DONE to some extent. CT really confounds the inference -- simple systems harder to infer -- *not convinced this is working*

Next: 
* make plots in R to match old document -- in progress
* test with direct time comparison -- see above -- CT confounds inference
* SGD
* better command line interface


Inference of evolutionary and progressive pathways, in discrete or continuous time, using combinations of longitudinal, cross-sectional, and phylogenetically-linked observations.

**Requirements:** Bash, Python 3 with ETE3 http://etetoolkit.org/download/ , GCC, Gnuplot

The general inference process involves several steps:
  1. Distilling transition data suitable for HyperTraPS from the raw data of interest
  2. Running HyperTraPS to infer posterior transition matrices from these data
  3. Post hoc analysis of these posteriors
  
This repo contains the code for (2) and (3), code to generate synthetic data for verification, and raw data and code for (1) for two case study systems: mitochondrion-related organelles (MROs) and drug-resistant tuberculosis (TB).

`./`
----
Scripts wrapping the curation, inference, and analysis process for different case. In each case, datasets are produced and set in their own directory. The analysis code is then called from `Inference/` to process these data, and plotting code follows.
  * `infer-verify.sh` -- generates verification datasets and runs HyperTraPS
  * `analyse-verify.sh` -- post hoc analysis of outputs from above
  * `plot-verify.sh` -- calls Gnuplot to produce graphical output of this analysis
  * `prepare-all.sh` -- Bash script using the code below in `Process/` to set up TB and MRO datasets. TB, and MRO with NCBI phylogeny, are straightforwardly processed using `cook-data.sh`. MRO with TimeTree phylogeny is a bit more involved and has its own script `mro-timetree-parse.sh`.
  * `infer-tb.sh` -- runs HyperTraPS for TB data
  * `infer-mro.sh` -- runs HyperTraPS for (different versions of) MRO data
  
`Inference/`
------------
Code for the inference and post hoc analysis.
  * `hypertraps-all.c` -- C code for (2), takes prepared data and outputs posteriors on model parameters. Can be run in continuous or discrete time.
  * `posteriors.c` -- C code for (3), performs post hoc analysis summarising posterior distributions and behaviour
  
`Verify/`
---------
Code to generate synthetic test datasets
  * `generate-easycube.c` -- C code to generate an easy L=3 case
  * `generate-hardcube.c` -- C code to generate a more difficult L=3 case
  * `generate-cross.c` -- C code to generate several L=5 cases supporting competing pathways
    
`RawData/`
----------
  * Phylogenies and feature lists for MRO and TB

`Process/`
----------
Code for (1), distilling transition data suitable for HyperTraPS, for raw MRO and TB data. 
  * `prune-tree.py` -- Python code to reduce a raw Newick tree to just those nodes contained a barcode datafile
  * `internal-labels.c` -- C code to introduce dummy internal node labels (for use in followup preparation)
  * `parse-new.py` -- Python code to infer internal node barcodes and produce transition datafiles ready for HyperTraPS (and summary graphic for checking)
  * `cook-data.sh` -- Bash script applying these steps to given data
  * `mro-timetree-parse.sh` -- Bash script taking care of some subtleties when combining MRO data with TimeTree phylogeny
