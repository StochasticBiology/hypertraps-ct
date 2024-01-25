Two examples of data simulated under models expressed as DAGs with XOR, AND, OR relationships.

Files:
- XOR-AND-OR-X3-example.R
  - Model specification using DAGs and data generation (simulation) under the models.
  - Analysis with HyperTraPS-ct.

- true_models.pdf:
   - Graphical representation of the DAGs of the true models and the transitions between genotypes under the true models.

- mhn_hesbcn_analyses.R:
   - Analysis of the simulated data using MHN and H-ESBCN.
   
- mhn_hesbcn_plots.pdf   
  - Plots of the analysis using MHN and H-ESBCN.

- *.RData: RData files of simulated data or analysis results
  - sim_XOR_AND_OR-X3.RData: data generated (simulated) under the true models and used as input for HyperTraPS, MHN, H-ESBCN.
  - output_mhn_hesbcn.RData: output of the MHN and H-ESBCN, run from evamtools.
  - xao.runs.RData, x3.runs.RData: output of the analysis using HyperTraPS-ct.
