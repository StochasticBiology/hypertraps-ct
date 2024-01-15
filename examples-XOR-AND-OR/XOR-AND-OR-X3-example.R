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
                            N = 10000) {
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
  # save(file = "sim_XOR_AND_OR-X3.RData", sim_X3, sim_XOR_AND_OR)
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











