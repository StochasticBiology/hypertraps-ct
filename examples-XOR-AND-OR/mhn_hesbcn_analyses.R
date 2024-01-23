## Analysis of the XOR-AND-OR and X3 examples
## using MHN and H-ESBCN
library(evamtools)
load("sim_XOR_AND_OR-X3.RData")
d_xao <- sim_XOR_AND_OR$data
d_x3  <- sim_X3$data

## Analyses take a few minutes each.
## Results from HESBCN for d_xao are not always the same
## i.e., there is some run-to-run variability despite the large
## number of iterations
r_xao <- evam(d_xao, methods = c("MHN", "HESBCN"),
              hesbcn_opts = list(MCMC_iter = 5e6))

r_x3 <- evam(d_x3, methods = c("MHN", "HESBCN"),
             hesbcn_opts = list(MCMC_iter = 5e6))

pdf("mhn_hesbcn_plots.pdf", height = 12, width = 12)
plot_evam(r_xao, top_paths = 10, methods = c("MHN", "HESBCN"))
plot_evam(r_x3, top_paths = 10, methods = c("MHN", "HESBCN"))
dev.off()

save(file = "output_mhn_hesbcn.RData", r_xao, r_x3)
