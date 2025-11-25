#######################################
# Main simulation results from "Joint Spectral Clustering in
# Multilayer Degree Corrected Blockmodles"
#######################################

# Load all methods for simulations
source("./Experiments/run_all_methods.R")

# Note: the  code excludes the method graph-tool by default.
# To run graph-tool, install the Python package and uncomment
# the corresponding lines in "R/run_all_methods.R"

#######################################
# Simulation settings
#######################################

num_replications <- 100
delta <- list(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
parameters_list <- delta
param_iter = parameters_list

results_simulationC2 <- iterate_parameters(sim_setting = simulation_identifiability_delta, parameters_list, param_iter, num_replications)

save(results_simulationC2, file = "./Added-Simu-Eva/Results-addC2-rep100-miscerror.RData")


#######################################
# Plot simulation results
#######################################
load("./Added-Simu-Eva/Results-addC2-rep100-miscerror.RData")
png("./Added-Simu-Eva/Figures/Simulation-C2-rep100-flipped.png", width = 1200, height = 1500, res = 200)
make_ggplot_single(results_simulationC2, "Proportion of A-type collapse layers", xbreaks = c(1, seq(10, 50, 10)),
                       methodnames = c("DC-MASE", "Sum of adj. matrices", "Bias-adjusted SoS",  "MASE", "OLMF", "graph-tool"))
dev.off()


