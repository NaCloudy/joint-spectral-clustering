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
num_layers <- list(1, 5, 10, 15, 20, 25,30,35, 40, 45, 50) #list(1, 2)

parameters_list <- num_layers
param_iter = parameters_list

results_simulationC1 <- iterate_parameters(sim_setting = simulation_identifiability_m, parameters_list, param_iter, num_replications)

save(results_simulationC1, file = "./Added-Simu-Eva/Results-addC1-rep100-miscerror.RData")


#######################################
# Plot simulation results
#######################################
load("./Added-Simu-Eva/Results-addC1-rep100-miscerror.RData")
png("./Added-Simu-Eva/Figures/Simulation-C1-rep100-flipped.png", width = 1200, height = 1500, res = 200)
make_ggplot_single(results_simulationC1, "Number of graphs", xbreaks = c(1, seq(10, 50, 10)),
                       methodnames = c("DC-MASE", "Sum of adj. matrices", "Bias-adjusted SoS",  "MASE", "OLMF", "graph-tool"))
dev.off()


