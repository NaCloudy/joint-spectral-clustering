#######################################
# Main simulation results from "Joint Spectral Clustering in
# Multilayer Degree Corrected Blockmodles"
#######################################

# Load all methods for simulations
source("Experiments/run_all_methods.R")

# Note: the  code excludes the method graph-tool by default.
# To run graph-tool, install the Python package and uncomment
# the corresponding lines in "R/run_all_methods.R"

#######################################
# Simulation settings
#######################################

num_replications <- 2 #100
num_layers <- list(1,2) #list(1, 5, 10, 15, 20, 25,30,35, 40, 45, 50)

parameters_list <- num_layers
param_iter = parameters_list


#######################################
# Run different scenarios
#######################################

######## added simulations 8 and 9 ###########
# Same B all the same theta
results_simulation8 <- iterate_parameters(simulation8, parameters_list, param_iter, num_replications)
# Diff B all the same theta
results_simulation9 <- iterate_parameters(simulation9, parameters_list, param_iter, num_replications)
###############################################

# Same B same theta
results_simulation1 <- iterate_parameters(sim_setting = simulation1, parameters_list, param_iter, num_replications)
# Different B same theta
results_simulation2 <- iterate_parameters(sim_setting = simulation2, parameters_list, param_iter, num_replications)

# Diff B diff theta
results_simulation3 <- iterate_parameters(simulation3, parameters_list, param_iter, num_replications)
# Same B different theta
results_simulation4 <- iterate_parameters(simulation4, parameters_list, param_iter, num_replications)

##############################################

parameters_list <- lapply(param_iter, 
                          function(x) c(x, 150, 3))
# Same B alternating theta
results_simulation6 <- iterate_parameters(simulation6, parameters_list, param_iter, num_replications)
# Different B alternating theta
results_simulation7 <- iterate_parameters(simulation7, parameters_list, param_iter, num_replications)


#######################################
# Combine results
#######################################
results_simulation8$scenarioB <- "Same B"   # added 
results_simulation9$scenarioB <- "Different B"  # added

results_simulation1$scenarioB <- "Same B"
results_simulation2$scenarioB <- "Different B"
results_simulation3$scenarioB <- "Different B"
results_simulation4$scenarioB <- "Same B"
results_simulation6$scenarioB <- "Same B"
results_simulation7$scenarioB <- "Different B"

results_simulation8$scenarioT <- "All Same \u0398"   # added
results_simulation9$scenarioT <- "All Same \u0398"  # added

results_simulation1$scenarioT <- "Same \u0398"
results_simulation2$scenarioT <- "Same \u0398"
results_simulation3$scenarioT <- "Different \u0398"
results_simulation4$scenarioT <- "Different \u0398"
results_simulation6$scenarioT <- "Alternating \u0398"
results_simulation7$scenarioT <- "Alternating \u0398"

different_scenarios <- rbind(
    results_simulation8, results_simulation9, # added
    results_simulation1, results_simulation2, results_simulation4, results_simulation3, 
    results_simulation6, results_simulation7)

different_scenarios$scenarioB <- factor(different_scenarios$scenarioB,
                                       levels = c("Same B", "Different B"))
different_scenarios$scenarioT <- factor(different_scenarios$scenarioT,
                                        levels = c("All Same \u0398", "Same \u0398", "Different \u0398", "Alternating \u0398"))

save(different_scenarios, file = "Results-allmethods-addA-rep2-miscerror.RData")


#######################################
# Plot simulation results
#######################################
load("./Results-allmethods-addA-rep2-miscerror.RData")
png("Simulation-rep2-8scenarios-flipped.png", width = 1200, height = 1500, res = 200)
make_ggplot_multipleBT2(different_scenarios, "Number of graphs", xbreaks = c(1, seq(10, 50, 10)),
                       methodnames = c("DC-MASE", "Sum of adj. matrices", "Bias-adjusted SoS",  "MASE", "OLMF", "graph-tool"))
dev.off()


