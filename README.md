# Code for "Joint spectral clustering in multilayer degree corrected blockmodels" by Joshua Agterberg, Zachary Lubberts and Jesus Arroyo

This repository contains all the code required to reproduce the experiments and data analysis from the paper "Joint spectral clustering in multilayer degree corrected blockmodels". In particular, this code implements the degree-corrected multiple adjacency spectral embedding (DC-MASE) to perform community detection in multilayer networks.

## Overview

Given a collection of *L* adjacency matrices representing graphs with $n$ aligned vertices (e.g., multilayer networks), DC-MASE obtains a joint embedding matrix with n rows, where each row corresponds to a vertex in the networks. After this, communities are obtained by applying a clustering procedure (such as K-means) to partition the latent positions into K groups in order to estimate the memberships in the multilayer degree-corrected stochastic block model. 

The joint embedding is calculated by performing a separate *adjacency spectral embedding*  (ASE) for each graph, which consists in computing the eigendecomposition of each adjacency matrix, with a possible eigenvalue scaling, followed by a row-normalization step (such as dividing the rows by its L2 norm), and then performs a joint singular value decomposition of the concatenated row-normalized ASEs. A pictorial representation is presented below.

![dc-mase embedding](https://github.com/jesusdaniel/dcmase/blob/main/img/DC-MASE2.png?raw=true)

# Data
The networks encode weighted edges representing the monthly number of flights between pairs of US airports. These data were obtained from the [T-100 Segment (US Carriers Only)](https://www.transtats.bts.gov/Fields.asp?gnoyr_VQ=GEE) database US Bureau of Transportation Statistics (BTS). A post-processed version of this dataset is included in this repository. The DC-MASE algorithm is illustrated using these data.

<img src="https://github.com/jesusdaniel/dcmase/blob/main/Figures/USmap-K4-dcmase.png" height="400" />
<img src="https://github.com/jesusdaniel/dcmase/blob/main/Figures/Bmatrices-K4-dcmase.png" height="400" />
<img src="https://github.com/jesusdaniel/dcmase/blob/main/Figures/degreecorrections-K4-dcmase.png" height="400" />

## File overview section
The folder contains scripts to implement the methods, implement the experiment, run the experiments and generate figures, and additional files containing the results of the experiments. Descriptions of each file are listed below.

### Method
- *R/dcmase.R*: implements the degree-corrected adjacency spectral embedding.
- *R/comdet-dcmase.R*: implements a community detection method based on DC-MASE.
- *R/comdetmethods.R*: code implementation of alternative methods for multilayer community detection (including method by Paul and Chen (2020) contained in *R/Codes_Spectral_Matrix_Paul_Chen_AOS_2020.r*, and the method from the Python package graphtool, which is called via the script *R/run_graph_tool.R*)
- *R/SpectralMethods.R* and *R/getElbows.R*: contain auxiliary functions to perform spectral embeddings.
- *R/make_ggplot.R* and *R/make_dcsbm_plots.R*: implement functions to create plots from experiments.

### Experiments
- *Experiments/run_all_methods.R*: wrapper functions for generating simulated data and running all community detection methods on these simulations
- *Experiments/simulations sparsity.R*: wrapper for data generation in simulations as function of network sparsity (Supplementary Materials of the paper)
- *Experiments/extrasimulations.R*: wrapper for data generation in additional simulations (Supplementary Materials of the paper)

### Code to generate figures
- *Main - simulations.R*: Run simulation experiments from Section 4 of the paper and generates Figure 4.1. Approximate running time: less than 20 minutes running in parallel on 10 nodes.
- *Main - US airport data analysis.R*: Run data analysis from Section 5 of the paper and generates Figures 5.1 and 5.2. Approximate running time: less than 30 seconds.
- *Supplement - Simulations sparsity.R*: Run simulation experiments from Section I of the supplementary materials and generates Figure I.1.
- *Supplement - Single clustering comparison.R*: Run simulation experiments from Section H of the paper and generates Figure H.1. 
- *Supplement - Spherical scaled vs unscaled.R*: Run simulation experiments from Section H of the paper and generates Figure H.2. 
- *Supplement - US airport out of sample error.R*: Run data experiments from Section J of the paper and generates Figure J.1.

### Other files
- *Experiments/Results-allmethods-rep100-miscerror.RData*: Results to generate Figure 4.1 of the paper.
- *Experiments/Results-allmethods-rep100-sparsity.RData*: Results to generate Figure I.1 of the supplementary materials..\
- *Python/graphtool-script.py*: wrapper for running graphtool method.

## Requirements

To run this project successfully, ensure you have the following:

### R Version
- R (>= 4.2.1)

### Required R Packages
Several packages are required to run the code and generate the figures. These are listed below, including the version of each package that was used. To install these versions, run the following code.
The code requires the package listed below, and was run using the versions

```r
install.packages("remotes") # If not already installed
remotes::install_version("igraph", version = "1.4.2")
remotes::install_version("mclust", version = "5.4.7")
remotes::install_version("ScorePlus", version = "0.1")
remotes::install_version("reshape2", version = "1.4.4")
remotes::install_version("ggplot2", version = "3.3.5")
#remotes::install_version("grid", version = "4.1.1")
remotes::install_version("plyr", version = "1.8.6")
remotes::install_version("ggthemes", version = "4.2.4")
remotes::install_version("maps", version = "3.4.0")
remotes::install_version("mapdata", version = "2.3.0")
#remotes::install_version("parallel", version = "4.1.1")

install.packages("ScorePlus")
install.packages("clue")
install.packages("psych")
install.packages("rARPACK")
install.packages("reticulate")
install.packages("grid")
install.packages("parallel")
install.packages("igraph")
install.packages("mclust")
install.packages("ScorePlus")
install.packages("reshape2")
install.packages("ggplot2")
install.packages("plyr")
install.packages("ggthemes")
install.packages("maps")
install.packages("mapdata")
```

## Data Access

This project includes access to publicly available data.  The following datasets are openly accessible and included in this repository:

### **Available Dataset**
- **Name:** US airport network (2016-2021)
- **Format:** RData
- **Location:** `Data/US_airport_data.RData`
- **Description:** This dataset contains weighted adjacency matrices representing the monthly number of flights between US airports. The data was obtained from the T-100 Segment (US Carriers Only) database of the US Bureau of Transportation Statistics (BTS).
- **Contents:**
  - `Adj_list`: A list of 69 adjacency matrices, where each matrix represents flights for a given month.
  - `airport_coordinates`: A 343×2 matrix with latitude and longitude of each airport (for visualization).
  - `airport_names`: A vector of airport codes matching the order in `Adj_list`.

Example code to load the data in R and fit 4 communities.

```{r}
source("R/dcmase.R")
source("R/comdet-dcmase.R")
source("R/SpectralMethods.R")
source("R/make_dcsbm_plots.R")

load("Data/US_airport_data.RData")

# Symmetrize networks (total number of flights between locations)
Adj_list <- lapply(Adj_list, function(A) ((A+t(A))))

# Fit a 4-community multilayer DCSBM
comdcmase.res <- comdet_dcmase(Adj_list, K = 4, clustering_method = "kmeans", scaled = T)

# Plot map with airports colored by community
plot_usmap(comdcmase.res$community_memberships)
```

# References
Agterberg, J., Lubberts, Z., & Arroyo, J. (2025). Joint spectral clustering in multilayer degree-corrected stochastic blockmodels. Journal of the American Statistical Association, (just-accepted), 1-23. 
[![arXiv shield](https://img.shields.io/badge/arXiv-2212.05053-red.svg?style=flat)](https://arxiv.org/abs/2212.05053)
