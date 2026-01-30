rm(list = ls())


# dependencies ------------------------------------------------------------

source("Functions/dependencies.R")


# Data preparation --------------------------------------------------------

source("1_data_preparation.R")

x <- readRDS("processed_data.RDS")
W <- readRDS("proximity_matrix.RDS")


# Gene selection: Spatially-variable --------------------------------------

# using log-counts
source("2_SVgenes_logcounts.R")

# using Poisson deviance residuals
source("2_SVgenes_PoissonResiduals.R")


# GLMPCA ------------------------------------------------------------------

source("2_GLMPCA.R")
x <- readRDS("Results/GLMPCA.RDS")
