rm(list = ls())


# dependencies ------------------------------------------------------------

source("Functions/dependencies.R")


# Data preparation --------------------------------------------------------

source("1_data_preparation.R")

x <- readRDS("processed_data.RDS")
W <- readRDS("proximity_matrix.RDS")


# Gene selection: Spatially-variable --------------------------------------

# using log-counts
source("2_SVgenes_logcountss.R")
x <- readRDS("Results/SVgenes_logcounts.RDS")
res <- readRDS("Results/perlaRun_SVgenes_logcounts.RDS")

