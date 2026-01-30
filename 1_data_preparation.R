rm(list = ls())
library(perla)
library(spatialDE) # it contains the data
library(SpatialExperiment) # to organize the data into a SpatialExperiment object
library(nnSVG)
library(scran)
library(scry)
library(ggplot2)
library(mclust)

set.seed(123)

# data collection ---------------------------------------------------------
count_matrix <- read.csv("Data/Rep11_MOB_0.csv", header = T, row.names = 1)
coordinates <- read.csv("Data/MOB_sample_info.csv", header = T, row.names = 1)

mismatching.spots <- c(54, 194)  # after manual search
count_matrix <- count_matrix[-mismatching.spots,]
all.equal(row.names(count_matrix), row.names(coordinates))
x <- SpatialExperiment(assays = SimpleList(counts = t(as.matrix(count_matrix))), 
                       spatialCoords = as.matrix(coordinates[,-3]))
rowData(x)$gene_names <- colnames(count_matrix)


# Create the adjacency matrix ---------------------------------------------

dists <- as.matrix(dist(coordinates[,-3]))
min_dists <- apply(dists, 1, min, na.rm = TRUE)  # we take 1.7 as it empirically works
W <- (dists < 1.7)*1
diag(W) <- 0
ggplot(coordinates, aes(x, y, col = as.factor(rowSums(W))))+geom_point(cex = 3)+theme_bw()
saveRDS(W, "proximity_matrix.RDS")


# filter genes ------------------------------------------------------------

# Remove genes expressed in fewer than 5% of spots
x <- x[rowSums(counts(x) > 0) > (ncol(x) * 0.05), ]
x <- logNormCounts(x)


# Rankings ----------------------------------------------------------------

# computing the binomial deviance
x <- devianceFeatureSelection(x, sorted = T)
rowData(x)$rank_binomial_deviance <- 1:nrow(x)
x <- x[1:500,]
x <- nullResiduals(x, assay = "counts",
                   fam = c("poisson"),
                   type = c("deviance"),
                   batch = NULL)

# computing the spatial expression raking
x <- nnSVG(input = x)

plot(rowData(x)$binomial_deviance)
plot(rowData(x)$rank_binomial_deviance, rowData(x)$rank, 
     main = round(cor(rowData(x)$rank_binomial_deviance, rowData(x)$rank),3))


# Dimensionality reduction ------------------------------------------------

x <- GLMPCA(object = x, L = 50)


# Export data -------------------------------------------------------------


saveRDS(x, file = "processed_data.RDS")

