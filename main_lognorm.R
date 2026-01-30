rm(list = ls())
library(perla)
library(spatialDE) # it contains the data
library(SpatialExperiment) # to organize the data into a SpatialExperiment object
library(nnSVG)
library(scran)
library(ggplot2)
library(mclust)

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



# filter genes ------------------------------------------------------------

# Remove genes expressed in fewer than 5% of spots
x <- x[rowSums(counts(x) > 0) > (ncol(x) * 0.05), ]

#x <- computeLibraryFactors(x)
x <- logNormCounts(x)
#assayNames(x)

library(scry)
x <- devianceFeatureSelection(x, sorted = T)
plot(rowData(x)$binomial_deviance)
selected_genes <- 20 #;abline(v = selected_genes, col = 2)
x <- x[1:selected_genes,]
x <- nullResiduals(x, assay = "counts",
                   fam = c("poisson"),
                   type = c("deviance"),
                   batch = NULL)
assay(x, i = "poisson_deviance_residuals")
boxplot(t(assay(x, i = "poisson_deviance_residuals")))
abline(h = 0, lty = 2)
boxplot(t(assay(x, i = "logcounts")))

# fitting perla on logcounts (centered) -----------------------------------------------------------
selected <- t(assay(x, 
                    i = "logcounts"))
selected <- selected - mean(selected)
boxplot(selected);abline(h = 0, lty = 2)
heatmap(selected)
R <- 5*10^3
K <- 3
nstart <- 5
runs_perla <- list()
for(i in 1:nstart){
 runs_perla[[i]] <- perla(y = selected, W = W, K = K, R = R,
                       mean.penalty = c("d","cd"), burnin = 1:(R/2), seed = 1234*i) 
 runs_perla[[i]] <- perla::recover.loglikelihood(runs_perla[[i]])
}
res <- runs_perla[[which.max(lapply(runs_perla, function(y) max(y$loglik)))]]
z <- factor(apply(res$Z[,,which.max(res$loglik)], 1, function(y) which(y == 1)))
# plot the clusters
ggplot(coordinates, aes(x, y, col = z))+geom_point(cex = 4)+
  theme_bw()

# preparing the ggplot
perla_draws <- res$Mu
dimnames(perla_draws) <- list(
  Cluster = factor(1:dim(perla_draws)[1]),
  Gene = colnames(selected),
  Iteration = 1:dim(perla_draws)[3])
perla_draws <- as.data.frame.table(perla_draws)
colnames(perla_draws)[colnames(perla_draws) == "Freq"] <- "Value"
ggplot(perla_draws, aes(x = Gene, y = Value))+
  geom_boxplot()+theme_bw()+
  facet_wrap(~Cluster, nrow = 3)+
  geom_abline(intercept = 0, slope = 0, lty = 2)


# k-means -----------------------------------------------------------------

km <- factor(kmeans(selected, centers = K, nstart = 10)$cluster)
ggplot(coordinates, aes(x, y, col = km))+geom_point(cex = 4)+
  theme_bw()
mclust::adjustedRandIndex(z, km)


# GMM ---------------------------------------------------------------------

gmm <- Mclust(selected, G = K)
ggplot(coordinates, aes(x, y, col = factor(gmm$classification)))+geom_point(cex = 4)+
  theme_bw()

mclust::adjustedRandIndex(z, gmm$classification)
table(z, gmm$classification)


