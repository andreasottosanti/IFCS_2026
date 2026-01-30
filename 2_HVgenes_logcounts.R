source("Functions/align_clusters.R")

# data selection ----------------------------------------------------------
selected_genes <- 100
x <- x[rowData(x)$rank_binomial_deviance <= selected_genes,]
plot(rowData(x)$rank_binomial_deviance, rowData(x)$rank,
     xlab = "deviance rank", ylab = "spatial rank")


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
saveRDS(res, file = "Results/perlaRun_HVgenes_logcounts.RDS")
z <- factor(apply(res$Z[,,which.max(res$loglik)], 1, function(y) which(y == 1)))
colData(x)$perla <- z

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

set.seed(123)
km <- factor(kmeans(selected, centers = K, nstart = 10)$cluster)
colData(x)$km <- align_clusters(ref = colData(x)$perla, target = factor(km))


# GMM ---------------------------------------------------------------------

set.seed(123)
gmm <- Mclust(selected, G = K)
colData(x)$gmm <- align_clusters(ref = colData(x)$perla, target = factor(gmm$classification))


# Plot --------------------------------------------------------------------

gg <- as.data.frame(cbind(spatialCoords(x), colData(x)[,c("perla","km","gmm")]))
gg_long <- reshape(gg, 
                   varying = c("perla", "km", "gmm"), # Le colonne da "srotolare"
                   v.names = "Cluster",               # Come chiamare la colonna dei valori
                   timevar = "Method",                # Come chiamare la colonna dei nomi
                   times = c("perla", "km", "gmm"),   # I nomi da assegnare
                   direction = "long")                # Direzione della trasformazione
ggplot(gg_long, aes(x, y, col = Cluster))+
  geom_point(cex = 3)+theme_bw()+
  facet_wrap(~Method)


# Save results ------------------------------------------------------------

saveRDS(x, file = "Results/HVgenes_logcounts.RDS")
