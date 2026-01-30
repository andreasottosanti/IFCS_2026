source("Functions/align_clusters.R")

# data selection ----------------------------------------------------------
selected_genes <- 20 
x <- x[rowData(x)$rank <= selected_genes,]
plot(rowData(x)$rank_binomial_deviance, rowData(x)$rank,
     xlab = "deviance rank", ylab = "spatial rank")


# fitting perla on logcounts (centered) -----------------------------------------------------------
selected <- t(assay(x, 
                    i = "logcounts"))
selected <- selected - mean(selected)
boxplot(selected);abline(h = 0, lty = 2)
heatmap(selected)
R <- 5*10^3
K <- 4
nstart <- 5
runs_perla <- list()
for(i in 1:nstart){
  runs_perla[[i]] <- perla(y = selected, W = W, K = K, R = R,
                           mean.penalty = c("d","cd"), burnin = 1:(R/2), seed = 1234*i) 
  runs_perla[[i]] <- perla::recover.loglikelihood(runs_perla[[i]])
}
res <- runs_perla[[which.max(lapply(runs_perla, function(y) max(y$loglik)))]]
saveRDS(res, file = "Results/perlaRun_SVgenes_logcounts.RDS")
z <- factor(apply(res$Z[,,which.max(res$loglik)], 1, function(y) which(y == 1)))
colData(x)$perla <- z

# computing the HPD intervals
hpd <- apply(res$Mu, c(1,2), function(y) HPDinterval(mcmc(y), prob = .99))

# preparing the ggplot
perla_draws <- res$Mu
dimnames(perla_draws) <- list(
  Cluster = factor(1:dim(perla_draws)[1]),
  Gene = colnames(selected),
  Iteration = 1:dim(perla_draws)[3])
perla_draws <- as.data.frame.table(perla_draws)
colnames(perla_draws)[colnames(perla_draws) == "Freq"] <- "Value"
levels(perla_draws$Cluster) <- paste("Cluster", levels(perla_draws$Cluster))
perla_draws <- cbind(perla_draws,
                     leftHPD = rep(as.vector(hpd[1,,]), dim(res$Mu)[3]),
                     rightHPD = rep(as.vector(hpd[2,,]), dim(res$Mu)[3]))
perla_draws$Contains_zero <- factor(as.numeric(perla_draws$leftHPD < 0 & perla_draws$rightHPD > 0))

ggplot(perla_draws, aes(x = Gene, y = Value, color = Contains_zero))+
  geom_boxplot()+theme_bw()+
  facet_wrap(~Cluster, nrow = 3)+
  geom_abline(intercept = 0, slope = 0, lty = 2)+
  labs(y = expression(mu))+
  scale_color_manual(values = c("0" = "black", "1" = "grey")) +
  theme( axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "none")
ggsave("Images/SVgenes_logcounts_posterior.pdf", width = 9, height = 7)


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
  geom_point(cex = 5)+theme_bw()+
  facet_wrap(~Method)

ggplot(gg_long[gg_long$Method == "perla",], aes(x, y, col = Cluster))+
  geom_point(cex = 5)+theme_bw()

img <- readPNG(source = "Images/tissue.png")
img_height <- dim(img)[1]+400
img_width  <- dim(img)[2]+640
spots <- as.data.frame(spatialCoords(x))
colnames(spots) <- c("x","y")
spots$x_img <- (spots$x - min(spots$x)) / (max(spots$x) - min(spots$x)) * img_width +220
spots$y_img <- img_height - ((spots$y - min(spots$y)) / (max(spots$y) - min(spots$y)) * img_height)-850
perla_viridis <- viridis::viridis(length(unique(colData(x)$perla)))[as.numeric(factor(colData(x)$perla))]
km_viridis <- viridis::viridis(length(unique(colData(x)$km)))[as.numeric(factor(colData(x)$km))]
gmm_viridis <- viridis::viridis(length(unique(colData(x)$gmm)))[as.numeric(factor(colData(x)$gmm))]

ggplot()+theme_minimal()
grid.raster(img)
grid.points(spots$x_img, -spots$y_img, pch=19, gp=gpar(col=perla_viridis), size = unit(1.5, "char"))

ggplot()+theme_minimal()
grid.raster(img)
grid.points(spots$x_img, -spots$y_img, pch=19, gp=gpar(col=km_viridis), size = unit(1.5, "char"))

ggplot()+theme_minimal()
grid.raster(img)
grid.points(spots$x_img, -spots$y_img, pch=19, gp=gpar(col=gmm_viridis), size = unit(1.5, "char"))


# plot some genes (estimated)
genes <- c("Apoe", "S100a5", "Kctd12", "Sash1")
gene_means <- apply(res$Mu, c(1,2), function(y) mean(mcmc(y)))[,rowData(x)$gene_names %in% genes]
colnames(gene_means) <- genes
vals <- t(gene_means[colData(x)$perla,])
vals_frame <- data.frame(values = as.vector(vals),
                         gene_names = rep(rownames(vals), ncol(x)),
                         x = rep(spatialCoords(x)[,1], each = length(genes)),
                         y = -rep(spatialCoords(x)[,2], each = length(genes)))
ggplot(vals_frame, aes(x = x, y = y, fill = values))+
  geom_point(cex = 5, shape = 21)+theme_minimal()+
  facet_wrap(~gene_names)+
  labs(fill = "", x = "", y = "")+
  theme(text = element_text(size = 18), axis.text = element_blank())+ 
  scale_fill_gradient2(midpoint=mean(vals_frame$values), mid="#999999", high="#E69F00",
                        low="#56B4E9", space ="Lab" )
ggsave(filename = "Images/SVgenes_logcounts_gene_expressions.pdf", width = 9, height = 7)

# Save results ------------------------------------------------------------

saveRDS(x, file = "Results/SVgenes_logcounts.RDS")
