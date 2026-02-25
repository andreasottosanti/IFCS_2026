source("Functions/align_clusters.R")

# data selection ----------------------------------------------------------
selected_genes <- 20 
x <- x[rowData(x)$rank <= selected_genes,]

# fitting perla on logcounts (centered) -----------------------------------------------------------
selected <- t(assay(x, 
                    i = "logcounts"))
selected <- selected - mean(selected)
boxplot(selected);abline(h = 0, lty = 2)
heatmap(selected)
R <- 10^4
K <- 4
nstart <- 5
runs_perla <- list()
runtime <- c()
for(i in 1:nstart){
  runtime[i] <- system.time(runs_perla[[i]] <- perla(y = selected, W = W, K = K, R = R,
                           mean.penalty = c("d","cd"), burnin = 1:(R/2), seed = 1234*i) )
  runs_perla[[i]] <- perla::recover.loglikelihood(runs_perla[[i]])
}
mean(runtime)

# resolve the label switching
p <- merge.perla(runs_perla)
p <- remove.label.switching2(values = p)
p$relabelling$similarity

# compute the potential scale reduction factor
posterior_4d <- array(p$results$ECR$M,
                      dim = c(dim(runs_perla[[1]]$Mu)[3], nstart, K, selected_genes))
rhat_coda <- array(NA, dim = c(K, selected_genes))
par(mfrow = c(2,2))
for(k in 1:K){
  for(j in 1:selected_genes){
    samples_matrix <- posterior_4d[, , k, j]
    mcmc_list <- mcmc.list(
      lapply(1:nstart, function(l) {
        mcmc(samples_matrix[, l])
      })
    )
    rhat_coda[k, j] <- gelman.diag(
      mcmc_list,
      autoburnin = FALSE
    )$psrf[1]
    if(rhat_coda[k, j]>1.1){
      matplot(samples_matrix, type = "l", 
              main = paste("Cluster",k, "Gene",j,",",round(rhat_coda[k, j],3)))
     abline(h = 0) 
    }
  }
}
summary(as.vector(rhat_coda))
quantile(as.vector(rhat_coda), .85)
mean(as.vector(rhat_coda) <= 1.05)

z <- factor(p$relabelling$clusters[2,])
colData(x)$perla <- z

# computing the HPD intervals
res_mu <- p$results$ECR$M
perla_draws <- res_mu <- aperm(res_mu, c(2, 3, 1))
hpd <- apply(perla_draws, c(1,2), function(y) HPDinterval(mcmc(y), prob = .99))
saveRDS(res_mu, file = "Results/perlaRun_SVgenes_logcounts.RDS")

# preparing the ggplot
dimnames(perla_draws) <- list(
  Cluster = factor(1:dim(perla_draws)[1]),
  Gene = colnames(selected),
  Iteration = 1:dim(perla_draws)[3])
perla_draws <- as.data.frame.table(perla_draws)
colnames(perla_draws)[colnames(perla_draws) == "Freq"] <- "Value"
levels(perla_draws$Cluster) <- paste("Cluster", levels(perla_draws$Cluster))
perla_draws <- cbind(perla_draws,
                     leftHPD = rep(as.vector(hpd[1,,]), dim(res_mu)[3]),
                     rightHPD = rep(as.vector(hpd[2,,]), dim(res_mu)[3]))
perla_draws$Contains_zero <- factor(as.numeric(perla_draws$leftHPD < 0 & perla_draws$rightHPD > 0))
perla_draws$Cluster_plot <- ifelse(perla_draws$Contains_zero == 1,
                                   "Contains zero",
                                   as.character(perla_draws$Cluster))
ggplot(perla_draws,
       aes(x = Value, y = Gene, 
           colour = Cluster_plot)) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_boxplot(width = 0.4) +
  scale_colour_manual(
    values = c(
      "Cluster 1" = "#440154FF",
      "Cluster 2" = "#31688EFF",
      "Cluster 3" = "#35B779FF",
      "Cluster 4" = "orange",
      "Contains zero" = "lightgrey"
    ),
    breaks = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4")  # <- legend only shows these
  ) +
  theme_bw() +
  labs(y = "", x = expression(mu), colour = "") +
  theme(
    text = element_text(size = 18),
    legend.position = "bottom"
  )
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
perla_viridis[perla_viridis == "#FDE725FF"] <- "orange"
km_viridis[km_viridis == "#FDE725FF"] <- "orange"
gmm_viridis[gmm_viridis == "#FDE725FF"] <- "orange"

# saved through Rstudio
ggplot()+theme_minimal()
grid.raster(img)
grid.points(spots$x_img, -spots$y_img, pch=19, gp=gpar(col=perla_viridis), size = unit(1.5, "char"))

# kmeans
ggplot()+theme_minimal()
grid.raster(img)
grid.points(spots$x_img, -spots$y_img, pch=19, gp=gpar(col=km_viridis), size = unit(1.5, "char"))

# gmm
ggplot()+theme_minimal()
grid.raster(img)
grid.points(spots$x_img, -spots$y_img, pch=19, gp=gpar(col=gmm_viridis), size = unit(1.5, "char"))


# plot some genes (estimated)
genes <- c("Apoe", "S100a5", "Kctd12", "Sash1")
gene_means <- apply(res_mu, c(1,2), function(y) mean(mcmc(y)))[,rowData(x)$gene_names %in% genes]
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
