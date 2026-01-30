align_clusters <- function(ref, target) {
  tab <- table(ref, target)
  mapping <- apply(tab, 2, which.max)
  aligned_target <- mapping[as.character(target)]
  return(aligned_target)
}

