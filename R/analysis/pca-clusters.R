require(optparse)

option_list <- list(
  make_option(
    c("-r", "--rel"),
    type="character",
    help="Relationship matrix."
  ),
  make_option(
    c("-i", "--ids"),
    type="character",
    help="Relationship matrix IDs."
  ),
  make_option(
    c("-o", "--output"),
    type="character",
    help="Output file base."
  ),
  make_option(
    c("-n", "--npc"),
    type="numeric",
    help="Number of PCs for clustering."
  ),
  make_option(
    c("-c", "--centres"),
    type="numeric",
    help="Centres for kmeans clustering."
  )
)

set.seed(123)

opt <-
  parse_args(
    OptionParser(option_list=option_list)
  )

kmat_ids <-
  read.csv(
    opt$ids,
    col.names = c("ID"),
  )$ID

kmat <-
  read.csv(
    opt$rel,
    header = FALSE,
    sep = "\t",
    col.names = kmat_ids,
    row.names = kmat_ids
  )

pca_res <-
  eigen(kmat)

pca_clusters <-
  kmeans(
    pca_res$vectors[, 1:opt$npc],
    centers = opt$centres
  )$cluster

write.table(
  file = paste0(opt$output, ".clusters.txt"),
  x = data.frame(
    ID = kmat_ids,
    CLUSTER = pca_clusters
  ),
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  sep = "\t"
)