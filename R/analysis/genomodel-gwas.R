set.seed(123)
source("envs/R/genomodel-env.R")

######################################################################
# GWAS using BLINK.
######################################################################
#
# 1) Load the data.
#
# If the blink data are not loaded,
# then load them.
if (
  !(
    exists("blink_snpmat") &&
      exists("blink_map") &&
      exists("blink_ids") &&
      exists("blink_snpmat_t")
  )
) {
  load_blink()
}

# Filter the field and glasshouse data
# to get the accessions in the blink
# SNP data.
data_fta24_gwasin <-
  data_fta24 |>
  filter(ID %in% blink_ids$ID)

data_glasshouse19_gwasin <-
  data_glasshouse19_filter |>
  filter(ID %in% blink_ids$ID)

######################################################################
#
# 2) Model phenotypes and extract emmeans.
#
# Model the stem lesion lengths across the field plot.
lm_fta24_gwasin <-
  lm(
    DPI15_norm ~ factor(REPLICATE) + ID, data = data_fta24_gwasin
  )

# Get the means for the GWAS.
emmeans_fta24_gwasin <-
  as.data.frame(emmeans(lm_fta24_gwasin, "ID")) |>
  mutate(ID = as.character(ID))

# Do the same for the glasshouse data.
# Note: spatial variables had no discernible
# effects, so this is just the means.
lm_glasshouse19_gwasin <-
  lm(DPI14_norm ~ 1 + ID, data = data_glasshouse19_gwasin)

emmeans_glasshouse19_gwasin <-
  as.data.frame(emmeans(lm_glasshouse19_gwasin, "ID")) |>
  mutate(ID = as.character(ID))

######################################################################

# 3) Subset the SNP and match the order of phenotypes to
# genotyped individuals.

# Get the ids from the two datasets.
fta24_ids <- emmeans_fta24_gwasin$ID
glasshouse19_ids <- emmeans_glasshouse19_gwasin$ID

# Keep individuals in the genotype data that
# are in the phenotype data.

# FTA
blink_keep_fta24 <- which(blink_ids$ID %in% fta24_ids)
blink_ids_fta24 <- blink_ids[blink_keep_fta24, ]
blink_snpmat_t_fta24 <- blink_snpmat_t[, blink_keep_fta24]

# Glasshouse
blink_keep_glasshouse19 <- which(blink_ids$ID %in% glasshouse19_ids)
blink_ids_glasshouse19 <- blink_ids[blink_keep_glasshouse19, ]
blink_snpmat_t_glasshouse19 <- blink_snpmat_t[, blink_keep_glasshouse19]

# Match the order of phenotypes to the order of genotypes.

# FTA
blink_y_fta24 <-
  emmeans_fta24_gwasin[
    match(blink_ids_fta24, emmeans_fta24_gwasin$ID),
  ] |>
  mutate(taxa = ID, FTALL = emmean) |>
  dplyr::select(taxa, FTALL)

# Glasshouse.
blink_y_glasshouse19 <-
  emmeans_glasshouse19_gwasin[
    match(blink_ids_glasshouse19, emmeans_glasshouse19_gwasin$ID),
  ] |>
  mutate(taxa = ID, GHLL = emmean) |>
  dplyr::select(taxa, GHLL)

# FTA non-winter lines

set.seed(123)
km_clusters <-
  kmeans(
    kinmat_pcs[blink_y_fta24$taxa, 1:2], centers = 3
  )$cluster

# Test plot.

test_plot_df <-
  data.frame(
    PC1 = kinmat_pcs[blink_y_fta24$taxa, 1],
    PC2 = kinmat_pcs[blink_y_fta24$taxa, 2],
    CLUSTER = factor(km_clusters)
  )

png("tmp/kinship-test-plot.png")
ggplot(test_plot_df, aes(x = PC1, y = PC2, color = CLUSTER)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(title = "K-means Clustering of Kinship PCs") +
  theme(legend.title = element_blank())
dev.off()

# Check plot before doing this.

nw_lines <-
  intersect(names(km_clusters)[km_clusters != 1], blink_ids_fta24)
blink_ids_fta24nw <-
  blink_ids_fta24[blink_ids_fta24 %in% nw_lines]

blink_y_fta24nw <-
  emmeans_fta24_gwasin |>
  filter(ID %in% blink_ids_fta24nw) |>
  mutate(taxa = ID, FTALL = emmean) |>
  dplyr::select(taxa, FTALL)

blink_y_fta24nw <-
  blink_y_fta24nw[match(nw_lines, blink_y_fta24nw$taxa), ]

blink_keep_fta24nw <- which(blink_ids$ID %in% blink_ids_fta24nw)
blink_snpmat_t_fta24nw <- blink_snpmat_t[, blink_keep_fta24nw]

# Write out the list of lines in each GWAS for later use.
gwas_lines_used <-
  bind_rows(
    blink_y_fta24 |>
      rename(SCORE = FTALL) |>
      mutate(GWAS = "fta24"),
    blink_y_glasshouse19 |>
      rename(SCORE = GHLL) |>
      mutate(GWAS = "glasshouse19"),
    blink_y_fta24nw |>
      rename(SCORE = FTALL) |>
      mutate(GWAS = "fta24nw")) |>
  dplyr::select(taxa, GWAS)

write.table(
  file = "data/gwas/out/gwas-lines-used.txt",
  x = gwas_lines_used, sep = "\t", row.names = FALSE,
  col.names = FALSE
)

######################################################################

# 4) Add PCA covariates.

# FTA
blink_y_fta24_pcs <-
  blink_y_fta24 |>
  bind_cols(kinmat_pcs[blink_y_fta24$taxa, 1:3]) |>
  rename("P1" = "...3", "P2" = "...4", "P3" = "...5")

# Glasshouse.
blink_y_glasshouse19_pcs <-
  blink_y_glasshouse19 |>
  bind_cols(kinmat_pcs[blink_y_glasshouse19$taxa, 1:3]) |>
  rename("P1" = "...3", "P2" = "...4", "P3" = "...5")

blink_y_fta24nw_pcs <-
  blink_y_fta24nw |>
  bind_cols(kinmat_pcs[blink_y_fta24nw$taxa, 1:3]) |>
  rename("P1" = "...3", "P2" = "...4", "P3" = "...5")

######################################################################

# 5) Fit GWAS models.

# Fit models with PCA covariates.

# FTA
blink_fta24_pca <-
  Blink(
    Y = blink_y_fta24_pcs[, 1:2],
    GD = blink_snpmat_t_fta24,
    GM = blink_map,
    CV = as.matrix(blink_y_fta24_pcs$P1),
    maf.threshold = 0.05,
    p.threshold = 0.01
  )

# Glasshouse.
blink_glasshouse19_pca <-
  Blink(
    Y = blink_y_glasshouse19_pcs[, 1:2],
    GD = blink_snpmat_t_glasshouse19,
    GM = blink_map,
    CV = as.matrix(blink_y_glasshouse19_pcs$P2),
    maf.threshold = 0.05,
    p.threshold = 0.01
  )

# Fit models without PCA covariates.

# FTA

blink_fta24_npca <-
  Blink(
    Y = blink_y_fta24_pcs[, 1:2],
    GD = blink_snpmat_t_fta24,
    GM = blink_map,
    maf.threshold = 0.05,
    p.threshold = 0.01
  )

# Glasshouse.

blink_glasshouse19_npca <-
  Blink(
    Y = blink_y_glasshouse19_pcs[, 1:2],
    GD = blink_snpmat_t_glasshouse19,
    GM = blink_map,
    maf.threshold = 0.05,
    p.threshold = 0.01
  )

# Fit models to the FTA data without winter lines.

blink_fta24nw_npca <-
  Blink(
    Y = blink_y_fta24nw_pcs[, 1:2],
    GD = blink_snpmat_t_fta24nw,
    GM = blink_map,
    maf.threshold = 0.05,
    p.threshold = 0.01
  )

# Write out results.
write.table(
  file = "data/gwas/out/fta24_pca_gwas.txt", 
  x = blink_fta24_pca$GWAS, sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(
  file = "data/gwas/out/glasshouse19_pca_gwas.txt", 
  x = blink_glasshouse19_pca$GWAS, sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(
  file = "data/gwas/out/fta24_npca_gwas.txt", 
  x = blink_fta24_npca$GWAS, sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(
  file = "data/gwas/out/glasshouse19_npca_gwas.txt",
  x = blink_glasshouse19_npca$GWAS, sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(
  file = "data/gwas/out/fta24nw_gwas.txt", 
  x = blink_fta24nw_npca$GWAS, sep = "\t", row.names = FALSE, quote = FALSE
)

save.image("R/sessions/genomodel-gwas.RData")

######################################################################

# 6) GWAS diagnostics.

# TODO
blink_gwas_res <-
  list(
    as.data.frame(blink_fta_pca$GWAS),
    as.data.frame(blink_glasshouse_pca$GWAS),
    as.data.frame(blink_ftanw_pca$GWAS),
    as.data.frame(blink_fta_npca$GWAS),
    as.data.frame(blink_glasshouse_npca$GWAS),
    as.data.frame(blink_ftanw_npca$GWAS)
  )

blink_qqs <-
  lapply(blink_gwas_res, function(df) {
    qq_plot_prep(df$P.value)
  })

model_names <- c(
  "FTA (PCA)", "Glasshouse (PCA)", "FTA NW (PCA)",
  "FTA (no PCA)", "Glasshouse (no PCA)", "FTA NW (no PCA)"
)

# Lambda summary table.
lambda_df <- data.frame(
  MODEL  = model_names,
  LAMBDA = sapply(blink_qqs, `[[`, 2)
)

# Combined QQ data frame for faceting.
qq_combined_df <- do.call(rbind, lapply(seq_along(blink_qqs), function(i) {
  df <- blink_qqs[[i]][[1]]
  df$MODEL  <- model_names[i]
  df$LAMBDA <- blink_qqs[[i]][[2]]
  df
})) |>
  mutate(
    MODEL = factor(MODEL, levels = model_names),
    LABEL = paste0(MODEL, "\nλ = ", round(LAMBDA, 3))
  )

# Make the qq plots.
png("R/plots/genomodel-gwas_qq.png", width = 2000, height = 3000)
ggplot(qq_combined_df, aes(x = expected, y = observed)) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, col = "red") +
  facet_wrap(~ LABEL, ncol = 3) +
  theme_classic() +
  labs(
    x = expression(Expected ~ -log[10](p)),
    y = expression(Observed ~ -log[10](p))
  ) +
  theme(axis.text = element_text(size = 11), text = element_text(size = 11))
dev.off()

save.image("R/sessions/genomodel-gwas.RData")