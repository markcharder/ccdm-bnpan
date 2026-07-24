source("envs/R/pedigree-env.R")

# Make a heatmap of the kinship matrix.
svg("R/plots/pedigree_heatmap.svg")

# Get ordered kinship matrix.
kinmat_founder_ord <-
  kinmat_founder[pedigree_key$ID, pedigree_key$ID]

# Name rows based on the common variety names.
rownames(kinmat_founder_ord) <-
  colnames(kinmat_founder_ord) <-
  pedigree_key$LINE_NAME

# Calculate euclidean distance.
g_diag <- diag(kinmat_founder_ord)
n <- nrow(kinmat_founder_ord)
gii <- matrix(g_diag, nrow = n, ncol = n, byrow = FALSE)
gjj <- matrix(g_diag, nrow = n, ncol = n, byrow = TRUE)
formal_distance <- sqrt(gii + gjj - 2 * kinmat_founder_ord)
formal_dist_obj <- as.dist(formal_distance)
dist_matrix_square <- data.frame(as.matrix(formal_dist_obj))

# Plot a heatmap.
p1 <- pheatmap(
  dist_matrix_square,
  clustering_distance_rows = formal_dist_obj,
  clustering_distance_cols = formal_dist_obj,
  clustering_method = "average",
  show_rownames = TRUE,
  show_colnames = FALSE,
  color = colorRampPalette(c("navy", "white", "firebrick"))(50)
)
p1
dev.off()

svg("R/plots/pedigree_pca-ve_combined.svg", width = 5, height = 10)

# Extract principal components from a PCA of the kinship matrix.
pedig_pca <-
  eigen(kinmat_founder)

pedig_pcs <-
  data.frame(pedig_pca$vectors[, 1:10])

pedig_pcs$ID <- rownames(kinmat_founder)

pedig_pcs <-
  pedig_pcs |>
  left_join(pedigree_key, by = "ID")

variance_explained <- (pedig_pca$values / sum(pedig_pca$values)) * 100

# Plot PCA and variance explained.
p2 <-
  ggplot(
    data = data.frame(
    PC = 1:10,
    VE = variance_explained[1:10]),
    aes(x = PC, y = VE)
  ) +
  geom_bar(stat = "identity", col = "black", fill = "deepskyblue") +
  geom_text(aes(label = sprintf("%.1f%%", VE)), vjust = -0.5) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 15),
    text = element_text(size = 15)
  ) +
  labs(x = "PC", y = "Variance explained (%)")

p3 <-
  ggplot(pedig_pcs, aes(x = X1, y = X2, fill = TYPE, label = LINE_NAME)) +
  geom_point(shape = 21, size = 4, alpha = 0.8) +
  geom_text_repel(size = 4, max.overlaps = 30) +
  labs(
    x = "PC1",
    y = "PC2",
    fill = "Type"
  ) +
  scale_fill_manual(
    values = c("darkolivegreen", "firebrick"),
    name = ""
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 15),
    text = element_text(size = 15),
    legend.position = "bottom"
  )

plot_grid(p2, p3, ncol = 1, rel_heights = c(1, 1.5))

dev.off()