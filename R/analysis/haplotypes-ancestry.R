source("R/analysis/haplotypes-env.R")

# Make a plot of anomalous regions.
svg("R/plots/pedigree_aas.svg", width = 20, height = 5)
p1 <-
  ggplot(
    pedig_aas,
    aes(x = POSITION, y = -log10(PVAL), col = CHROMOSOME)
  ) +
  geom_point(alpha = 0.5) +
  ylim(c(0, 20)) +
  geom_hline(
    aes(
      yintercept = -log10(0.01)
    ),
    col = "dodgerblue",
    linetype = 2
  ) +
  facet_grid(. ~ CHROMOSOME, scales = "free_x", space = "free_x", switch = "x") +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.spacing = unit(0, "lines"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text = element_text(size = 14),
    text = element_text(size = 14)
  ) +
  scale_color_manual(
    values = rep(
      c("darkred", "lightcoral"),
      length.out = 19
    )
  ) +
  labs(x = "Position (bp)", y = "-log10(P)")
p1
dev.off()