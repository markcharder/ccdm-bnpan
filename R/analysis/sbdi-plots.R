source("R/analysis/genomodel-env.R")
data_combined_filter_mean <-
  data_combined_filter |>
  group_by(ID, ENVIRONMENT) |>
  summarise(SCORE = median(SCORE), .groups = "drop")

p1_df <-
  data_fta24

p1 <-
  ggplot(
    p1_df,
    aes(x = SCORE)
  ) +
  geom_histogram(
    aes(y = after_stat(count)),  # Changed from density
    bins = 25,
    fill = "steelblue",
    color = "black"
  ) +
  stat_function(
    fun = function(x) {
      dnorm(x, 
            mean = mean(p1_df$SCORE, na.rm = TRUE), 
            sd = sd(p1_df$SCORE, na.rm = TRUE)) * 
        nrow(p1_df) * 
        diff(range(p1_df$SCORE, na.rm = TRUE)) / 25  # bin_width
    },
    color = "darkred",
    linewidth = 1
  ) +
  theme_minimal() +
  labs(x = "Lesion length (mm)", y = "Number of lines") +
  theme(
    axis.text = element_text(size = 15),
    text = element_text(size = 15)
  )

# Aggregate reps (mean score per line per environment)
emmeans_wide <- emmeans_combined |>
  pivot_wider(id_cols = ID, names_from = ENVIRONMENT, values_from = emmean)

cor_matrix <- cor(emmeans_wide[-1], use = "pairwise.complete.obs")
cor_matrix

aus_breeding_lines <-
  read.csv("doc/info/varieties_key_pedigree.txt", sep = "\t") |>
  filter(ID %in% data_fta24_filter$ID)

founder_ids <- aus_breeding_lines |>
  filter(TYPE == "founder") |>
  pull(ID)

cultivar_ids <- aus_breeding_lines |>
  filter(TYPE != "founder" & ID %in% names(base_fit_fta24$u[,1])) |>
  pull(ID)

founder_rel <- as.data.frame(kinmat_cultiv[, founder_ids])
founder_rel$ID <- rownames(kinmat_cultiv)

blups_fta24 <-
  data.frame(
    ID = cultivar_ids,
    BLUP = base_fit_fta24$u[cultivar_ids, 1]
  )

plot_data <- founder_rel |>
  pivot_longer(
    cols = all_of(founder_ids),
    names_to = "Founder",
    values_to = "Relatedness"
  )

plot_data <- plot_data |>
  inner_join(blups_fta24, by = "ID") |>
  filter(ID %in% cultivar_ids) |>
  inner_join(
    aus_breeding_lines |>
    mutate(Founder = ID),
    by = "Founder"
  )

svg("R/plots/blup_relatedness_fta24.svg")
ggplot(plot_data, aes(x = Relatedness, y = BLUP)) +
  geom_point(alpha = 0.4, color = "darkblue") +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  facet_wrap(~ LINE_NAME, scales = "free_x") +
  theme_minimal() +
  labs(
    x = "Proportion of Genome Shared with Founder (GRM)",
    y = "Disease Score (Lower = More Resistant)"
  )
dev.off()

p2 <- ggpairs(
  emmeans_wide |> mutate(field = fta24, glasshouse = glasshouse19, hoophouse = hh17), 
  columns = c("field", "glasshouse", "hoophouse"),
  lower = list(
    continuous = wrap(
      "smooth",
      method = "lm",
      se = FALSE,
      fill = "cadetblue",
      shape = 21,
      size = 3,
      alpha = 0.8)),
  upper = list(continuous = wrap("cor", size = 8))
  ) +
  theme_classic() +
  theme(axis.text = element_text(size = 15), text = element_text(size = 15))

svg("R/plots/lesion_lengths_histogram_fta24.svg")
p1
dev.off()

svg("R/plots/correlations_environments.svg")
p2
dev.off()

