set.seed(123)

# Source the environment.
source("envs/R/genomodel-env.R")

# If this has been run before, it should be saved.
if (file.exists("R/sessions/genomodel-gp.RData")) {
  load("R/sessions/genomodel-gp.RData")
}

# Start by making a plot of the PCA.
svg("R/plots/genomodel-gp_pca.svg", width = 6, height = 6)
ggplot(pca_plot_df, aes(x = X1, y = X2, fill = CLUSTER)) +
  geom_point(shape = 21, size = 4, alpha = 0.9) +
  labs(
    x = paste0("PC1 (", round(var_explained[1] * 100, 2), "%)"),
    y = paste0("PC2 (", round(var_explained[2] * 100, 2), "%)")
  ) +
  theme_minimal() +
  theme(
    legend.position = "none", text = element_text(size = 15),
    axis.text = element_text(size = 15)
  ) +
  scale_fill_manual(values = c("#4E79A7", "#E15759", "#76B7B2"))
dev.off()

# Use the full kinship matrix for genomic prediction.
kinmat_gp <- kinmat

# The field data are the unfiltered data.
data_fta24_gpin <-
  data_fta24 |>
  filter(ID %in% rownames(kinmat_gp))

# The glasshouse data are also unfiltered.
data_glasshouse19_gpin <-
  data_glasshouse19_filter |>
  filter(ID %in% rownames(kinmat_gp))

# Models to estimate narrow sense heritability.
# The field data. Add type so that the model
# is not just fitting broad population structure.
base_fit_fta24 <-
  mmes(
    DPI15_norm ~ 1 + factor(REPLICATE) + factor(TYPE),
    random = ~ vsm(ism(ID), Gu = kinmat_gp),
    rcov = ~ units,
    data = data_fta24_gpin
  )

# The glasshouse data - no covariates had any impacts.
# TYPE is not included as only 6 varieties were not spring.
base_fit_glasshouse19 <-
  mmes(
    DPI14_norm ~ 1,
    random = ~ vsm(ism(ID), Gu = kinmat_gp),
    rcov = ~ units,
    data = data_glasshouse19_gpin
  )

# Save the narrow sense heritability estimates.
# These are 0.33 and 0.26, respectively.
h2_fta24 <- vpredict(base_fit_fta24, h2 ~ V1 / (V1 + V2))
h2_glasshouse19 <- vpredict(base_fit_glasshouse19, h2 ~ V1 / (V1 + V2))

####### TRAIT CORRELATIONS #######

# Genetic correlation between flowering time and lesion length.
# The long format DF has the normalised DPI15 score and flowering week,
# both scaled and in the same numeric vector. This will be matched
# by a character vector mapping scores to either trait.
data_fta24_gpin_long <-
  data_fta24_gpin |>
  mutate(
    DPI15_norm = scale(DPI15_norm),
    FLOWERING_WEEK = scale(FLOWERING_WEEK)
  ) |>
  select(
    ID, DPI15_norm, FLOWERING_WEEK, REPLICATE, TYPE
  ) |>
  pivot_longer(
    names_to = "TRAIT",
    values_to = "SCORE",
    cols = c(DPI15_norm, FLOWERING_WEEK)
  ) |>
  mutate(TRAIT = factor(TRAIT, levels = c("DPI15_norm", "FLOWERING_WEEK")))

# Needs to be ordered for sommer.
data_fta24_gpin_long <-
  data_fta24_gpin_long[order(data_fta24_gpin_long$TRAIT), ]

# Fit the same basic model but with TRAIT is a fixed / random effect.
# The traits are allowed to covary but their errors are uncorrelated.
llft_fit_fta24 <- mmes(
  SCORE ~ 1 + TRAIT + factor(REPLICATE) + factor(TYPE),
  random = ~ vsm(usm(TRAIT), ism(ID), Gu = kinmat_gp),
  rcov = ~ vsm(dsm(TRAIT), ism(units)),
  data = data_fta24_gpin_long
)

# Do exactly the same steps for the glasshouse data.
data_glasshouse19_gpin_long <-
  data_glasshouse19_gpin |>
  mutate(
    DPI14_norm = scale(DPI14_norm),
    FLOWERING_WEEK = scale(FLOWERING_WEEK)
  ) |>
  select(
    ID, DPI14_norm, FLOWERING_WEEK, STEM_WIDTH
  ) |>
  pivot_longer(
    names_to = "TRAIT",
    values_to = "SCORE",
    cols = c(DPI14_norm, FLOWERING_WEEK)
  ) |>
  mutate(
    TRAIT = factor(TRAIT, levels = c("DPI14_norm", "FLOWERING_WEEK"))
  )

data_glasshouse19_gpin_long <-
  data_glasshouse19_gpin_long[order(data_glasshouse19_gpin_long$TRAIT), ]

llft_fit_glasshouse19 <- mmes(
  SCORE ~ 1 + TRAIT,
  random = ~ vsm(usm(TRAIT), ism(ID), Gu = kinmat_gp),
  rcov = ~ vsm(dsm(TRAIT), ism(units)),
  data = data_glasshouse19_gpin_long
)

# Get genetic correlations for flowering time
# and lesion length.
# -0.37 and -0.38, respectively.
g_cor_llft_fta24 <-
  vpredict(llft_fit_fta24, h2 ~ V2 / sqrt(V1 * V3))
g_cor_llft_glasshouse19 <-
  vpredict(llft_fit_glasshouse19, h2 ~ V2 / sqrt(V1 * V3))

##################################

#### ENVIRONMENT CORRELATIONS ####

# Make a combined data frame to test correlation across
# environments for flowering and disease.
data_combined_gpin_long <-
  bind_rows(
    data_fta24_gpin |>
      mutate(
        LESION = scale(DPI15_norm),
        FT = scale(FLOWERING_WEEK)
      ) |>
      select(ID, LESION, FT, ENVIRONMENT, REPLICATE, TYPE),
    data_glasshouse19_gpin |>
      mutate(
        LESION = scale(DPI14_norm),
        FT = scale(FLOWERING_WEEK),
        REPLICATE = 3
      ) |>
      select(ID, LESION, FT, ENVIRONMENT, REPLICATE, TYPE)
  )

# Genetic correlation for lesion length between environments.
# Replicate is dropped due to convergence failure.
llenv_fit_combined <-
  mmes(

    # Environment and type as fixed.
    LESION ~ 1 + ENVIRONMENT + factor(TYPE),

    # Environment and ID are random,
    # with a usm covariance matrix for environment.
    random = ~ vsm(
      usm(ENVIRONMENT),
      ism(ID),
      Gu = kinmat_gp # Kinship matrix.
    ),

    # Environment errors are uncorrelated.
    rcov = ~ vsm(dsm(ENVIRONMENT), ism(units)),
    data = data_combined_gpin_long
  )

# Exactly as above but for flowering time.
ftenv_fit_combined <-
  mmes(
    FT ~ 1 + ENVIRONMENT + factor(TYPE),
    random = ~ vsm(usm(ENVIRONMENT), ism(ID), Gu = kinmat_gp),
    rcov = ~ vsm(dsm(ENVIRONMENT), ism(units)),
    data = data_combined_gpin_long
  )

# Get the genetic correlations between environments.
# 0.54 and 0.81, respectively.
g_cor_ll_env <- vpredict(llenv_fit_combined, g_cor ~ V2 / sqrt(V1 * V3))
g_cor_ft_env <- vpredict(ftenv_fit_combined, g_cor ~ V2 / sqrt(V1 * V3))

##################################

# Then correct lesion length and flowering time for a new
# mt menv model based on emmeans.
lm_fta24ll <- lm(
  DPI15_norm ~ 1 + factor(REPLICATE) + ID,
  data = data_fta24_gpin
)

emmeans_fta24ll <-
  as.data.frame(
    emmeans(lm_fta24ll, "ID")
  ) |>
  mutate(
    SCORE = scale(emmean),
    ENVIRONMENT = "fta24",
    TRAIT = "LESION",
    W = 1 / (SE^2)
  ) |>
  select(ID, SCORE, ENVIRONMENT, TRAIT, W)

lm_glasshouse19ll <- lm(
  DPI14_norm ~ 1 + ID,
  data = data_glasshouse19_gpin
)

emmeans_glasshouse19ll <-
  as.data.frame(
    emmeans(lm_glasshouse19ll, "ID")
  ) |>
  mutate(
    SCORE = scale(emmean),
    ENVIRONMENT = "glasshouse19",
    TRAIT = "LESION",
    W = 1 / (SE^2)
  ) |>
  select(ID, SCORE, ENVIRONMENT, TRAIT, W)

lm_fta24ft <- lm(
  FLOWERING_WEEK ~ 1 + factor(REPLICATE) + ID,
  data = data_fta24_gpin
)

emmeans_fta24ft <-
  as.data.frame(
    emmeans(lm_fta24ft, "ID")
  ) |>
  mutate(
    SCORE = scale(emmean),
    ENVIRONMENT = "fta24",
    TRAIT = "FLOWERING",
    W = 1 / (SE^2)
  ) |>
  select(ID, SCORE, ENVIRONMENT, TRAIT, W)

lm_glasshouse19ft <- lm(
  FLOWERING_WEEK ~ 1 + ID,
  data = data_glasshouse19_gpin
)

emmeans_glasshouse19ft <-
  as.data.frame(
    emmeans(lm_glasshouse19ft, "ID")
  ) |>
  mutate(
    SCORE = scale(emmean),
    ENVIRONMENT = "glasshouse19",
    TRAIT = "FLOWERING",
    W = 1 / (SE^2)
  ) |>
  select(ID, SCORE, ENVIRONMENT, TRAIT, W)

# Combine these into a dataframe for the cross validations.
emmeans_combined <-
  bind_rows(
    emmeans_fta24ll, emmeans_glasshouse19ll,
    emmeans_fta24ft, emmeans_glasshouse19ft
  ) |>
  mutate(ENVTRAIT = paste0(ENVIRONMENT, TRAIT)) |>
  arrange(ENVTRAIT) |>
  inner_join(pca_plot_df |> select(ID, CLUSTER), by = "ID") |>
  mutate(TYPE = as.factor(CLUSTER), ENVTRAIT = factor(ENVTRAIT))

######## MODEL SELECTIONS ########

# Here, I will test different 1 and 2 stage models with different
# parameterisations to determine which is the most accurate
# in just the Australian cultivars. This will be the final model.

## 1) Fit 2-step emmeans models ##

# Test data are the emmeans with cultivar scores set to NA.
cultiv_test <-
  emmeans_combined |>
  mutate(
    SCORE = ifelse(
      ID %in% aus_cultivars & ENVIRONMENT == "fta24", NA, SCORE
    )
  )

# Model 1 - no covariates, just the FTA24 data.
cultiv_fit1 <- mmes(
  SCORE ~ 1,
  random = ~ vsm(ism(ID), Gu = kinmat_gp),
  rcov = ~ units,
  data = cultiv_test |> filter(ENVTRAIT == "fta24LESION")
)

# Model 2 - add the glasshouse data as a covariate.
cultiv_fit2 <- mmes(
  SCORE ~ 1 + ENVIRONMENT,
  random = ~ vsm(usm(ENVIRONMENT), ism(ID), Gu = kinmat_gp),
  rcov = ~ vsm(dsm(ENVIRONMENT), ism(units)),
  data = cultiv_test |> filter(TRAIT == "LESION")
)

# Model 3 - add flowering time as a covariate.
cultiv_fit3 <- mmes(
  SCORE ~ 1 + TRAIT,
  random = ~ vsm(usm(TRAIT), ism(ID), Gu = kinmat_gp),
  rcov = ~ vsm(dsm(TRAIT), ism(units)),
  data = cultiv_test |> filter(ENVIRONMENT == "fta24")
)

# Model 4 - envtrait model. Fits a 4 x 4 covariance matrix to
# all environment x trait combinations.
cultiv_fit4 <- mmes(
  SCORE ~ 1 + ENVTRAIT,
  random = ~ vsm(usm(ENVTRAIT), ism(ID), Gu = kinmat_gp),
  rcov = ~ vsm(dsm(ENVTRAIT), ism(units)),
  data = cultiv_test
)

## 2) Fit 1-step rawdata models ##

# Make an ENVTRAIT df from the full data
# and mask the cultivars.
data_combined_envtrait <-
  data_combined_gpin_long |>
  pivot_longer(
    names_to = "TRAIT",
    values_to = "SCORE",
    cols = c(LESION, FT)
  ) |>
  mutate(
    ENVTRAIT = paste0(ENVIRONMENT, TRAIT),
  ) |>
  arrange(ENVTRAIT)

cultiv_test_os_envtrait <-
  data_combined_envtrait |>
  mutate(
    SCORE = ifelse(
      ID %in% aus_cultivars & ENVIRONMENT == "fta24" & TRAIT == "LESION",
      NA, SCORE
    )
  )

# Model 5 - raw data, no covariates.
cultiv_fit5 <- mmes(
  SCORE ~ 1,
  random = ~ vsm(ism(ID), Gu = kinmat_gp),
  rcov = ~ units,
  data = cultiv_test_os_envtrait |> filter(ENVTRAIT == "fta24LESION")
)

# Model 6 - raw data, with environment as a covariate.
cultiv_fit6 <- mmes(
  SCORE ~ 1 + factor(ENVIRONMENT),
  random = ~ vsm(usm(ENVIRONMENT), ism(ID), Gu = kinmat_gp),
  rcov = ~ vsm(dsm(ENVIRONMENT), ism(units)),
  data = cultiv_test_os_envtrait |> filter(TRAIT == "LESION")
)

# Model 7 - raw data, with trait as a covariate.
cultiv_fit7 <- mmes(
  SCORE ~ 1 + factor(TRAIT),
  random = ~ vsm(usm(TRAIT), ism(ID), Gu = kinmat_gp),
  rcov = ~ vsm(dsm(TRAIT), ism(units)),
  data = cultiv_test_os_envtrait |> filter(ENVIRONMENT == "fta24")
)

# Model 8 - raw data, with envtrait as a covariate.
cultiv_fit8 <- mmes(
  SCORE ~ 1 + factor(ENVTRAIT),
  random = ~ vsm(usm(ENVTRAIT), ism(ID), Gu = kinmat_gp),
  rcov = ~ vsm(dsm(ENVTRAIT), ism(units)),
  data = cultiv_test_os_envtrait
)

## 3) Get predictive ability and accuracy from the models.

cultiv_emmeans <-
  emmeans_combined |>
  filter(ENVTRAIT == "fta24LESION") |>
  slice(match(aus_cultivars, ID)) |>
  pull(SCORE)

cultiv_means <-
  data_combined_envtrait |>
  filter(ENVTRAIT == "fta24LESION") |>
  group_by(ID) |>
  summarise(MEAN = mean(SCORE, na.rm = TRUE)) |>
  slice(match(aus_cultivars, ID)) |>
  pull(MEAN)

cultiv_pred1 <- cultiv_fit1$uList[[1]][, 1][aus_cultivars]
cultiv_pred2 <- cultiv_fit2$uList[[1]][, 1][aus_cultivars]
cultiv_pred3 <- cultiv_fit3$uList[[1]][, 2][aus_cultivars]
cultiv_pred4 <- cultiv_fit4$uList[[1]][, 2][aus_cultivars]
cultiv_pred5 <- cultiv_fit5$uList[[1]][, 1][aus_cultivars]
cultiv_pred6 <- cultiv_fit6$uList[[1]][, 1][aus_cultivars]
cultiv_pred7 <- cultiv_fit7$uList[[1]][, 2][aus_cultivars]
cultiv_pred8 <- cultiv_fit8$uList[[1]][, 2][aus_cultivars]

cultiv_cor1 <- cor(cultiv_pred1, cultiv_emmeans, use = "complete.obs")
cultiv_cor2 <- cor(cultiv_pred2, cultiv_emmeans, use = "complete.obs")
cultiv_cor3 <- cor(cultiv_pred3, cultiv_emmeans, use = "complete.obs")
cultiv_cor4 <- cor(cultiv_pred4, cultiv_emmeans, use = "complete.obs")
cultiv_cor5 <- cor(cultiv_pred5, cultiv_means, use = "complete.obs")
cultiv_cor6 <- cor(cultiv_pred6, cultiv_means, use = "complete.obs")
cultiv_cor7 <- cor(cultiv_pred7, cultiv_means, use = "complete.obs")
cultiv_cor8 <- cor(cultiv_pred8, cultiv_means, use = "complete.obs")

cultiv_cor_df <-
  data.frame(
    MODEL = c(
      "1-step emmeans no covariates",
      "1-step emmeans env",
      "1-step emmeans trait",
      "1-step emmeans envtrait",
      "1-step raw no covariates",
      "1-step raw env",
      "1-step raw trait",
      "1-step raw envtrait"
    ),
    CORRELATION = c(
      cultiv_cor1, cultiv_cor2, cultiv_cor3, cultiv_cor4,
      cultiv_cor5, cultiv_cor6, cultiv_cor7, cultiv_cor8
    )
  )

##################################

# Model 9 - Fit best model (ENVIRONMENT) with 'TYPE'
cultiv_fit9 <-
  mmes(
    SCORE ~ 1 + ENVIRONMENT + factor(TYPE),
    random = ~ vsm(usm(ENVIRONMENT), ism(ID), Gu = kinmat_gp),
    rcov = ~ vsm(dsm(ENVIRONMENT), ism(units)),
    data = cultiv_test_os_envtrait |> filter(TRAIT == "LESION")
  )

# Model 10 - the same but without environment.
cultiv_fit10 <-
  mmes(
    SCORE ~ 1 + factor(TYPE),
    random = ~ vsm(ism(ID), Gu = kinmat_gp),
    rcov = ~ units,
    data = cultiv_test_os_envtrait |> filter(ENVTRAIT == "fta24LESION")
  )

lm_fta24ll_t <-
  lm(
    DPI15_norm ~ 1 + factor(REPLICATE) + factor(TYPE) + ID,
    data = data_fta24_gpin
  )

emmeans_fta24ll_t_df <-
  as.data.frame(
    emmeans(lm_fta24ll_t, "ID")
  ) |>
  mutate(
    SCORE = scale(emmean)
  )

emmeans_fta24ll_t <-
  emmeans_fta24ll_t_df |>
  slice(match(aus_cultivars, ID)) |>
  pull(SCORE)

# Get predictive abilities.
cultiv_pred9 <- cultiv_fit9$uList[[1]][, 1][aus_cultivars]
cultiv_pred10 <- cultiv_fit10$uList[[1]][, 1][aus_cultivars]

cultiv_cor9 <- cor(cultiv_pred9, emmeans_fta24ll_t, use = "complete.obs")
cultiv_cor10 <- cor(cultiv_pred10, emmeans_fta24ll_t, use = "complete.obs")

# Stratified cross-validation.
cluster_for_cv <-
  pca_plot_df |>
  filter(ID %in% lines_shared) |>
  arrange(match(ID, lines_shared)) |>
  pull(CLUSTER)

# Cross-validate single and multi env models.
# Get a list of lines in both environments.
glasshouse19_ids <-
  data_combined_envtrait |> filter(ENVTRAIT == "glasshouse19LESION") |> pull(ID)

fta24_ids <-
  data_combined_envtrait |> filter(ENVTRAIT == "fta24LESION") |> pull(ID)

# Mean kinship of each shared line to the rest of the panel.
# Used as a per-line proxy for how relatedness to the training
# set affects prediction accuracy.
lines_shared <-
  intersect(
    as.character(glasshouse19_ids),
    as.character(fta24_ids)
  )

mean_kin_panel <-
  setNames(
    vapply(lines_shared, function(id) {
      mean(kinmat_gp[id, setdiff(lines_shared, id)])
    }, numeric(1)),
    lines_shared
  )

# Repeated stratified 5-fold CV (20 repeats x 5 folds = 100 estimates).
# The outer lapply runs repeats sequentially; within each repeat the
# 5 fold-models run in parallel.
n_repeats <- 20
n_folds <- 3
plan(multisession, workers = n_folds)

cv_strat_results <-
  lapply(seq_len(n_repeats), function(rep_i) {
    fold_ids_rep <- stratified_fold_ids(cluster_for_cv, n_folds)
    future_map(
      seq_len(n_folds),
      function(fold) {
        res <- cv_function(
          fold       = fold,
          line_names = lines_shared,
          data       = data_combined_envtrait,
          kinmat     = kinmat_gp,
          fold_ids   = fold_ids_rep
        )
        res$repeat_i <- rep_i
        res
      },
      .options = furrr_options(seed = TRUE)
    )
  }) |>
  unlist(recursive = FALSE)

plan(sequential)

# Collect per-line predictions across all repeats and folds,
# attaching each line's mean kinship to the panel.
cv_strat_linedata <-
  do.call(rbind, lapply(cv_strat_results, function(x) {
    data.frame(
      repeat_i = x$repeat_i,
      fold     = x$fold,
      ID       = x$ID,
      y_gblup1 = x$y_gblup1,
      y_gblup2 = x$y_gblup2,
      y_obs    = x$y_obs,
      MEAN_KIN = mean_kin_panel[x$ID]
    )
  }))

# Per-repeat correlation (100 values → stable mean and CI).
cv_strat_summary <-
  cv_strat_linedata |>
  group_by(repeat_i) |>
  summarise(
    COR_SINGLE = cor(y_gblup1, y_obs, use = "complete.obs"),
    COR_MULTI  = cor(y_gblup2, y_obs, use = "complete.obs"),
    .groups = "drop"
  )

# Breakdown by kinship quartile: does higher mean kinship to the
# training panel improve prediction accuracy?
cv_strat_kin_summary <-
  cv_strat_linedata |>
  mutate(KIN_QUARTILE = ntile(MEAN_KIN, 4)) |>
  group_by(repeat_i, KIN_QUARTILE) |>
  summarise(
    COR_SINGLE = cor(y_gblup1, y_obs, use = "complete.obs"),
    COR_MULTI  = cor(y_gblup2, y_obs, use = "complete.obs"),
    .groups = "drop"
  ) |>
  group_by(KIN_QUARTILE) |>
  summarise(
    MEAN_COR_SINGLE = mean(COR_SINGLE, na.rm = TRUE),
    MEAN_COR_MULTI  = mean(COR_MULTI,  na.rm = TRUE),
    .groups = "drop"
  )

### PLOTS ###

# Plot the results for the Australian cultivar-based models.
svg("R/plots/genomodel-gp_predvobs-aus-cultivars.svg", width = 4, height = 4)
plot_df <- data.frame(
  ID = aus_cultivars,
  observed = cultiv_emmeans,
  predicted = cultiv_pred6
  ) |>
  filter(!is.na(observed)) |>
  inner_join(emmeans_fta24ft |>
  select(ID, SCORE), by = "ID") |>
  inner_join(varieties_key, by = "ID")

ggplot(
  plot_df,
  aes(x = predicted, y = observed, fill = SCORE)
) +
  geom_smooth(method = "lm", col = "darkred", fill = "lightgrey") +
  geom_point(size = 5, alpha = 0.8, shape = 21) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 15),
    text = element_text(size = 15)
  ) +
  geom_text_repel(
    data = plot_df |> filter(LINE_NAME %in% c(
      "YICKADEE", "SIREN", "CLANCY", "Sapphire", "ATR-Summit",
      "Beacon", "EUREKA", "AG Emblem", "DRUM", "SHIRALEE"
    )),
    aes(label = LINE_NAME),
    size = 5,
    segment.color = "black",
    segment.size = 0.5,
    min.segment.length = 0,
    box.padding = 0.5,
    point.padding = 0.5
  ) +
  scale_fill_gradient2(
    low = "#4E79A7", high = "#E15759", name = "",
    midpoint = 0, limits = c(-2, 2)
  ) +
  labs(x = "Predicted", y = "Observed")
dev.off()

# Make some plots.
fold_colors <- c(
  "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
  "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC"
)

svg("R/plots/genomodel-gp_predvobs-allmods.svg", width = 6, height = 10)

p1 <- ggplot(
  cv_fta24_predictions,
  aes(x = predicted, y = observed)
  ) +
  geom_smooth(method = "lm", col = "black", fill = "lightgrey") +
  geom_point(
    aes(fill = as.factor(fold)),
    size = 4, alpha = 0.6, shape = 21,
    show.legend = FALSE
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 15),
    text = element_text(size = 15)
  ) +
  scale_fill_manual(values = fold_colors) +
  labs(x = "Predicted", y = "Observed") +
  ggtitle("FTA single") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.margin = margin(
      t = 5, r = 15, b = 5, l = 5, unit = "pt"
    )
  )

p2 <- ggplot(
  cv_glasshouse19_predictions,
  aes(x = predicted, y = observed)
  ) +
  geom_smooth(method = "lm", col = "black", fill = "lightgrey") +
  geom_point(
    aes(fill = as.factor(fold)),
    size = 4, alpha = 0.6, shape = 21,
    show.legend = FALSE
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 15),
    text = element_text(size = 15)
  ) +
  scale_fill_manual(values = fold_colors) +
  labs(x = "Predicted", y = "Observed") +
  ggtitle("GH single") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.margin = margin(
      t = 5, r = 15, b = 5, l = 5, unit = "pt"
    )
  )

p3 <- ggplot(
  cv2_combined_predictions,
  aes(x = predicted, y = observed)
  ) +
  geom_smooth(method = "lm", col = "black", fill = "lightgrey") +
  geom_point(
    aes(fill = as.factor(fold)),
    size = 4, alpha = 0.6, shape = 21,
    show.legend = FALSE
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 15),
    text = element_text(size = 15)
  ) +
  scale_fill_manual(values = fold_colors) +
  labs(x = "Predicted", y = "Observed") +
  ggtitle("All lines FTA+GH") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.margin = margin(
      t = 5, r = 15, b = 5, l = 5, unit = "pt"
    )
  )

p4 <- ggplot(
  cv2_combined_predictions,
  aes(x = predicted_s, y = observed)
  ) +
  geom_smooth(method = "lm", col = "black", fill = "lightgrey") +
  geom_point(
    aes(fill = as.factor(fold)),
    size = 4, alpha = 0.6, shape = 21,
    show.legend = FALSE
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 15),
    text = element_text(size = 15)
  ) +
  scale_fill_manual(values = fold_colors) +
  labs(x = "Predicted", y = "Observed") +
  ggtitle("All lines FTA") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.margin = margin(
      t = 5, r = 15, b = 5, l = 5, unit = "pt"
    )
  )

p5 <- ggplot(
  cv1_combined_predictions,
  aes(x = predicted, y = observed)
  ) +
  geom_smooth(method = "lm", col = "black", fill = "lightgrey") +
  geom_point(
    aes(fill = as.factor(fold)),
    size = 4, alpha = 0.6, shape = 21,
    show.legend = FALSE
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 15),
    text = element_text(size = 15)
  ) +
  scale_fill_manual(values = fold_colors) +
  labs(x = "Predicted", y = "Observed") +
  ggtitle("Shared lines FTA+GH") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.margin = margin(
      t = 5, r = 15, b = 5, l = 5, unit = "pt"
    )
  )

p6 <- ggplot(
  cv1_combined_predictions,
  aes(x = predicted_s, y = observed)
  ) +
  geom_smooth(method = "lm", col = "black", fill = "lightgrey") +
  geom_point(
    aes(fill = as.factor(fold)),
    size = 4, alpha = 0.6, shape = 21,
    show.legend = FALSE
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 15),
    text = element_text(size = 15)
  ) +
  scale_fill_manual(values = fold_colors) +
  labs(x = "Predicted", y = "Observed") +
  ggtitle("Shared lines FTA") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.margin = margin(
      t = 5, r = 15, b = 5, l = 5, unit = "pt"
    )
  )

plot_grid(
  p1, p2, p3, p4, p5, p6, ncol = 2,
  align = 'hv'
  )

dev.off()

combined_plot_df <-
  cv1_combined_summary |>
  pivot_longer(
    cols = c(3,5),
    names_to = "MODEL",
    values_to = "CORRELATION"
  ) |> mutate(DATA = "Shared lines") |>
  rbind(
    cv2_combined_summary |>
    pivot_longer(
      cols = c(3,5),
      names_to = "MODEL",
      values_to = "CORRELATION"
    ) |> mutate(DATA = "All lines")
  ) |>
  group_by(MODEL, DATA) |>
  summarise(
    MEDIAN_COR = median(CORRELATION),
    IQR_COR = IQR(CORRELATION),
    .groups = "drop"
  ) |>
  ungroup()

svg("R/plots/genomodel-gp_cor-improve.svg")
p7 <- ggplot(
    combined_plot_df,
    aes(x = MODEL, fill = DATA, y = MEDIAN_COR)
  ) +
  geom_bar(
    stat = "identity",
    position = "dodge",
    col = "black"
  ) +
  geom_errorbar(
    aes(ymin = MEDIAN_COR - IQR_COR, ymax = MEDIAN_COR + IQR_COR),
    width = 0.1,
    position = position_dodge(0.9)
  ) +
  theme_classic() +
  scale_fill_manual(values = c("#4E79A7", "#F28E2B"), name = "") +
  scale_x_discrete(labels = c("FTA+GH", "FTA")) +
  labs(x = "", y = "Median Pearson") +
  theme(axis.text = element_text(size = 15), text = element_text(size = 15))
p7
dev.off()

##### WORK IN PROGRESS ###
# Here, I want to make an input file for helium with gblups
# instead of lesion lengths for phenotypes.
# Then, I want to assess correlation between relationship
# to different founders and level of resistance.
svg("R/plots/genomodel-gp_program-gblups.svg", width = 6, height = 4)
japvcan_gblup <-
  ggplot(
    helifull_info,
    aes(x = TYPE, y = GBLUP, fill = TYPE)
  ) +
  geom_boxplot(varwidth = TRUE, fill = "white") +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.8) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12),
    text = element_text(size = 14),
    legend.position = "none"
  ) +
  scale_fill_manual(
    values = c("darkblue", "orange", "deepskyblue", "darkolivegreen1")
  ) +
  xlab("") +
  ylab("Estimated breeding value")
japvcan_gblup
dev.off()

# The helium input file with data from just the Shiga program.
helishiga <-
  read.csv("heliinputshiga.txt", header = TRUE, sep = "\t", skip = 1) |>
  mutate(
    ID = gsub(" $", "", ID),
    P1 = gsub(" $", "", P1),
    P2 = gsub(" $", "", P2)
  )

aus_cultiv <- varieties_key_pedigree |>
  filter(TYPE != "founder" & ID %in% rownames(kinmat)) |>
  pull(ID)

helifull_ids <-
  helifull_info |>
  filter(!is.na(GBLUP)) |>
  select(ID, TYPE) |>
  rename("HELIID" = ID) |>
  left_join(varieties_key, by = c("HELIID" = "SIMPLE_NAME"))

pca_heli <-
  pca_plot_df |>
  inner_join(helifull_ids, by = "ID")

svg("R/plots/genomodel-gp_heli-pca.svg", width = 5, height = 5)
pca_heli_plot <-
  ggplot(
    pca_heli,
    aes(x = X1, y = X2, fill = TYPE)
  ) +
  geom_point(size = 5, shape = 21) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12),
    text = element_text(size = 14),
    legend.position = "right"
  ) +
  scale_fill_manual(
    values = c("darkblue", "orange", "deepskyblue", "darkolivegreen1"),
    name = "Program"
  ) +
  labs(x = "PC1", y = "PC2")
pca_heli_plot
dev.off()

group_centroids <- pca_heli |>
  group_by(TYPE) |>
  summarise(across(starts_with("X"), mean))

# Distance of each individual from, e.g., the Japanese centroid
japanese_centroid <-
  group_centroids |>
  filter(TYPE == "Japanese") |>
  select(starts_with("X"))

pc_dist_japanese <-
  pca_heli |>
  mutate(
    DIST_JP = sqrt(
      (X1 - japanese_centroid$X1)^2 +
        (X2 - japanese_centroid$X2)^2 +
        (X3 - japanese_centroid$X3)^2 +
        (X4 - japanese_centroid$X4)^2
    ),
    DIST_CA = sqrt(
      (X1 - group_centroids$X1[group_centroids$TYPE == "Canadian"])^2 +
        (X2 - group_centroids$X2[group_centroids$TYPE == "Canadian"])^2 +
        (X3 - group_centroids$X3[group_centroids$TYPE == "Canadian"])^2 +
        (X4 - group_centroids$X4[group_centroids$TYPE == "Canadian"])^2 +
        (X5 - group_centroids$X5[group_centroids$TYPE == "Canadian"])^2
    )
  ) |>
  filter(
    TYPE == "Australian"
  ) |>
  select(ID, HELIID, DIST_JP, DIST_CA) |>
  left_join(helifull_info, by = c("HELIID" = "ID")) |>
  select(ID, GBLUP, DIST_JP, DIST_CA) |>
  pivot_longer(
    cols = c(DIST_JP, DIST_CA),
    names_to = "DIST_TYPE",
    values_to = "DISTANCE"
  ) |>
  mutate(
    DIST_TYPE = ifelse(
      DIST_TYPE == "DIST_JP", "Japanese", "Canadian"
    )
  )

svg("R/plots/genomodel-gp_japvblup.svg", width = 5, height = 5)
japvblup_plot <-
  ggplot(
    pc_dist_japanese,
    aes(x = DISTANCE, y = GBLUP, fill = DIST_TYPE)
  ) +
  geom_smooth(method = "lm", col = "darkred", fill = "lightgrey") +
  geom_point(size = 2, shape = 21) +
  facet_wrap(~ DIST_TYPE, ncol = 1) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12),
    text = element_text(size = 14),
    legend.position = "right"
  ) +
  scale_fill_manual(
    values = c("orange", "deepskyblue"),
    name = "Program"
  ) +
  labs(x = "Distance", y = "GBLUP")
japvblup_plot
dev.off()

pc_dist_japanese |>
  group_by(DIST_TYPE) |>
  summarise(
    COR = cor(DISTANCE, GBLUP, use = "complete.obs"),
    P_VALUE = cor.test(DISTANCE, GBLUP)$p.value,
    N = sum(!is.na(DISTANCE) & !is.na(GBLUP))
  )

aus_founders <-
  c(
    "V30SR",
    "Chikuzen_verkko",
    "Chisaya_verkko",
    "Haya_verkko",
    "Mutu_flye",
    "Bronowski_flye"
    )

# aus_founders <- varieties_key_pedigree |>
#  filter(TYPE == "founder" & ID %in% rownames(kinmat)) |>
#  pull(ID)

# aus_founders <- c(aus_founders, c("V53SR", "V24SR"))

aus_cultiv_blups <-
  rowMeans(llenv_fit_combined$uList[[1]][aus_cultiv, 1:2])

all_blups <-
  data.frame(
    rowMeans(llenv_fit_combined$uList[[1]][, 1:2])
  ) |>
  mutate(
    ID = names(llenv_fit_combined$uList[[1]][, 1])
  )

aus_cultiv_frel <-
  kinmat[aus_cultiv, aus_founders] |>
  as.data.frame() |>
  mutate(BLUP = aus_cultiv_blups, ID = names(aus_cultiv_blups)) |>
  pivot_longer(
    cols = -c(ID, BLUP),
    names_to = "FOUNDER",
    values_to = "RELATEDNESS"
  ) |>
  left_join(
    varieties_key_pedigree |> select(ID, LINE_NAME),
    by = c("ID" = "ID")
  ) |>
  left_join(
    varieties_key_pedigree |> select(ID, LINE_NAME),
    by = c("FOUNDER" = "ID")
  ) |>
  rename(CULTIVAR = LINE_NAME.x, FOUNDER_NAME = LINE_NAME.y)

# Compute per-founder Pearson correlation between genomic relatedness and GBLUP.
aus_cultiv_frel_cor <-
  aus_cultiv_frel |>
  group_by(FOUNDER, FOUNDER_NAME) |>
  summarise(
    COR     = cor(RELATEDNESS, BLUP, use = "complete.obs"),
    P_VALUE = cor.test(RELATEDNESS, BLUP)$p.value,
    N       = sum(!is.na(RELATEDNESS) & !is.na(BLUP)),
    .groups = "drop"
  ) |>
  mutate(
    SIG   = P_VALUE < 0.05,
    LABEL = ifelse(is.na(FOUNDER_NAME), FOUNDER, FOUNDER_NAME)
  ) |>
  arrange(COR)

# Bar plot of per-founder correlations, coloured by significance.
svg("R/plots/genomodel-gp_founder-relatedness-cor.svg", width = 8, height = 5)
ggplot(
  aus_cultiv_frel_cor,
  aes(x = reorder(LABEL, COR), y = COR, fill = SIG)
) +
  geom_col(col = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(
    values = c("TRUE" = "#E15759", "FALSE" = "#BAB0AC"),
    name   = "p < 0.05"
  ) +
  labs(x = "Founder", y = "Spearman r (relatedness vs GBLUP)") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    text        = element_text(size = 14)
  )
dev.off()

# Scatter plots of relatedness vs GBLUP, faceted by founder.
# Founders with significant correlations are highlighted.
sig_founders <- aus_cultiv_frel_cor |> filter(SIG) |> pull(FOUNDER)

svg(
  "R/plots/genomodel-gp_founder-relatedness-scatter.svg",
  width  = 4 * ceiling(sqrt(length(sig_founders))),
  height = 4 * ceiling(length(sig_founders) / ceiling(sqrt(length(sig_founders))))
)
aus_cultiv_frel |>
  filter(FOUNDER %in% sig_founders) |>
  mutate(LABEL = ifelse(is.na(FOUNDER_NAME), FOUNDER, FOUNDER_NAME)) |>
  ggplot(aes(x = RELATEDNESS, y = BLUP)) +
    geom_point(size = 2, alpha = 0.7, shape = 21, fill = "#4E79A7") +
    geom_smooth(method = "lm", col = "darkred", fill = "lightgrey") +
    facet_wrap(~ LABEL, scales = "free_x") +
    labs(
      x = "Genomic relatedness to founder",
      y = "GBLUP (resistance)"
    ) +
    theme_classic() +
    theme(
      axis.text = element_text(size = 11),
      text      = element_text(size = 13),
      strip.text = element_text(size = 11, face = "bold")
    )
dev.off()

# ---- Multivariate founder contribution ----
# Because cultivars related to one founder tend to be related to several
# others (shared ancestry), univariate correlations confound contributions.
# OLS multiple regression gives partial coefficients (each founder controlling
# for all others); ridge regression stabilises these under multicollinearity.

# Wide format: one row per cultivar, one column per founder.
aus_cultiv_frel_wide <-
  aus_cultiv_frel |>
  mutate(LABEL = ifelse(is.na(FOUNDER_NAME), FOUNDER, FOUNDER_NAME)) |>
  pivot_wider(
    id_cols     = c(ID, BLUP),
    names_from  = LABEL,
    values_from = RELATEDNESS
  )

blup_vec    <- aus_cultiv_frel_wide$BLUP
founder_mat <- aus_cultiv_frel_wide |> select(-ID, -BLUP) |> as.matrix()

# 1. Inter-founder Spearman correlation matrix — quantifies collinearity.
founder_cor_mat <-
  cor(founder_mat, method = "spearman", use = "pairwise.complete.obs")

svg("R/plots/genomodel-gp_founder-cor-heatmap.svg", width = 7, height = 6)
founder_cor_mat |>
  as.data.frame() |>
  tibble::rownames_to_column("FOUNDER_A") |>
  pivot_longer(-FOUNDER_A, names_to = "FOUNDER_B", values_to = "COR") |>
  mutate(
    FOUNDER_A = factor(FOUNDER_A, levels = rownames(founder_cor_mat)),
    FOUNDER_B = factor(FOUNDER_B, levels = rownames(founder_cor_mat))
  ) |>
  ggplot(aes(x = FOUNDER_A, y = FOUNDER_B, fill = COR)) +
    geom_tile(col = "white") +
    scale_fill_gradient2(
      low      = "#4E79A7", mid = "white", high = "#E15759",
      midpoint = 0, limits = c(-1, 1), name = "Spearman r"
    ) +
    labs(x = NULL, y = NULL, title = "Inter-founder relatedness correlations") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
      axis.text.y = element_text(size = 11),
      text        = element_text(size = 13),
      plot.title  = element_text(hjust = 0.5)
    )
dev.off()

# 2. OLS multiple regression: standardised BLUP ~ all founder relatedness.
# Standardising makes coefficients directly comparable across founders.
founder_mat_sc <- scale(founder_mat)
blup_sc        <- as.numeric(scale(blup_vec))

# Use make.names() in case any founder labels contain spaces or punctuation.
colnames(founder_mat_sc) <- make.names(colnames(founder_mat_sc))

ols_df       <- as.data.frame(founder_mat_sc)
ols_df$BLUP  <- blup_sc

ols_fit  <- lm(BLUP ~ ., data = ols_df)
ols_summ <- summary(ols_fit)$coefficients[-1, , drop = FALSE]

ols_coef_df <- data.frame(
  LABEL    = colnames(founder_mat_sc),
  ESTIMATE = ols_summ[, "Estimate"],
  SE       = ols_summ[, "Std. Error"],
  P_VALUE  = ols_summ[, "Pr(>|t|)"],
  stringsAsFactors = FALSE
) |>
  mutate(SIG = P_VALUE < 0.05) |>
  arrange(ESTIMATE)

# 3. Ridge regression with GCV-selected lambda.
# More stable than OLS when founders are correlated.
ridge_seq <- seq(0, 100, 0.1)
ridge_fit <- MASS::lm.ridge(BLUP ~ ., data = ols_df, lambda = ridge_seq)
best_lam  <- ridge_fit$lambda[which.min(ridge_fit$GCV)]

ridge_coef_vec <- as.numeric(
  coef(MASS::lm.ridge(BLUP ~ ., data = ols_df, lambda = best_lam))
)

ridge_coef_df <- data.frame(
  LABEL    = colnames(founder_mat_sc),
  ESTIMATE = ridge_coef_vec[-1],          # drop intercept
  stringsAsFactors = FALSE
)

# 4. Side-by-side coefficient plot: OLS vs Ridge.
founder_order <- ols_coef_df$LABEL[order(ols_coef_df$ESTIMATE)]

coef_combined <-
  bind_rows(
    ols_coef_df   |> select(LABEL, ESTIMATE, SIG) |> mutate(MODEL = "OLS (partial)"),
    ridge_coef_df |> mutate(SIG = FALSE,           MODEL = "Ridge (GCV)")
  ) |>
  mutate(LABEL = factor(LABEL, levels = founder_order))

svg("R/plots/genomodel-gp_founder-multivar-coef.svg", width = 11, height = 5)
ggplot(
  coef_combined,
  aes(x = LABEL, y = ESTIMATE, fill = SIG)
) +
  geom_col(col = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ MODEL, scales = "free_y") +
  scale_fill_manual(
    values = c("TRUE" = "#E15759", "FALSE" = "#BAB0AC"),
    name   = "OLS p < 0.05"
  ) +
  labs(x = "Founder", y = "Standardised coefficient (BLUP ~ relatedness)") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    text        = element_text(size = 13),
    strip.text  = element_text(size = 12, face = "bold")
  )
dev.off()

colnames(all_blups) <- c("GBLUP", "ID")

write.table(
  file = "doc/info/genomodel-gp_glubs.txt",
  x = all_blups,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)
