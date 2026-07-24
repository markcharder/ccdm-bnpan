# Environment setup.
require(this.path)
require(nlme)
require(ggplot2)
require(emmeans)
require(dplyr)
require(tidyr)
require(Hmisc)
require(SpATS)
require(lme4)
require(sommer)
require(BEDMatrix)
require(ggraph)
require(tidygraph)
require(data.table)
require(bestNormalize)
require(cowplot)
require(GGally)
library(future)
library(furrr)
require(dplyr)
require(data.table)
require(vcfR)
require(reshape2)

# QQ-plot function.
qq_plot_prep <- function(p_obs) {

  plen <- length(p_obs)
  pexp <- (1:plen) / (plen + 1)

  # Genomic inflation factor
  lambda <- median(qchisq(1 - p_obs, df = 1)) / qchisq(0.5, df = 1)

  # QQ plot
  qq_df <- data.frame(
    expected = -log10(pexp),
    observed = -log10(sort(p_obs))
  )

  # Return a list.
  list(qq_df, lambda)

}

# Stratified fold generation function.
stratified_fold_ids <- function(clusters, n_folds) {
  fold_vec <- integer(length(clusters))
  for (cl in unique(clusters)) {
    idx <- which(clusters == cl)
    fold_vec[idx] <- sample(rep(seq_len(n_folds), length.out = length(idx)))
  }
  fold_vec
}

# Function for genomic prediction cross-validation.
cv_function <-

  function(
    fold, line_names,
    data, kinmat, fold_ids
  ) {

    require(dplyr)
    lines_mask <- line_names[fold_ids == fold]

    data_cv <-
      data |>
      dplyr::mutate(
        SCORE = ifelse(
          ID %in% lines_mask & ENVIRONMENT == "fta24", NA, SCORE
        )
      )

    # Model 1 - no covariates, just the FTA24 data.
    cv1_fit <- mmes(
      SCORE ~ 1,
      random = ~ vsm(ism(ID), Gu = kinmat),
      rcov = ~ units,
      data = data_cv |> filter(ENVTRAIT == "fta24LESION")
    )

    # Model 2 - add the glasshouse data as a covariate.
    cv2_fit <- mmes(
      SCORE ~ 1 + ENVIRONMENT,
      random = ~ vsm(usm(ENVIRONMENT), ism(ID), Gu = kinmat),
      rcov = ~ vsm(dsm(ENVIRONMENT), ism(units)),
      data = data_cv |> filter(TRAIT == "LESION")
    )

    y_obs   <- data |>
      dplyr::filter(ENVTRAIT == "fta24LESION") |>
      dplyr::group_by(ID) |>
      dplyr::summarise(
        MEAN = mean(SCORE, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::filter(ID %in% lines_mask) |>
      dplyr::arrange(match(ID, lines_mask)) |>
      dplyr::pull(MEAN)

    y_gblup1 <- cv1_fit$uList[[1]][, 1][lines_mask]
    y_gblup2 <- cv2_fit$uList[[1]][, 1][lines_mask]

    list(
      fold = fold,
      n_test = length(lines_mask),
      y_gblup1 = y_gblup1,
      y_gblup2 = y_gblup2,
      y_obs = y_obs,
      ID = lines_mask
    )
  }

# A function to read in the binary GCTA GRM file
# make it readable as a kinship matrix.
read_grm_bin <-

  function(
    prefix, all_n = FALSE, size = 4
  ) {

    sum_i <- function(i){
      sum(1:i)
    }

    # The bin file.
    bin_file_name <- paste(prefix, ".grm.bin", sep = "")

    # The 'N' file.
    n_file_name <- paste(prefix, ".grm.N.bin", sep = "")

    # The ID file.
    id_file_name <- paste(prefix, ".grm.id", sep = "")

    # The ids.
    id <- read.table(id_file_name)

    # The number of individuals.
    n <- dim(id)[1]

    # Read the bin file.
    bin_file <- file(bin_file_name, "rb")

    # Make a grm from the bin file.
    grm <- readBin(
      # n individuals * (n + 1) / 2
      # because the GRM is stored as a
      # lower triangular matrix.
      bin_file, n = n * (n + 1) / 2,
      what = numeric(0), size = size
    )

    close(bin_file)

    # The n file.
    n_file <- file(n_file_name, "rb")

    # If all_n is true, read all the values.
    if (all_n == TRUE) {

      n_read <- readBin(
        n_file,
        n = n * (n + 1) / 2, what = numeric(0), size = size
      )

    } else {

      n_read = readBin(n_file, n = 1, what = numeric(0), size = size)

    }

    close(n_file)

    i <- sapply(1:n, sum_i)
    list(diag = grm[i], off = grm[-i], id = id, N = n_read)

  }

grm_to_matrix <- function(grm) {

  # The number of individuals.
  n <- length(grm$diag)

  # Square n x n matrix of zeros.
  mat <- matrix(0, n, n)

  # Start k at 1.
  k <- 1

  # For i up to 2 x n.
  for (i in 2:n) {

    # Row i, columns 1 to i - 1 = the offset[k:(k + i - 2)].
    mat[i, 1:(i - 1)] <- grm$off[k:(k + i - 2)]

    # Set k to k + i - 1.
    k <- k + i - 1
  }

  mat <- mat + t(mat)
  diag(mat) <- grm$diag
  rownames(mat) <- grm$id$V2
  colnames(mat) <- grm$id$V2
  mat

}

# Read the grm and make the kinship matrices for genomic prediction
# and PCA.
grm <- read_grm_bin("data/variants/kinship/merge_ldfiltgp_grm")
kinmat <- grm_to_matrix(grm)
grm_pca <- read_grm_bin("data/variants/kinship/merge_ldfiltpca_grm")
kinmat_pca <- grm_to_matrix(grm_pca)

# Run the PCA and assess population stratification.
kinmat_pca_eig <- eigen(kinmat_pca)
kinmat_pcs <- kinmat_pca_eig$vectors[, 1:10]
rownames(kinmat_pcs) <- rownames(kinmat_pca)
var_explained <- kinmat_pca_eig$values / sum(kinmat_pca_eig$values)
k <- 3
clusters <- kmeans(kinmat_pcs[, 1:2], centers = k)$cluster
cluster_df <-
  data.frame(
    ID = rownames(kinmat_pcs),
    CLUSTER = as.factor(clusters)
  )

# Make a data frame for plotting the PCA, for convenience.
pca_plot_df <-
  cluster_df |>
  inner_join(
    data.frame(kinmat_pcs) |> mutate(ID = rownames(kinmat_pcs)),
    by = "ID"
  )

# Variety IDs - maps the bcf / kinship IDs to the line names in the AGG.
varieties_key <-
  read.csv(
    "doc/info/varieties_key.txt",
    header = TRUE, sep = "\t"
  )

# Get Australian cultivar IDs from the info file.
aus_cultivars <-
  read.csv("doc/info/varieties_key_pedigree.txt", sep = "\t") |>
  filter(TYPE != "founder") |>
  pull(ID)

# Read the FTA24 data and add the TYPE spring / winter column.
data_fta24 <-
  read.csv(
    "data/gp/inputs/data_fta24.txt",
    sep = "\t", header = TRUE
  ) |>
  inner_join(pca_plot_df |> select(ID, CLUSTER), by = "ID") |>
  mutate(
    UNIT = as.factor(UNIT), DPI15 = as.numeric(DPI15),
    TYPE = as.factor(CLUSTER),
    FLOWERING_WEEK = as.numeric(FLOWERING_WEEK)
  ) |>
  select(-CLUSTER)

# Same for the glasshouse 19 data.
data_glasshouse19 <-
  read.csv(
    "data/gp/inputs/data_glasshouse19.txt",
    sep = "\t", header = TRUE
  ) |>
  inner_join(pca_plot_df |> select(ID, CLUSTER), by = "ID") |>
  mutate(
    UNIT = as.factor(UNIT), BENCH = as.factor(BENCH),
    COLUMN = as.factor(COLUMN), DPI14 = as.numeric(DPI14),
    TYPE = as.factor(CLUSTER), FLOWERING_WEEK = as.numeric(FLOWERING_WEEK),
    STEM_WIDTH = as.numeric(STEM_WIDTH)
  ) |>
  select(-CLUSTER)

# Same for the hoophouse data.
data_hoophouse <-
  read.csv(
    "data/gp/inputs/GWAS_Poly_2016_2017_data.csv"
  ) |>
  mutate(
    LINE_NAME = sub("^\\d+ ", "", Variety, useBytes = TRUE),
    DPI14 = as.numeric(Lesion_14dpi), ENVIRONMENT = paste0("hh", Year),
    POLY = Location
  ) |>
  inner_join(varieties_key, by = "LINE_NAME") |>
  select(ID, LINE_NAME, ENVIRONMENT, DPI14, POLY) |>
  inner_join(pca_plot_df |> select(ID, CLUSTER), by = "ID") |>
  mutate(
    TYPE = as.factor(CLUSTER),
    ENVIRONMENT = as.factor(ENVIRONMENT),
    POLY = as.factor(POLY)
  ) |>
  select(-CLUSTER)

# Subset hoophouse data for the 2 different trials for convenience.
data_hh16 <-
  data_hoophouse |>
  filter(ENVIRONMENT == "hh2016") |>
  mutate(COLUMN = "dummyhh16")

data_hh17 <-
  data_hoophouse |>
  filter(ENVIRONMENT == "hh2017") |>
  mutate(COLUMN = "dummyhh17")

# Normalise the DPI15 and DPI14 data for modelling.
# Doing this improved lm fit and h2.
data_fta24$DPI15_norm <- bestNormalize(data_fta24$DPI15)$x.t
data_glasshouse19$DPI14_norm <- bestNormalize(data_glasshouse19$DPI14)$x.t
data_hh16$DPI14_norm <- bestNormalize(data_hh16$DPI14)$x.t
data_hh17$DPI14_norm <- bestNormalize(data_hh17$DPI14)$x.t

# Data filtering for self-consistency.
# Outlier reps are removed and lines without enough replication are
# removed. The different trials can have different thresholds
# but it is best to just set one universal threshold rather than
# fiddling with them to force a better correlation.
fnum1 <- 2
fnum2 <- 2
fnum3 <- 2

# Start with the FTA24 data.
# Add a filter column to data_fta24 based on difference in z score between reps.
data_fta24 <- data_fta24 |>
  mutate(ZSCORE = scale(DPI15_norm)) |>
  group_by(ID) |>
  mutate(
    MEANZ = mean(ZSCORE, na.rm = TRUE),
    RESID = abs(ZSCORE - MEANZ),
  ) |>
  ungroup() |>
  group_by(ID) |>
  mutate(
    FILTER = max(RESID >= fnum1 | is.na(ZSCORE), na.rm = TRUE)
  ) |>
  ungroup() |>
  mutate(FILTER = as.logical(FILTER)) |>
  dplyr::select(-MEANZ, -RESID, -ZSCORE)

# Add the same kind of column to the glasshouse 19 data.
# Filter aberrant infections and low replication.
data_glasshouse19 <-
  data_glasshouse19 |>
  mutate(ZSCORE = scale(DPI14_norm)) |>
  group_by(ID) |>
  mutate(
    MEANZ = mean(ZSCORE, na.rm = TRUE),
    RESID = abs(ZSCORE - MEANZ),
    N = sum(!is.na(DPI14_norm) & RESID < fnum2)
  ) |>
  ungroup() |>
  mutate(FILTER = as.logical(RESID >= fnum2 | is.na(RESID) | N < 3)) |>
  dplyr::select(-MEANZ, -RESID, -N, -ZSCORE)

# And for the two hoophouse trials.
data_hh16 <-
  data_hh16 |>
  mutate(ZSCORE = scale(DPI14_norm)) |>
  group_by(ID) |>
  mutate(
    MEANZ = mean(ZSCORE, na.rm = TRUE),
    RESID = abs(ZSCORE - MEANZ),
    N = sum(!is.na(DPI14_norm) & RESID < fnum2)
  ) |>
  ungroup() |>
  mutate(FILTER = as.logical(RESID >= fnum2 | is.na(RESID) | N < 4)) |>
  dplyr::select(-MEANZ, -RESID, -N, -ZSCORE)

# Hoophouse 17
data_hh17 <-
  data_hh17 |>
  mutate(ZSCORE = scale(DPI14_norm)) |>
  group_by(ID) |>
  mutate(
    MEANZ = mean(ZSCORE, na.rm = TRUE),
    RESID = abs(ZSCORE - MEANZ),
    N = sum(!is.na(DPI14_norm) & RESID < fnum3)
  ) |>
  ungroup() |>
  mutate(FILTER = as.logical(RESID >= fnum3 | is.na(RESID) | N < fnum3)) |>
  dplyr::select(-MEANZ, -RESID, -N, -ZSCORE)

# For convenience, subset the data based on the filter columns.
data_fta24_filter <-
  data_fta24 |>
  filter(!FILTER)

data_glasshouse19_filter <-
  data_glasshouse19 |>
  filter(!FILTER)

data_hh16_filter <-
  data_hh16 |>
  filter(!FILTER)

data_hh17_filter <-
  data_hh17 |>
  filter(!FILTER)

load_blink <- function(
  blink_prefix = "data/gwas/inputs/merge_basnpsmaf_iprune"
) {

  # Load blink.
  source("http://zzlab.net/GAPIT/gapit_functions.txt")
  source("http://zzlab.net/FarmCPU/FarmCPU_functions.txt")
  library(BLINK)

  print("Loading full blink SNP data...")
  blink_snpmat <<-
    read.big.matrix(
      paste0(blink_prefix, ".dat"),
      sep = "\t",
      head = FALSE,
      type = "char"
    )

  print("Loading full blink map...")
  blink_map <<-
    fread(
      paste0(blink_prefix, ".map"),
      head = TRUE,
      data.table = FALSE
    )

  print("Loading full blink IDs...")
  blink_ids <<-
    read.table(
    paste0(blink_prefix, ".id"),
    head = FALSE,
    col.names = c("ID")
  )

  print("Transposing full blink SNP data...")
  blink_snpmat_t <<-
    t(as.matrix(blink_snpmat))

}
