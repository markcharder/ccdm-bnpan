source("envs/R/pedigree-env.R")

######################################################################
#
# Overview:
#
#   There are 3 different tests for selection that will be
#   applied. These are:
#
#   1. Allele frequency changes from founders to cultivars.
#
#       This looks for alleles that have increased in
#       frequency from the founder lines to the cultivars. It
#       aggregates scores across 100 kb sliding windows and
#       identifies 95th percentile outliers.
#
#   2. iHS (integrated haplotype score) in cultivars.
#
#       This looks for evidence of recent selection in the
#       cultivars. It aggregates scores across 100 kb sliding
#       windows and identifies 95th percentile outliers. The
#       iHS score will be high if there is one high frequency
#       haplotype that is exceptionally long. In a pedigree,
#       this could be evidence of selection for a trait from
#       a single or several related founders.
#
#   3. Nei's haplotype diversity.
#
#       This tests the diversity of haplotypes in cultivars and
#       and founders across 100 kb sliding windows. It finds
#       windows where the ratio of diversity in cultivars to
#       founders is an outlier (99th percentile). It does this
#       first whilst ignoring windows with a ratio of 0 as these
#       could be due to low diversity in the founders and drift.
#       It adds to the list of candidates windows with a founder
#       diversity of >= 0.5 and cultivar diversity of 0.
#
#   Basic statistics like pi and tajima's D are also assessed for
#   context and further interpretation. The final set of candidates
#   is computed this way:
#
#   1. Identify outliers for each of the 3 tests.
#   2. Identify windows that are:
#       a. Outliers for at least one of test 1 or 2.
#       b. Outliers for test 3.
#
# Plots:
#
#   Manhatten-style plots of the 3 tests and the pi / tajima's D scans
#   will be generated. The following plots will be made:
#
#   1. 3 panels, allele freq test, iHS, and Nei's ratio. Colours:
#     a. Alternating black / grey by chromosome.
#     b. AF outliers = red
#     c. iHS outliers = blue
#     d. Nei's outliers = green
#     e. Nei's + AF outliers = goldenrod
#     f. Nei's + iHS outliers = darksalmon
#     g. AF + iHS outliers = purple
#     h. Nei's + AF + iHS outliers = cyan
#
######################################################################

# General parameters.

# This is to order chromosomes for plotting.
chrom_order <- paste0(
  rep(c("A", "C"), c(10, 9)),
  c(paste0("0", 1:9), "10", paste0("0", 1:9))
)

######################################################################

# Test 1: allele frequency change (prefix: selscan_af).

# Define sliding window size and step size.
window_size <- 100000L
step_size <-  50000L

# Generate sliding windows across the pre-calculated allele frequency
# data.
selscan_af_win <-
  selscan_af[,
    .(WIN_START = seq(0L, max(POS), by = step_size)),
    by = CHROM
  ][, WIN_END := WIN_START + window_size - 1L]

# Set individual AF scores as point intervals.
selscan_af_point <- selscan_af[,
  .(CHROM, WIN_START = POS, WIN_END = POS, DEVIATION, P_BASE, P_MODERN)
]

# Set keys for foverlaps.
setkey(selscan_af_win, CHROM, WIN_START, WIN_END)
setkey(selscan_af_point, CHROM, WIN_START, WIN_END)

# Match each site to its window.
selscan_af_overlaps <- foverlaps(
  selscan_af_point,
  selscan_af_win,
  type = "within",
  nomatch = 0L
)

# Calculate maximum deviation per window - a window can have
# one or more variants with high deviation.
# Deviation is absolute because allele frequency can go from
# high to low or low to high depending on whether it is called
# as reference or alt.
selscan_af_windev <-
  selscan_af_overlaps[,
    .(
      MAX_DEVIATION = max(abs(DEVIATION), na.rm = TRUE),
      MEAN_P_BASE    = mean(P_BASE,    na.rm = TRUE),
      MEAN_P_MODERN  = mean(P_MODERN,  na.rm = TRUE),
      N_SNP = .N
    ),
    by = .(CHROM, WIN_START, WIN_END)
  ][order(CHROM, WIN_START)]

# Identify outliers as windows with a maximum deviation in the top
# 5% of windows. Keep Z scores to assess extremity of change.
selscan_af_windev[, Z_DEVIATION := as.numeric(scale(MAX_DEVIATION))]

selscan_af_thresh <- quantile(
  selscan_af_windev$Z_DEVIATION,
  0.95,
  na.rm = TRUE
)

selscan_af_windev[, OUTLIER := Z_DEVIATION >= selscan_af_thresh]

######################################################################

# Test 2: iHS (prefix: selscan_ihs).

# Format and filter the selscan output.
selscan_ihs_filt <-
  selscan_ihs |>
  filter(!is.na(score)) |>
  mutate(
    CHROM = factor(
      chr, levels = chrom_order
    ), WIN_START = start, WIN_END = end
  ) |>
  as.data.table()

# Get a threshold score.
selscan_ihs_thresh <-
  quantile(selscan_ihs_filt$score, 0.95)

# Identify outliers.
selscan_ihs_filt[, OUTLIER := score > selscan_ihs_thresh]

# Test 3: Nei's haplotype diversity (prefix: selscan_nei).

# Create a temporary directory for chromosome VCFs.
dir.create("R/tmp", showWarnings = FALSE)

# Set number of SNPs per window and BCF input.
pedig_ba_bcf  <- "data/pedigree/inputs/acultivar_basnpsmaf.bcf.gz"
win_snps_gw   <- 100L
step_snps_gw  <-  10L

# Get the founder IDs from the key file.
founder_ids_gw  <-
  varieties_key_pedigree |>
  filter(TYPE == "founder") |>
  pull(ID)

# Get the cultivar IDs from the key file.
cultivar_ids_gw <-
  varieties_key_pedigree |>
  filter(TYPE == "cultivar") |>
  pull(ID)

# Calculate Nei's across the genome.
selscan_nei <- rbindlist(
  lapply(
    chrom_order,
    function(chr) {

      # Temporary storage for the chromosome VCF.
      vcf_tmp <- file.path("R/tmp", paste0("hapdiv_", chr, ".vcf.gz"))

      # The bcftools command to extract the chromosome VCF.
      system(
        sprintf(
          paste0(
            "bcftools view --regions %s ",
            "--output-type z --output %s ",
            "%s && bcftools index --force %s"
          ),
          chr, vcf_tmp, pedig_ba_bcf, vcf_tmp
        )
      )
      print(paste("Processing chromosome", chr, "complete"))

      # Read the temporary VCF.
      vcf <- read.vcfR(vcf_tmp, verbose = FALSE)

      # Extract genotypes and positions from the VCF.
      gt_chr <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
      pos_chr <- as.integer(vcf@fix[, "POS"])

      # Print progress.
      print(paste("VCF for", chr, "read"))

      # If there are fewer SNPs in the chromosome than 
      # the window size, skip it.
      if (nrow(gt_chr) < win_snps_gw) return(NULL)

      # Set missing genotypes.
      is_missing <- is.na(gt_chr) | gt_chr %in% c(".", "./.", ".|.")

      # Get genotype dosages (0, 1, 2) from the GT field.
      gt_dos <- matrix(
        ifelse(is_missing, NA_integer_, nchar(gsub("[^1]", "", gt_chr))),
        nrow     = nrow(gt_chr),
        ncol     = ncol(gt_chr),
        dimnames = dimnames(gt_chr)
      )
      storage.mode(gt_dos) <- "integer"

      # Print progress.
      print("gt_dos calculated")

      # Get the genotype dose for founders and cultivars.
      gt_f <- gt_dos[, colnames(gt_dos) %in% founder_ids_gw,  drop = FALSE]
      gt_c <- gt_dos[, colnames(gt_dos) %in% cultivar_ids_gw, drop = FALSE]

      # The starts are from SNP 1 to the last SNP - window size, by step size.
      starts <- seq(1L, nrow(gt_dos) - win_snps_gw + 1L, by = step_snps_gw)

      # Calculate Nei's haplotype diversity for founders and cultivars in each
      # window using the functions defined in the env file.
      rbindlist(lapply(starts, function(s) {
        idx <- s:(s + win_snps_gw - 1L)
        data.table(
          CHROM = chr,
          POS = as.integer(median(pos_chr[idx])),
          WIN_START = pos_chr[idx[1]],
          WIN_END = pos_chr[idx[length(idx)]],
          DIV_FOUNDER  = nei_hap_div(make_haps(gt_f, idx)),
          DIV_CULTIVAR = nei_hap_div(make_haps(gt_c, idx))
        )

      }))

    }
  )
)

# Calculate ratio of cultivar to founder haplotype diversity.
selscan_nei[, RATIO := DIV_CULTIVAR / DIV_FOUNDER]
selscan_nei[, CHROM := factor(CHROM, levels = chrom_order)]

# Get outlier windows for Nei's.
selscan_nei_thresh <-
  quantile(
    selscan_nei$RATIO[selscan_nei$DIV_CULTIVAR != 0],
    0.01,
    na.rm = TRUE
  )

selscan_nei[,
  OUTLIER := !is.na(RATIO) &
    ((RATIO <= selscan_nei_thresh & DIV_CULTIVAR != 0) |
       (DIV_FOUNDER >= 0.5 & DIV_CULTIVAR == 0))
]

# Clean up temporary VCF files.
file.remove(list.files("R/tmp", full.names = TRUE))

######################################################################

# Save all intermediate results.
write.table(
  x = selscan_af_windev,
  file = "data/pedigree/out/selscan_af_windev.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  x = selscan_ihs_filt,
  file = "data/pedigree/out/selscan_ihs_filt.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  selscan_nei,
  file = "data/pedigree/out/selscan_nei.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

######################################################################

# Identify selection candidates.
# Set keys for finding foverlaps.
setkey(selscan_af_windev, CHROM, WIN_START, WIN_END)
setkey(selscan_ihs_filt, CHROM, WIN_START, WIN_END)
setkey(selscan_nei, CHROM, WIN_START, WIN_END)

# Add an ID column to selscan_nei.
selscan_nei[, NEI_ID := paste(CHROM, WIN_START, WIN_END, sep = "_")]

# Find Nei outliers that overlap with either test.
selscan_nei_afol <-
  unique(
    foverlaps(
      selscan_nei[OUTLIER == TRUE],
      selscan_af_windev[OUTLIER == TRUE],
      type = "any",
      nomatch = 0L
    )[["NEI_ID"]]
  )

selscan_nei_ihsol <-
  unique(
    foverlaps(
      selscan_nei[OUTLIER == TRUE],
      selscan_ihs_filt[OUTLIER == TRUE],
      type = "any",
      nomatch = 0L
    )[["NEI_ID"]]
  )

selscan_nei[, CATEGORY := fcase(
  !OUTLIER, "none",
  NEI_ID %chin% selscan_nei_afol &
    NEI_ID %chin% selscan_nei_ihsol, "nei + af + ihs",
  NEI_ID %chin% selscan_nei_afol, "nei + af",
  NEI_ID %chin% selscan_nei_ihsol, "nei + ihs",
  default = "nei"
)]

selscan_nei_candidates <-
  selscan_nei[!CATEGORY %in% c("nei", "none")]

# Then, merge the selection candidate windows.
selscan_nei_candmerge <-
  selscan_nei_candidates[order(CHROM, WIN_START),
    GRP := cumsum(
      WIN_START > cummax(data.table::shift(WIN_END, fill = WIN_END[1])) + 1L
    ),
    by = CHROM
  ][,
    .(
      WIN_START = min(WIN_START),
      WIN_END = max(WIN_END),
      POS = POS[which.min(RATIO)],
      RATIO_MIN = min(RATIO, na.rm = TRUE),
      DIV_FOUNDER = mean(DIV_FOUNDER, na.rm = TRUE),
      DIV_CULTIVAR = mean(DIV_CULTIVAR, na.rm = TRUE),
      CATEGORY = CATEGORY[which.min(RATIO)],
      ALL_CATEGORIES = paste(sort(unique(CATEGORY)), collapse = " + ")
    ),
    by = .(CHROM, GRP)
  ][, GRP := NULL][order(CHROM, WIN_START)]

# Write out unmerged candidate windows with their categories.
write.table(
  x = selscan_nei_candidates,
  file = "data/pedigree/out/selection_candidates.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Write out the selection candidate windows.
write.table(
  x = selscan_nei_candmerge,
  file = "data/pedigree/out/selection_candidates_merged.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

######################################################################

# Plotting.

# Generate chromosome offsets for plotting.
selscan_nei_coff <-
  selscan_nei[,
    .(CHR_LEN = max(POS)),
    by = CHROM
  ][order(CHROM)
  ][, OFFSET := c(0L, head(cumsum(CHR_LEN), -1))]

# Prepare data frames and plots.

# The selscan AF test.
selscan_af_plot <-
  merge(
    selscan_af_windev,
    selscan_nei_coff[, .(CHROM, OFFSET)],
    by = "CHROM"
  ) |>
  mutate(CHROM = factor(CHROM, levels = chrom_order)) |>
  transform(
    POS_CUM = WIN_START + OFFSET,
    ALT_COL = as.integer(CHROM) %% 2 == 0
  ) |>
  mutate(AF_ID = paste(CHROM, WIN_START, WIN_END, sep = "_")) |>
  as.data.table()

# Get windows overlapping iHS.
setkey(selscan_af_plot, CHROM, WIN_START, WIN_END)
selscan_af_ihsol <-
  foverlaps(
    selscan_af_plot,
    selscan_ihs_filt[OUTLIER == TRUE],
    type = "any",
    nomatch = 0L
  ) |>
  pull(AF_ID) |>
  unique()

# Get windows overlapping Nei's.
selscan_af_neiol <-
  foverlaps(
    selscan_af_plot,
    selscan_nei[OUTLIER == TRUE],
    type = "any",
    nomatch = 0L
  ) |>
  pull(AF_ID) |>
  unique()

# Set category column.
selscan_af_plot[, CATEGORY := fcase(
  (OUTLIER & !AF_ID %chin% selscan_af_neiol & !AF_ID %chin% selscan_af_ihsol), "AF",
  (!OUTLIER & !AF_ID %chin% selscan_af_neiol & AF_ID %chin% selscan_af_ihsol), "iHS",
  (!OUTLIER & AF_ID %chin% selscan_af_neiol & !AF_ID %chin% selscan_af_ihsol), "Nei",
  (OUTLIER & AF_ID %chin% selscan_af_neiol & !AF_ID %chin% selscan_af_ihsol), "AF + Nei",
  (!OUTLIER & AF_ID %chin% selscan_af_neiol & AF_ID %chin% selscan_af_ihsol), "iHS + Nei",
  (OUTLIER & !AF_ID %chin% selscan_af_neiol & AF_ID %chin% selscan_af_ihsol), "AF + iHS",
  (OUTLIER & AF_ID %chin% selscan_af_neiol & AF_ID %chin% selscan_af_ihsol), "AF + iHS + Nei",
  default = "none"
)]

selscan_af_mids <-
  selscan_af_plot |>
  as.data.frame() |>
  group_by(CHROM) |>
  summarise(
    MID = (min(POS_CUM) + max(POS_CUM)) / 2, .groups = "drop"
  )

p_sel_af <-
  ggplot(
    selscan_af_plot,
    aes(x = POS_CUM, y = Z_DEVIATION, colour = ALT_COL)
  ) +
  geom_hline(
    yintercept = selscan_af_thresh,
    linetype = "dashed",
    colour = "black"
  ) +
  geom_point(
    data = ~ subset(.x, CATEGORY == "none"),
    size = 1.5, alpha = 0.5
  ) +
  scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "grey60"),
                      guide = "none") +
  geom_point(
    data = ~ subset(.x, CATEGORY == "AF"),
    colour = "firebrick3",
    size = 1.5, alpha = 0.5
  ) +
  geom_point(
    data = ~ subset(.x, CATEGORY == "iHS"),
    colour = "steelblue",
    size = 1.5, alpha = 0.5
  ) +
  geom_point(
    data = ~ subset(.x, CATEGORY == "Nei"),
    colour = "chartreuse3",
    size = 1.5, alpha = 0.5
  ) +
  geom_point(
    data = ~ subset(.x, CATEGORY == "AF + Nei"),
    colour = "goldenrod2",
    size = 1.5, alpha = 0.5
  ) +
  geom_point(
    data = ~ subset(.x, CATEGORY == "iHS + Nei"),
    colour = "darksalmon",
    size = 1.5, alpha = 0.5
  ) +
  geom_point(
    data = ~ subset(.x, CATEGORY == "AF + iHS"),
    colour = "purple",
    size = 1.5, alpha = 0.5
  ) +
  geom_point(
    data = ~ subset(.x, CATEGORY == "AF + iHS + Nei"),
    colour = "cyan",
    size = 1.5, alpha = 0.5
  ) +
  scale_x_continuous(
    breaks = selscan_af_mids$MID,
    labels = selscan_af_mids$CHROM
  ) +
  labs(
    x = NULL, y = "AF Z-score"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12),
    text = element_text(size = 15)
  )

# The selscan iHS test.
selscan_ihs_plot <-
  merge(
    selscan_ihs_filt,
    selscan_nei_coff[, .(CHROM, OFFSET)],
    by = "CHROM"
  ) |>
  mutate(CHROM = factor(CHROM, levels = chrom_order)) |>
  transform(
    POS_CUM = WIN_START + OFFSET,
    ALT_COL = as.integer(CHROM) %% 2 == 0
  ) |>
  mutate(IHS_ID = paste(CHROM, WIN_START, WIN_END, sep = "_")) |>
  as.data.table()

# Get windows overlapping AF.
setkey(selscan_ihs_plot, CHROM, WIN_START, WIN_END)
setkey(selscan_af_windev, CHROM, WIN_START, WIN_END)
selscan_ihs_afol <-
  foverlaps(
    selscan_ihs_plot,
    selscan_af_windev[OUTLIER == TRUE],
    type = "any",
    nomatch = 0L
  ) |>
  pull(IHS_ID) |>
  unique()

# Get windows overlapping Nei's.
selscan_ihs_neiol <-
  foverlaps(
    selscan_ihs_plot,
    selscan_nei[OUTLIER == TRUE],
    type = "any",
    nomatch = 0L
  ) |>
  pull(IHS_ID) |>
  unique()

# Set category column.
selscan_ihs_plot[, CATEGORY := fcase(
  (OUTLIER & !IHS_ID %chin% selscan_ihs_neiol & !IHS_ID %chin% selscan_ihs_afol), "iHS",
  (!OUTLIER & !IHS_ID %chin% selscan_ihs_neiol & IHS_ID %chin% selscan_ihs_afol), "AF",
  (!OUTLIER & IHS_ID %chin% selscan_ihs_neiol & !IHS_ID %chin% selscan_ihs_afol), "Nei",
  (OUTLIER & IHS_ID %chin% selscan_ihs_neiol & !IHS_ID %chin% selscan_ihs_afol), "iHS + Nei",
  (!OUTLIER & IHS_ID %chin% selscan_ihs_neiol & IHS_ID %chin% selscan_ihs_afol), "AF + Nei",
  (OUTLIER & !IHS_ID %chin% selscan_ihs_neiol & IHS_ID %chin% selscan_ihs_afol), "AF + iHS",
  (OUTLIER & IHS_ID %chin% selscan_ihs_neiol & IHS_ID %chin% selscan_ihs_afol), "AF + iHS + Nei",
  default = "none"
)]

selscan_ihs_mids <-
  selscan_ihs_plot |>
  as.data.frame() |>
  group_by(CHROM) |>
  summarise(
    MID = (min(POS_CUM) + max(POS_CUM)) / 2, .groups = "drop"
  )

p_sel_ihs <-
  ggplot(
    selscan_ihs_plot,
    aes(x = POS_CUM, y = score, colour = ALT_COL)
  ) +
  geom_hline(
    yintercept = selscan_ihs_thresh,
    linetype = "dashed",
    colour = "black"
  ) +
  geom_point(
    data = ~ subset(.x, CATEGORY == "none"),
    size = 1.5, alpha = 0.5
  ) +
  scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "grey60"),
                      guide = "none") +
  geom_point(
    data = ~ subset(.x, CATEGORY == "AF"),
    colour = "firebrick3",
    size = 1.5, alpha = 0.5
  ) +
  geom_point(
    data = ~ subset(.x, CATEGORY == "iHS"),
    colour = "steelblue",
    size = 1.5, alpha = 0.5
  ) +
  geom_point(
    data = ~ subset(.x, CATEGORY == "Nei"),
    colour = "chartreuse3",
    size = 1.5, alpha = 0.5
  ) +
  geom_point(
    data = ~ subset(.x, CATEGORY == "AF + Nei"),
    colour = "goldenrod2",
    size = 1.5, alpha = 0.5
  ) +
  geom_point(
    data = ~ subset(.x, CATEGORY == "iHS + Nei"),
    colour = "darksalmon",
    size = 1.5, alpha = 0.5
  ) +
  geom_point(
    data = ~ subset(.x, CATEGORY == "AF + iHS"),
    colour = "purple",
    size = 1.5, alpha = 0.5
  ) +
  geom_point(
    data = ~ subset(.x, CATEGORY == "AF + iHS + Nei"),
    colour = "cyan",
    size = 1.5, alpha = 0.5
  ) +
  scale_x_continuous(
    breaks = selscan_ihs_mids$MID,
    labels = selscan_ihs_mids$CHROM
  ) +
  labs(
    x = NULL, y = "iHS"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12),
    text = element_text(size = 15)
  )

# The Nei's haplotype diversity test.
selscan_nei_plot <-
  merge(
    selscan_nei,
    selscan_nei_coff[, .(CHROM, OFFSET)],
    by = "CHROM"
  ) |>
  mutate(CHROM = factor(CHROM, levels = chrom_order)) |>
  transform(
    POS_CUM = WIN_START + OFFSET,
    ALT_COL = as.integer(CHROM) %% 2 == 0
  ) |>
  mutate(NEI_ID = paste(CHROM, WIN_START, WIN_END, sep = "_")) |>
  as.data.table()

# Get windows overlapping AF.
setkey(selscan_nei_plot, CHROM, WIN_START, WIN_END)
setkey(selscan_af_windev, CHROM, WIN_START, WIN_END)
selscan_nei_afol <-
  foverlaps(
    selscan_nei_plot,
    selscan_af_windev[OUTLIER == TRUE],
    type = "any",
    nomatch = 0L
  ) |>
  pull(NEI_ID) |>
  unique()

# Get windows overlapping Nei's.
selscan_nei_ihsol <-
  foverlaps(
    selscan_nei_plot,
    selscan_ihs_filt[OUTLIER == TRUE],
    type = "any",
    nomatch = 0L
  ) |>
  pull(NEI_ID) |>
  unique()

# Set category column.
selscan_nei_plot[, CATEGORY := fcase(
  (OUTLIER & !NEI_ID %chin% selscan_nei_ihsol & !NEI_ID %chin% selscan_nei_afol), "Nei",
  (!OUTLIER & !NEI_ID %chin% selscan_nei_ihsol & NEI_ID %chin% selscan_nei_afol), "AF",
  (!OUTLIER & NEI_ID %chin% selscan_nei_ihsol & !NEI_ID %chin% selscan_nei_afol), "iHS",
  (OUTLIER & NEI_ID %chin% selscan_nei_ihsol & !NEI_ID %chin% selscan_nei_afol), "iHS + Nei",
  (!OUTLIER & NEI_ID %chin% selscan_nei_ihsol & NEI_ID %chin% selscan_nei_afol), "AF + iHS",
  (OUTLIER & !NEI_ID %chin% selscan_nei_ihsol & NEI_ID %chin% selscan_nei_afol), "AF + Nei",
  (OUTLIER & NEI_ID %chin% selscan_nei_ihsol & NEI_ID %chin% selscan_nei_afol), "AF + iHS + Nei",
  default = "none"
)]

selscan_nei_mids <-
  selscan_nei_plot |>
  as.data.frame() |>
  group_by(CHROM) |>
  summarise(
    MID = (min(POS_CUM) + max(POS_CUM)) / 2, .groups = "drop"
  )

p_sel_nei <-
  ggplot(
    selscan_nei_plot,
    aes(x = POS_CUM, y = RATIO, colour = ALT_COL)
  ) +
  geom_hline(
    yintercept = selscan_nei_thresh,
    linetype = "dashed",
    colour = "black"
  ) +
  geom_point(
    data = ~ subset(.x, CATEGORY == "none"),
    size = 1.5, alpha = 0.5
  ) +
  scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "grey60"),
                      guide = "none") +
  geom_point(
    data = ~ subset(.x, CATEGORY == "AF"),
    colour = "firebrick3",
    size = 1.5, alpha = 0.5
  ) +
  geom_point(
    data = ~ subset(.x, CATEGORY == "iHS"),
    colour = "steelblue",
    size = 1.5, alpha = 0.5
  ) +
  geom_point(
    data = ~ subset(.x, CATEGORY == "Nei"),
    colour = "chartreuse3",
    size = 1.5, alpha = 0.5
  ) +
  geom_point(
    data = ~ subset(.x, CATEGORY == "AF + Nei"),
    colour = "goldenrod2",
    size = 1.5, alpha = 0.5
  ) +
  geom_point(
    data = ~ subset(.x, CATEGORY == "iHS + Nei"),
    colour = "darksalmon",
    size = 1.5, alpha = 0.5
  ) +
  geom_point(
    data = ~ subset(.x, CATEGORY == "AF + iHS"),
    colour = "purple",
    size = 1.5, alpha = 0.5
  ) +
  geom_point(
    data = ~ subset(.x, CATEGORY == "AF + iHS + Nei"),
    colour = "cyan",
    size = 1.5, alpha = 0.5
  ) +
  scale_x_continuous(
    breaks = selscan_nei_mids$MID,
    labels = selscan_nei_mids$CHROM
  ) +
  labs(
    x = NULL, y = "Nei's H ratio"
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 12),
    text = element_text(size = 15)
  )
# https://www.ncbi.nlm.nih.gov/sra/SRX26643859[accn]
# https://www.ncbi.nlm.nih.gov/sra/ERX10014205[accn]
# https://www.ncbi.nlm.nih.gov/sra/SRX3423551[accn]
# https://www.ncbi.nlm.nih.gov/sra/SRX25980045[accn]
# https://www.ncbi.nlm.nih.gov/sra/SRX25979878[accn]
# https://www.ncbi.nlm.nih.gov/sra/SRX1638418[accn]
# https://www.ncbi.nlm.nih.gov/sra/SRX1553476[accn]
# https://www.ncbi.nlm.nih.gov/sra/SRX26643868[accn]
# https://www.ncbi.nlm.nih.gov/sra/SRX26643809[accn]

png(
  "R/plots/selscan_combined.png",
  width = 15, height = 5,
  units = "in", res = 300
)
plot_grid(p_sel_af, p_sel_ihs, p_sel_nei, ncol = 1, align = "v")
dev.off()

# Look for inversion on Chromosome A08.
pafs_a8 <-
  list.files(
    path = "data/pedigree/out/minimap_a8",
    pattern = "\\.paf$",
    full.names = TRUE
  )

chromrename_key <-
  read.table("doc/info/chromosomes_rename_key_lens.txt", header = TRUE)
accs_a <- chromrename_key |> filter(REAL_ID %like% "A") |> pull(ACCESSION)
accs_c <- chromrename_key |> filter(REAL_ID %like% "C") |> pull(ACCESSION)

paf <-
  fread(pafs_a8[[1]], sep = "\t", header = FALSE, select = 1:12, fill = TRUE)

colnames(paf) <-
  c("qname", "qlen", "qstart", "qend", "strand",
    "tname", "tlen", "tstart", "tend", "nmatch",
    "alen", "mapq"
  )

paf_a8_region <- paf |>
  filter(tname %in% accs_a,
         tstart >= 1e6, tend <= 2e6)

paf <- paf |>
  mutate(
    IDENTITY = nmatch / alen,
    GENOME_GROUP = ifelse(tname %in% accs_a, "A", "C")
  )

pafinv <- paf |>
  filter(tname %in% accs_a) |>
  group_by(qname) |>
  filter(n_distinct(strand) > 1) |>
  arrange(qname, qstart)

ggplot(paf, aes(x = qstart / 1e6, y = IDENTITY, col = tname)) +
  geom_point() +
  facet_wrap(~GENOME_GROUP, ncol = 1) +
  labs(x = "A08 position (Mb)", y = "Alignment identity",
       colour = "X5A chromosome") +
  theme_classic()