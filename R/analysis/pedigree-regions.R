source("envs/R/pedigree-env.R")


# ALLELE FREQUENCY TEST

# Calculate across sliding windows.
window_size <- 100000L  # 100 kb
step_size   <-  50000L  # 50 kb step → 50 % overlap

# Build all window intervals per chromosome.
windows_dt <-
  selscan_af[,
    .(WIN_START = seq(0L, max(POS), by = step_size)),
    by = CHROM
  ][, WIN_END := WIN_START + window_size - 1L]

snps_dt <- selscan_af[,
  .(CHROM, POS_START = POS, POS_END = POS, DEVIATION, P_BASE, P_MODERN)
]

setkey(windows_dt, CHROM, WIN_START, WIN_END)
setkey(snps_dt, CHROM, POS_START, POS_END)

# Each SNP is matched to every window it falls within.
overlaps <- foverlaps(snps_dt, windows_dt, type = "within", nomatch = 0L)

# Mean DEVIATION per window.
window_dev <-
  overlaps[,
    .(
      MEAN_DEVIATION = mean(DEVIATION, na.rm = TRUE),
      MEAN_P_BASE    = mean(P_BASE,    na.rm = TRUE),
      MEAN_P_MODERN  = mean(P_MODERN,  na.rm = TRUE),
      N_SNP = .N
    ),
    by = .(CHROM, WIN_START, WIN_END)
  ][order(CHROM, WIN_START)]

# Add a Z score.
window_dev[, Z_DEVIATION := as.numeric(scale(MEAN_DEVIATION))]

# Get the 95th percentile and label candidate windows.
thresh_upper <- quantile(
  window_dev$Z_DEVIATION,
  0.95,
  na.rm = TRUE
)
window_dev[, OUTLIER := Z_DEVIATION > thresh_upper]

# Make a plot of the deviations for windows.
# Manhattan plot of Z_DEVIATION across windows.
chrom_order <- paste0(
  rep(c("A", "C"), c(10, 9)),
  c(paste0("0", 1:9), "10", paste0("0", 1:9))
)

window_dev[, CHROM := factor(CHROM, levels = chrom_order)]

# Compute cumulative offsets so chromosomes sit end-to-end on the x-axis.
chrom_offsets_af <-
  window_dev[,
    .(CHR_LEN = max(WIN_END)),
    by = CHROM
  ][order(CHROM)
  ][, OFFSET := c(0L, head(cumsum(CHR_LEN), -1))]

window_dev_plot <-
  merge(window_dev, chrom_offsets_af[, .(CHROM, OFFSET)], by = "CHROM") |>
  transform(
    POS_CUM = WIN_START + OFFSET,
    ALT_COL = as.integer(CHROM) %% 2 == 0
  )

chr_mids_af <-
  window_dev_plot |>
  as.data.frame() |>
  group_by(CHROM) |>
  summarise(MID = (min(POS_CUM) + max(POS_CUM)) / 2, .groups = "drop")

# SELSCAN iHS TEST

# Format and filter the selscan output.
selscan_out_filt <-
  selscan_out |>
  filter(!is.na(score)) |>
  mutate(
    CHROM = factor(
      chr, levels = chrom_order
    ), WIN_START = start, WIN_END = end
  ) |>
  as.data.table()

# Take the 95th percentile.
selscan_out_quant <-
  quantile(selscan_out_filt$score, 0.95)

selscan_out_filt[, OUTLIER := score > selscan_out_quant]

chrom_offsets_ihs <-
  selscan_out_filt[,
    .(CHR_LEN = max(WIN_END)),
    by = CHROM
  ][order(CHROM)
  ][, OFFSET := c(0L, head(cumsum(CHR_LEN), -1))]

selscan_plot <-
  merge(
    selscan_out_filt,
    chrom_offsets_ihs[, .(CHROM, OFFSET)],
    by = "CHROM"
  ) |>
  transform(
    POS_CUM = WIN_START + OFFSET,
    ALT_COL = as.integer(CHROM) %% 2 == 0
  )

# Chromosome midpoints for x-axis labels.
chr_mids_ihs <-
  selscan_plot |>
  as.data.frame() |>
  group_by(CHROM) |>
  summarise(
    MID = (
      min(POS_CUM) + max(POS_CUM)
    ) / 2, .groups = "drop"
  )

chr_labels_ihs <- selscan_plot |>
  group_by(CHROM) |>
  summarise(mid = (min(POS_CUM) + max(POS_CUM)) / 2, .groups = "drop")

# Now do the same for the windowed nucleotide diversity (pi) scan output.
window_pi <-
  window_pi |>
  mutate(
    CHROM = factor(CHROM, levels = chrom_order),
    WIN_START = BIN_START
  ) |>
  as.data.table()

window_pi_quant <-
  quantile(window_pi$PI, 0.05)

window_pi[, OUTLIER := PI < window_pi_quant]

chrom_offsets_pi <-
  window_pi[,
    .(CHR_LEN = max(WIN_START)),
    by = CHROM
  ][order(CHROM)
  ][, OFFSET := c(0L, head(cumsum(CHR_LEN), -1))]

window_pi_plot <-
  merge(window_pi, chrom_offsets_pi[, .(CHROM, OFFSET)], by = "CHROM") |>
  transform(
    POS_CUM = WIN_START + OFFSET,
    ALT_COL = as.integer(CHROM) %% 2 == 0
  )

chr_mids_pi <-
  window_pi_plot |>
  as.data.frame() |>
  group_by(CHROM) |>
  summarise(MID = (min(POS_CUM) + max(POS_CUM)) / 2, .groups = "drop")

# Chromosome midpoints for x-axis labels.
chr_labels_pi <- window_pi_plot |>
  group_by(CHROM) |>
  summarise(mid = (min(POS_CUM) + max(POS_CUM)) / 2, .groups = "drop")

# Find iHS windows overlapping selscan_af outlier windows.
af_outlier_wins <- as.data.table(window_dev[OUTLIER == TRUE])[,
  .(CHROM, WIN_START, WIN_END)
]

setkey(af_outlier_wins, CHROM, WIN_START, WIN_END)
setkey(selscan_out_filt, CHROM, WIN_START, WIN_END)

ihs_af_overlaps <- foverlaps(
  selscan_out_filt,
  af_outlier_wins,
  type    = "any",
  nomatch = 0L
)

ihs_af_overlaps_outliers <- ihs_af_overlaps[score > selscan_out_quant]

# Tidy up the overlaps dataframe.
ihs_af_overlaps_outliers <-
  ihs_af_overlaps_outliers[, CHROM := factor(CHROM, levels = chrom_order)][
    order(CHROM, WIN_START)
  ] |>
  mutate(
    WIN_START_AF = WIN_START,
    WIN_END_AF   = WIN_END,
    WIN_START_IHS = i.WIN_START,
    WIN_END_IHS   = i.WIN_END
  ) |>
  select(-i.WIN_START, -i.WIN_END, -WIN_START, -WIN_END) |>
  as.data.frame()

# Tag windows in each plot data frame with a three-level category.
both_af_wins <- ihs_af_overlaps_outliers |>
  select(CHROM, WIN_START = WIN_START_AF, WIN_END = WIN_END_AF) |>
  distinct() |>
  mutate(CHROM = factor(CHROM, levels = chrom_order), BOTH = TRUE)

window_dev_plot <- window_dev_plot |>
  left_join(both_af_wins, by = c("CHROM", "WIN_START", "WIN_END")) |>
  mutate(
    CATEGORY = case_when(
      BOTH    ~ "both",
      OUTLIER ~ "single",
      TRUE    ~ "normal"
    )
  ) |>
  select(-BOTH)

both_ihs_wins <- ihs_af_overlaps_outliers |>
  select(CHROM, WIN_START_IHS) |>
  mutate(WIN_START = WIN_START_IHS) |>
  distinct() |>
  mutate(CHROM = factor(CHROM, levels = chrom_order), BOTH = TRUE)

selscan_plot <- selscan_plot |>
  left_join(both_ihs_wins, by = c("CHROM", "WIN_START")) |>
  mutate(
    CATEGORY = case_when(
      BOTH ~ "both",
      OUTLIER ~ "single",
      TRUE ~ "normal"
    )
  ) |>
  select(-BOTH)

# Build the set of "both" regions from ihs_af_overlaps_outliers (AF coords).
both_wins <- ihs_af_overlaps_outliers |>
  select(CHROM, WIN_START = WIN_START_AF, WIN_END = WIN_END_AF) |>
  distinct() |>
  as.data.table()
both_wins[, CHROM := factor(CHROM, levels = chrom_order)]
setkey(both_wins, CHROM, WIN_START, WIN_END)

# Build the set of "both" regions from ihs_af_overlaps_outliers (AF coords).
both_wins <- ihs_af_overlaps_outliers |>
  select(CHROM, WIN_START = WIN_START_AF, WIN_END = WIN_END_AF) |>
  distinct() |>
  as.data.table()
both_wins[, CHROM := factor(CHROM, levels = chrom_order)]
setkey(both_wins, CHROM, WIN_START, WIN_END)

# Give pi calculations an end window coordinate.
pi_win_size <- window_size
pi_dt <- as.data.table(window_pi_plot)[,
  WIN_END := WIN_START + pi_win_size - 1L
]
setkey(pi_dt, CHROM, WIN_START, WIN_END)

# Overlap pi windows against the "both" regions.
pi_both <- foverlaps(pi_dt, both_wins, type = "any", nomatch = 0L)[,
  .(CHROM, WIN_START = i.WIN_START, BOTH = TRUE)
] |> unique()

# Join back and add CATEGORY (based on test overlap, not pi outlier).
af_wins_pi <- as.data.table(window_dev[OUTLIER == TRUE])[,
  .(CHROM, WIN_START, WIN_END)
]

ihs_wins_pi <-
  selscan_out_filt[OUTLIER == TRUE, .(CHROM, WIN_START, WIN_END)]

pi_dt2 <-
  as.data.table(window_pi_plot)[, WIN_END_PI := WIN_START + pi_win_size - 1L]

setkey(pi_dt2, CHROM, WIN_START, WIN_END_PI)
setkey(af_wins_pi,  CHROM, WIN_START, WIN_END)
setkey(ihs_wins_pi, CHROM, WIN_START, WIN_END)

pi_af_hit <- foverlaps(
  pi_dt2, af_wins_pi, by.x = c("CHROM", "WIN_START", "WIN_END_PI"),
  type = "any", nomatch = 0L
)[, .(CHROM, WIN_START = i.WIN_START, AF_HIT = TRUE)] |> unique()

pi_ihs_hit <- foverlaps(
  pi_dt2, ihs_wins_pi, by.x = c("CHROM", "WIN_START", "WIN_END_PI"),
  type = "any", nomatch = 0L
)[, .(CHROM, WIN_START = i.WIN_START, IHS_HIT = TRUE)] |> unique()

window_pi_plot <- window_pi_plot |>
  left_join(pi_af_hit,  by = c("CHROM", "WIN_START")) |>
  left_join(pi_ihs_hit, by = c("CHROM", "WIN_START")) |>
  mutate(
    AF_HIT  = !is.na(AF_HIT),
    IHS_HIT = !is.na(IHS_HIT),
    CATEGORY = case_when(
      AF_HIT & IHS_HIT ~ "both",
      AF_HIT ~ "af",
      IHS_HIT ~ "ihs",
      TRUE ~ "none"
    ),
    CATEGORY = factor(CATEGORY, levels = c("none", "af", "ihs", "both"))
  ) |>
  select(-AF_HIT, -IHS_HIT)

# Add CATEGORY to Tajima using the same logic.
window_tajima <- window_tajima |>
  mutate(
    CHROM     = factor(CHROM, levels = chrom_order),
    WIN_START = BIN_START,
    WIN_END   = BIN_START + 100000L - 1L
  ) |>
  as.data.table()

taj_dt <- copy(window_tajima)
setkey(taj_dt, CHROM, WIN_START, WIN_END)

taj_af_hit <- foverlaps(
  taj_dt, af_wins_pi, type = "any", nomatch = 0L
)[, .(CHROM, WIN_START = i.WIN_START, AF_HIT = TRUE)] |> unique()

taj_ihs_hit <- foverlaps(
  taj_dt, ihs_wins_pi, type = "any", nomatch = 0L
)[, .(CHROM, WIN_START = i.WIN_START, IHS_HIT = TRUE)] |> unique()

chrom_offsets_taj <- window_tajima[,
  .(CHR_LEN = max(WIN_START)),
  by = CHROM
][order(CHROM)
][, OFFSET := c(0L, head(cumsum(CHR_LEN), -1))]

window_tajima_plot <- merge(
  window_tajima, chrom_offsets_taj[, .(CHROM, OFFSET)], by = "CHROM"
) |>
  transform(
    POS_CUM = WIN_START + OFFSET, ALT_COL = as.integer(CHROM) %% 2 == 0
  ) |>
  left_join(taj_af_hit,  by = c("CHROM", "WIN_START")) |>
  left_join(taj_ihs_hit, by = c("CHROM", "WIN_START")) |>
  mutate(
    AF_HIT  = !is.na(AF_HIT),
    IHS_HIT = !is.na(IHS_HIT),
    CATEGORY = case_when(
      AF_HIT & IHS_HIT ~ "both",
      AF_HIT ~ "af",
      IHS_HIT ~ "ihs",
      TRUE ~ "none"
    ),
    CATEGORY = factor(CATEGORY, levels = c("none", "af", "ihs", "both"))
  ) |>
  select(-AF_HIT, -IHS_HIT)

chr_mids_taj <- window_tajima_plot |>
  as.data.frame() |>
  group_by(CHROM) |>
  summarise(MID = (min(POS_CUM) + max(POS_CUM)) / 2, .groups = "drop")

af_out  <- window_dev[OUTLIER == TRUE, .(CHROM, WIN_START, WIN_END)]
ihs_out <- selscan_out_filt[OUTLIER == TRUE, .(CHROM, WIN_START, WIN_END)]
all_out <- unique(rbindlist(list(af_out, ihs_out)))[order(CHROM, WIN_START)]

all_out[,
  GRP := cumsum(
    WIN_START > cummax(data.table::shift(WIN_END, fill = WIN_END[1]))
  ),
  by = CHROM
]

# GLOBAL DIVERSITY LOSS

# Plot global theta W decline.
theta_w_global_long <- melt(
  theta_w_global,
  id.vars = "CHROM",
  variable.name = "GROUP",
  value.name = "THETA_W"
)


# For SVs
theta_w_global_sv_long <- melt(
  theta_w_global_sv,
  id.vars = "CHROM",
  variable.name = "GROUP",
  value.name = "THETA_W"
)

theta_w_global$DPC <- (theta_w_global$PFOUNDER - theta_w_global$PCULTIVAR) / theta_w_global$PFOUNDER * 100
theta_w_global$DOC <- (theta_w_global$PFOUNDER - theta_w_global$OCULTIVAR) / theta_w_global$PFOUNDER * 100
theta_w_global$GENOME <-
  ifelse(
    theta_w_global$CHROM %like% "A", "A", "C"
  )

svg("R/plots/diversity_loss_genome.svg", width = 6, height = 6)
p_theta_w_global_genome <-
  ggplot(theta_w_global, aes(x = GENOME, y = DPC)) +
  geom_boxplot(fill = "lightblue") +
  geom_jitter(aes(fill = CHROM), shape = 21, size = 2, col = "black") +
  theme_classic() +
  labs(x = "Genome", y = "% diversity reduction") +
  theme(text = element_text(size = 15), axis.text = element_text(size = 12)) +
  scale_fill_discrete(name = "")
p_theta_w_global_genome
dev.off()

theta_w_global_sv$DPC <- (theta_w_global_sv$PFOUNDER - theta_w_global_sv$PCULTIVAR) / theta_w_global_sv$PFOUNDER * 100
theta_w_global_sv$DOC <- (theta_w_global_sv$PFOUNDER - theta_w_global_sv$OCULTIVAR) / theta_w_global_sv$PFOUNDER * 100
theta_w_global_sv$GENOME <-
  ifelse(
    theta_w_global_sv$CHROM %like% "A", "A", "C"
  )

alleles_retained$DPC <- alleles_retained[,3] / alleles_retained[,2] * 100
alleles_retained_sv$DPC <- alleles_retained_sv[,3] / alleles_retained_sv[,2] * 100

# NEI'S HAPLOTYPE DIVERSITY

# Create a temporary directory for chromosome VCFs.
dir.create("R/plots/tmp", showWarnings = FALSE)

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
hapdiv_genome <- rbindlist(lapply(chrom_order, function(chr) {

  vcf_tmp <- file.path("R/plots/tmp", paste0("hapdiv_", chr, ".vcf.gz"))

  system(sprintf(
    "bcftools view --regions %s --output-type z --output %s %s && bcftools index --force %s",
    chr, vcf_tmp, pedig_ba_bcf, vcf_tmp
  ))

  print(paste("Processing chromosome", chr, "complete"))

  vcf     <- read.vcfR(vcf_tmp, verbose = FALSE)
  gt_chr  <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
  pos_chr <- as.integer(vcf@fix[, "POS"])

  print(paste("VCF for", chr, "read"))

  if (nrow(gt_chr) < win_snps_gw) return(NULL)

  is_missing <- is.na(gt_chr) | gt_chr %in% c(".", "./.", ".|.")
  gt_dos <- matrix(
    ifelse(is_missing, NA_integer_, nchar(gsub("[^1]", "", gt_chr))),
    nrow     = nrow(gt_chr),
    ncol     = ncol(gt_chr),
    dimnames = dimnames(gt_chr)
  )
  storage.mode(gt_dos) <- "integer"

  print("gt_dos calculated")

  gt_f <- gt_dos[, colnames(gt_dos) %in% founder_ids_gw,  drop = FALSE]
  gt_c <- gt_dos[, colnames(gt_dos) %in% cultivar_ids_gw, drop = FALSE]

  starts <- seq(1L, nrow(gt_dos) - win_snps_gw + 1L, by = step_snps_gw)

  rbindlist(lapply(starts, function(s) {
    idx <- s:(s + win_snps_gw - 1L)
    data.table(
      CHROM        = chr,
      POS          = as.integer(median(pos_chr[idx])),
      WIN_START = pos_chr[idx[1]],
      WIN_END = pos_chr[idx[length(idx)]],
      DIV_FOUNDER  = nei_hap_div(make_haps(gt_f, idx)),
      DIV_CULTIVAR = nei_hap_div(make_haps(gt_c, idx))
    )
  }))
}))

# Calculate ratio of cultivar to founder haplotype diversity.
hapdiv_genome[, RATIO := DIV_CULTIVAR / DIV_FOUNDER]
hapdiv_genome[, CHROM := factor(CHROM, levels = chrom_order)]

# For the low threshold,
# only consider instances where cultivar diversity is not zero.
# Where cultivar diversity = 0, we will keep regions where
# founder diversity >= 0.5.
thresh_ratio_low <-
  quantile(
    hapdiv_genome$RATIO[hapdiv_genome$DIV_CULTIVAR != 0],
    0.01,
    na.rm = TRUE
  )

ratio_low_dt <-
  hapdiv_genome[!is.na(RATIO) & (
    (RATIO <= thresh_ratio_low & DIV_CULTIVAR != 0) |
    (DIV_FOUNDER >= 0.5 & DIV_CULTIVAR == 0)
  ), ]

setkey(ratio_low_dt, CHROM, WIN_START, WIN_END)

thresh_af_top1 <-
  quantile(window_dev$Z_DEVIATION,  0.95, na.rm = TRUE)

thresh_ihs_top1 <-
  quantile(selscan_out_filt$score,  0.95, na.rm = TRUE)

af_top1_dt <- window_dev[Z_DEVIATION >= thresh_af_top1,
  .(CHROM, WIN_START, WIN_END, Z_DEVIATION)
]

setkey(af_top1_dt, CHROM, WIN_START, WIN_END)

ihs_top1_dt <- selscan_out_filt[score >= thresh_ihs_top1,
  .(CHROM, WIN_START, WIN_END, score)
]

setkey(ihs_top1_dt, CHROM, WIN_START, WIN_END)

ratio_af_ovlp <- foverlaps(
  ratio_low_dt, af_top1_dt,
  by.x = c("CHROM", "WIN_START", "WIN_END"),
  type = "any", nomatch = 0L
)[, .(
  PEAK_AF_Z = max(Z_DEVIATION)),
  by = .(CHROM, WIN_START = i.WIN_START, WIN_END = i.WIN_END
  )]

ratio_ihs_ovlp <- foverlaps(
  ratio_low_dt, ihs_top1_dt,
  by.x    = c("CHROM", "WIN_START", "WIN_END"),
  type    = "within",
  nomatch = 0L
)[, .(
  PEAK_IHS = max(score)),
  by = .(CHROM, WIN_START = i.WIN_START, WIN_END = i.WIN_END
  )]

final_candidates <- Reduce(
  function(a, b) merge(
    a, b, by = c("CHROM", "WIN_START", "WIN_END"),
    all.x = TRUE
  ),
  list(ratio_low_dt, ratio_af_ovlp, ratio_ihs_ovlp)
)

print(final_candidates)

write.table(
  x         = final_candidates,
  file      = "data/pedigree/out/positive_selection_candidates.txt",
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

# Merge overlapping windows.
final_candidates[order(CHROM, WIN_START),
  GRP := cumsum(
    WIN_START > cummax(data.table::shift(WIN_END, fill = WIN_END[1]
  ))),
  by = CHROM
]

final_candidates_merged <- final_candidates[,
  .(
    WIN_START = min(WIN_START),
    WIN_END = max(WIN_END),
    POS = POS[which.min(RATIO)],
    RATIO_MIN = min(RATIO, na.rm = TRUE),
    DIV_FOUNDER = mean(DIV_FOUNDER, na.rm = TRUE),
    DIV_CULTIVAR = mean(DIV_CULTIVAR, na.rm = TRUE),
    PEAK_AF_Z = max(PEAK_AF_Z, na.rm = TRUE),
    PEAK_IHS = max(PEAK_IHS, na.rm = TRUE)
  ),
  by = .(CHROM, GRP)
][, GRP := NULL]

final_candidates_merged[, SUPPORT := fcase(
  is.finite(PEAK_AF_Z) & is.finite(PEAK_IHS), "both",
  is.finite(PEAK_AF_Z), "af",
  is.finite(PEAK_IHS), "ihs"
)]

final_candidates_merged[, CHROM := factor(CHROM, levels = chrom_order)]
final_candidates_merged <- final_candidates_merged[order(CHROM, WIN_START)]

# PLOTS

# Selection scan plot.
# Reproject AF and iHS data onto the hapdiv coordinate system
# so all three panels share the same x-axis.

# Cumulative positions for Manhattan-style plot.

hapdiv_genome_plot <-
  merge(hapdiv_genome, chrom_offsets_hd[, .(CHROM, OFFSET)], by = "CHROM") |>
  transform(
    POS_CUM = POS + OFFSET,
    ALT_COL = as.integer(CHROM) %% 2 == 0
  )

chr_mids_hd <-
  hapdiv_genome_plot |>
  as.data.frame() |>
  group_by(CHROM) |>
  summarise(
    MID = (min(POS_CUM) + max(POS_CUM)) / 2, .groups = "drop"
  )

chrom_offsets_hd <-
  hapdiv_genome[,
    .(CHR_LEN = max(POS)),
    by = CHROM
  ][order(CHROM)
  ][, OFFSET := c(0L, head(cumsum(CHR_LEN), -1))]

window_dev_sel <-
  merge(
    window_dev,
    chrom_offsets_hd[, .(CHROM, OFFSET)],
    by = "CHROM"
  ) |>
  transform(POS_CUM = WIN_START + OFFSET, ALT_COL = as.integer(CHROM) %% 2 == 0)

selscan_sel <-
  merge(
    selscan_out_filt, chrom_offsets_hd[, .(CHROM, OFFSET)], by = "CHROM"
  ) |>
  transform(POS_CUM = WIN_START + OFFSET, ALT_COL = as.integer(CHROM) %% 2 == 0)

ihs_af_overlaps <- foverlaps(
  selscan_out_filt,
  af_outlier_wins,
  type    = "any",
  nomatch = 0L
)

p_sel_af <- ggplot(window_dev_sel, aes(x = POS_CUM, y = Z_DEVIATION, colour = CAT)) +
  geom_point(
    data  = ~ subset(.x, startsWith(CAT, "bg_")),
    size = 0.8, alpha = 0.5
  ) +
  geom_point(
    data  = ~ subset(.x, !startsWith(CAT, "bg_")),
    size = 1.5, alpha = 0.9
  ) +
  geom_hline(yintercept = thresh_af_top1, linetype = "dashed", colour = "red") +
  scale_colour_manual(values = sel_colours, guide = "none") +
  scale_x_continuous(breaks = chr_mids_hd$MID, labels = chr_mids_hd$CHROM) +
  labs(x = NULL, y = "AF Z-score") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12), text = element_text(size = 15)
  )

p_sel_ihs <- ggplot(selscan_sel, aes(x = POS_CUM, y = score, colour = CAT)) +
  geom_point(
    data  = ~ subset(.x, startsWith(CAT, "bg_")),
    size = 0.8, alpha = 0.5
  ) +
  geom_point(
    data  = ~ subset(.x, !startsWith(CAT, "bg_")),
    size = 1.5, alpha = 0.9
  ) +
  geom_hline(yintercept = thresh_ihs_top1, linetype = "dashed", colour = "red") +
  scale_colour_manual(values = sel_colours, guide = "none") +
  scale_x_continuous(breaks = chr_mids_hd$MID, labels = chr_mids_hd$CHROM) +
  labs(x = NULL, y = "iHS score") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12), text = element_text(size = 15)
  )

p_sel_ratio <- ggplot(hapdiv_genome_plot, aes(x = POS_CUM, y = RATIO, colour = CAT)) +
  geom_point(
    data  = ~ subset(.x, startsWith(CAT, "bg_")),
    size = 0.8, alpha = 0.5
  ) +
  geom_point(
    data  = ~ subset(.x, !startsWith(CAT, "bg_")),
    size = 1.5, alpha = 0.9
  ) +
  geom_hline(yintercept = thresh_ratio_low, linetype = "dashed", colour = "red") +
  scale_colour_manual(values = sel_colours_ratio, guide = "none") +
  scale_x_continuous(breaks = chr_mids_hd$MID, labels = chr_mids_hd$CHROM) +
  labs(x = "Chromosome", y = "H ratio") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12), text = element_text(size = 15)
  )

png("R/plots/selection_combined.png", width = 4200, height = 2700, res = 300)
plot_grid(
  p_sel_af, p_sel_ihs, p_sel_ratio,
  ncol = 1, align = "v", axis = "lr",
  rel_heights = c(1, 1, 1)
)
dev.off()

### MARK FOR REMOVAL

# Find runs of low haplotype diversity.
thresh_hd_low <- quantile(
  hapdiv_genome$DIV_CULTIVAR,
  0.01,
  na.rm = TRUE
)

hapdiv_low_runs <-
  hapdiv_genome[order(CHROM, POS)][!is.na(DIV_CULTIVAR)][,
    {
      low  <- DIV_CULTIVAR <= thresh_hd_low
      r    <- rle(low)
      ends <- cumsum(r$lengths)
      strt <- c(1L, head(ends, -1L) + 1L)
      ki   <- which(r$values == TRUE)
      if (length(ki) == 0L) return(NULL)
      data.table(
        WIN_COUNT        = r$lengths[ki],
        POS_START        = POS[strt[ki]],
        POS_END          = POS[ends[ki]],
        DIV_CULTIVAR_MIN  =
          vapply(
            ki,
            function(k) min(DIV_CULTIVAR[strt[k]:ends[k]]),
            numeric(1)
          ),
        DIV_CULTIVAR_MEAN =
          vapply(
            ki,
            function(k) mean(DIV_CULTIVAR[strt[k]:ends[k]]),
            numeric(1)
          ),
        DIV_FOUNDER_MIN =
          vapply(
            ki,
            function(k) min(DIV_FOUNDER[strt[k]:ends[k]],  na.rm = TRUE),
            numeric(1)
          ),
        DIV_FOUNDER_MEAN =
          vapply(
            ki,
            function(k) mean(DIV_FOUNDER[strt[k]:ends[k]], na.rm = TRUE),
            numeric(1)
          )
      )
    },
    by = CHROM
  ][WIN_COUNT > 1L][order(CHROM, -WIN_COUNT)]

print(hapdiv_low_runs)
###

# Runs where founders are diverse (H >= 0.5) but cultivar diversity is
# at the genome-wide bottom 1st percentile — the primary selection candidates.
hapdiv_low_sel <- hapdiv_low_runs[DIV_FOUNDER_MEAN >= 0.5]
print(hapdiv_low_sel)

# Convert run coordinates to cumulative x for genome-wide shading.
hapdiv_shade <-
  merge(
    as.data.frame(hapdiv_low_sel),
    as.data.frame(chrom_offsets_hd[, .(CHROM, OFFSET)]),
    by = "CHROM"
  ) |>
  transform(
    XMIN = POS_START + OFFSET,
    XMAX = POS_END   + OFFSET
  )

# Redraw genome-wide plots with pastel-yellow shading over selection windows.
p_hd_cultivar2 <-
  ggplot(
    hapdiv_genome_plot,
    aes(
      x = POS_CUM, y = DIV_CULTIVAR, colour = ALT_COL
    )
  ) +
  geom_rect(
    data = hapdiv_shade,
    aes(
      xmin = XMIN, xmax = XMAX, ymin = -Inf, ymax = Inf
    ),
    fill = "#FFEC8B", inherit.aes = FALSE
  ) +
  geom_point(
    size = 0.8, alpha = 0.5
  ) +
  scale_colour_manual(
    values = c("FALSE" = "black", "TRUE" = "grey60"),
    guide = "none"
  ) +
  scale_x_continuous(breaks = chr_mids_hd$MID, labels = chr_mids_hd$CHROM) +
  labs(x = NULL, y = "Nei's H\n(cultiv.)") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12), text = element_text(size = 15)
  )

p_hd_ratio2 <-
  ggplot(hapdiv_genome_plot, aes(x = POS_CUM, y = RATIO, colour = ALT_COL)) +
  geom_rect(
    data        = hapdiv_shade,
    aes(xmin = XMIN, xmax = XMAX, ymin = -Inf, ymax = Inf),
    fill = "#FFEC8B", inherit.aes = FALSE
  ) +
  geom_point(size = 0.8, alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "red") +
  scale_colour_manual(
    values = c(
      "FALSE" = "steelblue", "TRUE" = "coral3"
    ),
    guide = "none"
  ) +
  scale_x_continuous(breaks = chr_mids_hd$MID, labels = chr_mids_hd$CHROM) +
  labs(x = "Chromosome", y = "H ratio") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12), text = element_text(size = 15)
  )

png("R/plots/hapdiv_genome_shaded.png", width = 4200, height = 900, res = 300)
plot_grid(p_hd_cultivar2, p_hd_ratio2, ncol = 1, align = "hv")
dev.off()

png("R/plots/A08_hapdiv.png", width = 4200, height = 900, res = 300)
plot_grid(p_a8_cultivar, p_a8_founder, ncol = 1, align = "v", axis = "lr")
dev.off()


# Genome-wide RATIO plot with candidate windows highlighted.
final_cands_plot <-
  merge(
    final_candidates,
    chrom_offsets_hd[, .(CHROM, OFFSET)],
    by = "CHROM"
  ) |>
  transform(POS_CUM = POS + OFFSET)

p_cand_ratio <-
  ggplot(hapdiv_genome_plot, aes(x = POS_CUM, y = RATIO)) +
  geom_point(aes(colour = ALT_COL), size = 0.8, alpha = 0.4) +
  geom_point(
    data  = final_cands_plot,
    aes(x = POS_CUM, y = RATIO, fill = SUPPORT),
    shape = 21, size = 2.5, colour = "black", alpha = 0.9
  ) +
  geom_hline(
    yintercept = thresh_ratio_low, linetype = "dashed", colour = "red"
  ) +
  scale_colour_manual(
    values = c("FALSE" = "grey80", "TRUE" = "grey60"),
    guide  = "none"
  ) +
  scale_fill_manual(
    values = c("af" = "lightblue", "ihs" = "red", "both" = "goldenrod"),
    name   = "Support"
  ) +
  scale_x_continuous(breaks = chr_mids_hd$MID, labels = chr_mids_hd$CHROM) +
  labs(x = "Chromosome", y = "H ratio (cultiv. / found.)") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    text        = element_text(size = 15)
  )

png(
  "R/plots/positive_selection_candidates.png",
  width = 4200, height = 900, res = 300
)
print(p_cand_ratio)
dev.off()


# Loop over all windows to inspect and save a plot for each
dir.create("R/plots/haplotypes", showWarnings = FALSE)
dir.create("R/plots/tmp", showWarnings = FALSE)

for (i in seq_len(nrow(windows_to_inspect))) {

  w <- windows_to_inspect[i]
  slug <- sprintf("%s_%s_%d_%d", w$SOURCE, w$CHROM, w$WIN_START, w$WIN_END)
  vcf_path <- file.path("R/plots/tmp", paste0(slug, ".vcf.gz"))

  win_start_ex <- max(0L, w$WIN_START - 100000L)
  max_chrom <-
    windows_to_inspect |>
    filter(CHROM == w$CHROM) |>
    summarise(MAX_END = max(WIN_END)) |>
    pull(MAX_END)

  win_end_ex <- min(max_chrom, w$WIN_END + 100000L)

  extract_window_vcf(
    bcf       = "data/pedigree/inputs/acultivar_basnpsmaf.bcf.gz",
    chrom     = as.character(w$CHROM),
    win_start = win_start_ex,
    win_end   = win_end_ex,
    out_vcf   = vcf_path
  )

  p <- plot_haplotypes(
    vcf_path,
    title = sprintf(
      "%s  %s:%d–%d", w$SOURCE, w$CHROM, w$WIN_START, w$WIN_END
    ),
    keyfile   = varieties_key_pedigree,
    win_start = w$WIN_START,
    win_end   = w$WIN_END
  )

  if (is.null(p)) {
    message(sprintf("Skipping %s — no variants in window", slug))
    next
  }

  png(
    file.path(
      "R/plots/haplotypes/",
      paste0(slug, ".png")
    ),
    width = 1000, height = 1000, res = 150
  )
  print(p)
  dev.off()
}

# Melt to long format, keeping row identity

eu_a8_haplotypes <- eu_a8_haplotypes[, c(1, 14840:15240)]
eu_a8_haplotypes[, row := V1]

eu_a8_haplotypes_long <- melt(
  eu_a8_haplotypes,
  id.vars = "row",
  variable.name = "col",
  value.name = "val"
)

labels <- c(
  "missing", "shared", "Haya", "Zephyr",
  "Westar", "Tower", "Erglu", "Oro"
)

# Encode categorical values as integers for clustering
val_levels <- c("missing", "shared", "Haya_verkko", "Zephyr_flye",
                 "Westar_flye", "V30SR", "V17SR", "V269SR")

# Wide matrix of encoded values (rows = samples, cols = variants)
wide_mat <- dcast(eu_a8_haplotypes_long, row ~ col, value.var = "val")
row_ids  <- wide_mat$row
enc_mat  <- wide_mat[, -1] |>
  lapply(\(x) as.integer(factor(x, levels = val_levels))) |>
  as.data.frame() |>
  as.matrix()
rownames(enc_mat) <- row_ids

hc       <- hclust(dist(enc_mat, method = "manhattan"), method = "ward.D2")
row_order <- hc$labels[hc$order]

eu_a8_haplotypes_long$row <- factor(eu_a8_haplotypes_long$row, levels = row_order)

eu_a8_haplotypes_long$val <-
  factor(eu_a8_haplotypes_long$val,
    levels = c(
      "missing", "shared", "Haya_verkko",
      "Zephyr_flye", "Westar_flye", "V30SR",
      "V17SR", "V269SR"
    )
  )

png("R/plots/eu_unique_haps.png", width = 2000, height = 500)
p_a8_haplotypes <- ggplot(
  eu_a8_haplotypes_long,
  aes(x = col, y = row, fill = val)
) +
geom_tile() +
scale_fill_manual(
  values = c(
    "darkgrey", "lightgrey", "steelblue",
    "coral3", "chartreuse3", "magenta",
    "goldenrod", "darkcyan"
  ),
  labels = labels
) 
p_a8_haplotypes
dev.off()

p_hd_founder <-
  ggplot(
    hapdiv_genome_plot,
    aes(x = POS_CUM, y = DIV_FOUNDER, colour = ALT_COL)
  ) +
  geom_point(size = 0.8, alpha = 0.5) +
  scale_colour_manual(
    values = c(
      "FALSE" = "black", "TRUE" = "grey60"
    ),
    guide = "none"
  ) +
  scale_x_continuous(breaks = chr_mids_hd$MID, labels = chr_mids_hd$CHROM) +
  labs(x = NULL, y = "Nei's H (founders)") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12), text = element_text(size = 15)
  )

p_hd_cultivar <-
  ggplot(
    hapdiv_genome_plot, aes(x = POS_CUM, y = DIV_CULTIVAR, colour = ALT_COL)
  ) +
  geom_point(
    size = 0.8, alpha = 0.5
  ) +
  scale_colour_manual(
    values = c(
      "FALSE" = "black", "TRUE" = "grey60"
    ),
    guide = "none"
  ) +
  scale_x_continuous(breaks = chr_mids_hd$MID, labels = chr_mids_hd$CHROM) +
  labs(x = NULL, y = "Nei's H") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12), text = element_text(size = 15)
  )

p_hd_ratio <-
  ggplot(
    hapdiv_genome_plot, aes(x = POS_CUM, y = RATIO, colour = ALT_COL)
  ) +
  geom_point(size = 0.8, alpha = 0.5) +
  geom_hline(
    yintercept = 1, linetype = "dashed", colour = "red"
  ) +
  scale_colour_manual(
    values = c("FALSE" = "steelblue", "TRUE" = "coral3"), guide = "none"
  ) +
  scale_x_continuous(breaks = chr_mids_hd$MID, labels = chr_mids_hd$CHROM) +
  labs(x = "Chromosome", y = "H ratio") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12), text = element_text(size = 15)
  )

svg("R/plots/hapdiv_genome.svg", width = 14, height = 3)
plot_grid(p_hd_cultivar, p_hd_ratio, ncol = 1, align = "hv")
dev.off()


# Redraw p1 and p2 with the three-category colour scheme.
p1 <- ggplot(window_dev_plot, aes(x = POS_CUM, y = Z_DEVIATION)) +
  geom_point(
    data = ~ subset(.x, CATEGORY == "normal"),
    aes(colour = ALT_COL),
    size = 1.5, alpha = 0.7
  ) +
  geom_point(
    data = ~ subset(.x, CATEGORY == "single"),
    colour = "lightblue",
    size = 2, alpha = 0.9
  ) +
  geom_point(
    data = ~ subset(.x, CATEGORY == "both"),
    colour = "goldenrod",
    size = 2.5, alpha = 0.9
  ) +
  geom_hline(
    yintercept = thresh_upper, linetype = "dashed", colour = "red"
  ) +
  scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "grey60"),
                      guide = "none") +
  scale_x_continuous(breaks = chr_mids_af$MID, labels = chr_mids_af$CHROM) +
  labs(x = "", y = "Z-score deviation") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    text = element_text(size = 15)
  )

p2 <- ggplot(selscan_plot, aes(x = POS_CUM, y = score)) +
  geom_point(
    data = ~ subset(.x, CATEGORY == "normal"),
    aes(colour = ALT_COL),
    size = 1.5, alpha = 0.7
  ) +
  geom_point(
    data = ~ subset(.x, CATEGORY == "single"),
    colour = "red",
    size = 2, alpha = 0.9
  ) +
  geom_point(
    data = ~ subset(.x, CATEGORY == "both"),
    colour = "goldenrod",
    size = 2.5, alpha = 0.9
  ) +
  geom_hline(
    yintercept = selscan_out_quant, linetype = "dashed", colour = "red"
  ) +
  scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "grey60"),
                      guide = "none") +
  scale_x_continuous(breaks = chr_mids_ihs$MID, labels = chr_mids_ihs$CHROM) +
  labs(x = "", y = "iHS score") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    text = element_text(size = 15)
  )

cat_colours <- c(
  "none" = "black", "af" = "lightblue", "ihs" = "red", "both" = "goldenrod"
)

p3 <- ggplot(window_pi_plot, aes(x = POS_CUM, y = PI)) +
  geom_point(
    data = ~ subset(.x, CATEGORY == "none"),
    aes(colour = ALT_COL),
    size = 1.5, alpha = 0.7
  ) +
  geom_point(
    data = ~ subset(.x, CATEGORY != "none"),
    aes(colour = CATEGORY),
    size = 2, alpha = 0.9
  ) +
  scale_colour_manual(
    values = c("FALSE" = "black", "TRUE" = "grey60",
               "af" = "lightblue", "ihs" = "red", "both" = "goldenrod"),
    guide = "none"
  ) +
  scale_x_continuous(breaks = chr_mids_pi$MID, labels = chr_mids_pi$CHROM) +
  labs(x = "", y = "Nucleotide diversity (pi)") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12),
    text = element_text(size = 15)
  )

p_tajima <- ggplot(window_tajima_plot, aes(x = POS_CUM, y = TajimaD)) +
  geom_point(
    data = ~ subset(.x, CATEGORY == "none"),
    aes(colour = ALT_COL),
    size = 1.5, alpha = 0.7
  ) +
  geom_point(
    data = ~ subset(.x, CATEGORY != "none"),
    aes(colour = CATEGORY),
    size = 2, alpha = 0.9
  ) +
  scale_colour_manual(
    values = c("FALSE" = "black", "TRUE" = "grey60",
               "af" = "lightblue", "ihs" = "red", "both" = "goldenrod"),
    guide = "none"
  ) +
  geom_hline(aes(yintercept = 0), colour = "darkblue", size = 0.8) +
  scale_x_continuous(breaks = chr_mids_taj$MID, labels = chr_mids_taj$CHROM) +
  labs(x = "", y = "Tajima's D") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    text = element_text(size = 15)
  )

p4 <- ggplot(window_pi_plot, aes(x = CATEGORY, y = PI)) +
  geom_boxplot(varwidth = TRUE, fill = "coral3") +
  labs(x = "", y = "Nucleotide diversity (pi)") +
  theme_classic() +
  scale_x_discrete(labels = c("none", "AF", "iHS", "both")) +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    text = element_text(size = 15)
  )

p5 <- ggplot(window_tajima_plot, aes(x = CATEGORY, y = TajimaD)) +
  geom_boxplot(varwidth = TRUE, fill = "coral3") +
  labs(x = "", y = "Absolute Tajima's D") +
  theme_classic() +
  scale_x_discrete(labels = c("none", "AF", "iHS", "both")) +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    text = element_text(size = 15)
  )

png("R/plots/selscan_combined.png", width = 2000, height = 4000, res = 300)
plot_grid(p1, p2, p3, p_tajima, ncol = 1, align = "v", axis = "lr",
          rel_heights = c(1, 1, 1, 1))
dev.off()

svg("R/plots/pi_tajima_box.svg", width = 12, height = 6)
plot_grid(p4, p5, ncol = 2, align = "hv", axis = "lr")
dev.off()


svg("R/plots/diversity_loss.svg", width = 12, height = 6)
p_theta_w_global <-
  ggplot(theta_w_global_long, aes(x = CHROM, y = THETA_W, fill = GROUP)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Chromosome", y = expression(theta[W])) +
  scale_fill_manual(
    values = c(
      "PFOUNDER" = "steelblue",
      "PCULTIVAR" = "coral3",
      "OCULTIVAR" = "magenta"
    ),
    labels = c(
      "PFOUNDER" = "Founders",
      "PCULTIVAR" = "Cultivars (pedigree)",
      "OCULTIVAR" = "Cultivars (other)"
    ),
    name   = NULL
  ) +
  theme_classic() +
  theme(text = element_text(size = 15), axis.text = element_text(size = 12))
p_theta_w_global

p_alleles_retained <-
  ggplot(
    alleles_retained, aes(x = CHROM, y = SRETAINED / SFOUNDER * 100)
  ) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Chromosome", y = "% founder alleles") +
  theme_classic() +
  theme(text = element_text(size = 15), axis.text = element_text(size = 12))
plot_grid(
  p_theta_w_global, p_alleles_retained, ncol = 1, align = "v", axis = "lr"
  )
dev.off()


svg("R/plots/diversity_loss_sv.svg", width = 12, height = 6)
p_theta_w_global_sv <-
  ggplot(theta_w_global_sv_long, aes(x = CHROM, y = THETA_W, fill = GROUP)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Chromosome", y = expression(theta[W])) +
  scale_fill_manual(
    values = c(
      "PFOUNDER" = "steelblue",
      "PCULTIVAR" = "coral3",
      "OCULTIVAR" = "magenta"
    ),
    labels = c(
      "PFOUNDER" = "Founders",
      "PCULTIVAR" = "Cultivars (pedigree)",
      "OCULTIVAR" = "Cultivars (other)"
    ),
    name   = NULL
  ) +
  theme_classic() +
  theme(text = element_text(size = 15), axis.text = element_text(size = 12))

p_alleles_retained_sv <-
  ggplot(
    alleles_retained_sv, aes(x = CHROM, y = SRETAINED / SFOUNDER * 100)
  ) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Chromosome", y = "% founder alleles") +
  theme_classic() +
  theme(text = element_text(size = 15), axis.text = element_text(size = 12))
plot_grid(
  p_theta_w_global_sv, p_alleles_retained_sv, ncol = 1, align = "v", axis = "lr"
  )
dev.off()
