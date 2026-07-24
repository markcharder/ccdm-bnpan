source("envs/R/pedigree-env.R")

######################################################################
#
# Overview:
#
#   Here, we investigate candidate selction regions, focusing on some
#   known major selection targets. We do the following:
#
#   1. Make syteny dot-plots of region surrounding FAE on A08.
#
#       It looks like the A08 region is most closely related to the one
#       in Haya, in some derivative lines, suggesting this is where the
#       FAE allele came from, rather than Zephyr. To help confirm this,
#       synteny dot plots will be made comparing derivative lines, Haya
#       and Zephyr with the reference X5A. An alignment of an extracted
#       50 Kb region will also be plotted.
#
#
#  2. Plot Nei's H for cultivars and founders across A08 and C03.
#  3. Plot Nei's H for cultivars and founders in the candidate regions on C07.
#
######################################################################

# General variables
selection_candidates_merged <-
  read.table(
    "data/pedigree/out/selection_candidates_merged.txt",
    header = TRUE,
    sep = "\t"
  )

######################################################################

# 1. Syntey plots of the region surrounding FAE on A08.

# List the PAF files and varieties of interest.
pafs <- list.files(
  "data/pedigree/out//minimap2/a8/",
  pattern = "paf$",
  full.names = TRUE
)

varieties <- c()
for (paf in pafs) {
  var <- strsplit(basename(paf), "_")[[1]][1]
  varieties <- c(varieties, var)
}

# These are the N50s of the assemblies that were aligned
# in Mb (see PEDIGREE.md). "R" denotes a reference assembly.
varieties_n50s <-
  c(
    "R", "R", 23.17, 2.01, 1.63,
    29.10, 21.15, 21.53, 12.89, 1.42
  )

# This is the order I would like the varieties to be plotted in.
varieties_plot <- c(
  "Darmorbzh", "zs11", "Chikuzen", "Chisaya", "Haya",
  "Mystic", "Trooper", "Insignia", "Range", "Rainbow"
)

# Order the PAF list according to the desired plot order.
pafs <- pafs[match(varieties_plot, varieties)]

varieties_plot_rf <-
  gsub("Darmorbzh", "Darmor-bzh", varieties_plot)
varieties_plot_rf <-
  gsub("zs11", "ZS11", varieties_plot_rf)

# Extract a shared identity legend from the first PAF.
shared_legend <-
  get_legend(
    plot_synteny(pafs[1]) +
      theme(
        legend.position = "right",
        legend.background = element_rect(fill = "white", colour = NA),
        legend.key = element_rect(fill = "white", colour = NA)
      )
  )

svg("R/plots/A08_synteny.svg", width = 12, height = 5)
plot_grid(
  plot_synteny_panel(
    paf_files = pafs,
    titles = paste(varieties_plot_rf, varieties_n50s),
    ncol = 5,
    byrow = FALSE
  ),
  shared_legend,
  ncol = 2,
  rel_widths = c(1, 0.08)
)
dev.off()

# Read in the alignment of the 50 Kb region.
library(ggmsa)
library(Biostrings)

aln <- readDNAMultipleAlignment("snps_only.fasta")

ggmsa(aln, color = "Chemistry_NT", seq_name = TRUE) +
    theme_minimal()

######################################################################

# 2. Plot Nei's H for cultivars and founders across A08 and C03.

# Read the nei selection scan.
selscan_nei <- fread("data/pedigree/out/selscan_nei.txt")

# Read in FAE regions for C03 and A08.
eu_regions <- read.table("data/pedigree/out/Darmor_EU_regions.txt")

# Make a plot of Nei's across C03 and A08.
# Extract data tables for A08 and C03 and add a column to indicate
# where the gene lies.
selscan_nei_a8 <- selscan_nei[CHROM == "A08"]
selscan_nei_c3 <- selscan_nei[CHROM == "C03" & POS <= 40e6]

selscan_nei_plot <-
  bind_rows(selscan_nei_a8, selscan_nei_c3) |>
  select(CHROM, POS, DIV_CULTIVAR, DIV_FOUNDER) |>
  rename(
    DIV_CULTIVAR = "Cultivar",
    DIV_FOUNDER   = "Founder"
  ) |>
  pivot_longer(
    cols = c(Cultivar, Founder),
    names_to = "GROUP", values_to = "DIV"
  )

selscan_nei_vlines <-
  bind_rows(
    eu_regions |>
      mutate(
        CHROM = V1,
        POS = V2 + (V3 - V2) / 2,
        GROUP = "Cultivar",
        DIV = 0.5
      ),
    eu_regions |>
      mutate(
        CHROM = V1,
        POS = V3,
        GROUP = "Founder",
        DIV = 0.5
      )
  ) |>
  select(CHROM, POS, GROUP, DIV)

# Shade the selection windows.
selscan_nei_shade <-
  selection_candidates_merged |>
  filter(CHROM %in% c("C03", "A08")) |>
  mutate(
    GROUP = "Cultivar",
    DIV = 0.5,
    POS = WIN_START + (WIN_END - WIN_START) / 2
  ) |>
  select(CHROM, POS, WIN_START, WIN_END, GROUP, DIV)

selscan_nei_shade <-
  selscan_nei_shade |>
  bind_rows(
    selscan_nei_shade |>
      mutate(GROUP = "Founder")
  )

# Make the plot.
p_fae <-
  ggplot(
    selscan_nei_plot, aes(x = POS / 1e6, y = DIV, colour = GROUP)
  ) +
  geom_rect(
    data = selscan_nei_shade,
    aes(xmin = WIN_START / 1e6, xmax = WIN_END / 1e6, ymin = -Inf, ymax = Inf),
    fill = "darkolivegreen1", alpha = 0.8, inherit.aes = FALSE
  ) +
  facet_grid(
    GROUP ~ CHROM, scales = "free_x"
  ) +
  geom_point(show.legend = FALSE) +
  scale_color_manual(
    values = c(
      "Cultivar" = "coral3", "Founder" = "steelblue"
    )
  ) +
  labs(x = "Position (Mb)", y = "Nei's H") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12),
    text = element_text(size = 15),
    strip.background = element_blank()
  ) +
  geom_hline(
    aes(yintercept = 0.5),
    size = 1,
    linetype = "dashed",
    colour = "grey50"
  ) +
  geom_vline(
    data = selscan_nei_vlines,
    aes(xintercept = POS / 1e6),
    size = 0.8,
    colour = "goldenrod"
  )

png(
  "R/plots/nei_fae.png",
  width = 12, height = 3,
  units = "in",
  res = 300
)
p_fae
dev.off()

######################################################################

# 4. Plot the haplotypes across the A08 region.

# Haplotype plot for the erucic acid locus on A08.

# Create a temporary directory for the extracted VCF.
dir.create("R/plots/tmp", showWarnings = FALSE)

# Get the region of interest for BnaA08.FAE1.
ea_regions <-
  read.table("data/pedigree/out/Darmor_EU_regions.txt")
ea_chrom <- "A08"
ea_start <- ea_regions[1, 2]
ea_end <- ea_regions[1, 3]
ec_chrom <- "C03"
ec_start <- ea_regions[2, 2]
ec_end <- ea_regions[2, 3]

# Set the extracted VCF path.
ea_vcf_path <- "R/plots/tmp/ea_locus.vcf.gz"
ec_vcf_path <- "R/plots/tmp/ec_locus.vcf.gz"

# Extract A08 for the samples of interest.
varieties <-
  c(
    "HAYAC", "CHIKUZEN", "CHISAYA",
    "MUTU", "NORIN 40", "NORIN 41", "NORIN 8",
    "NORIN 21", "NORIN 25", "NORIN 27", "NORIN 28",
    "NORIN 2", "NORIN 4", "Nakate Chousen", "GENKAI",
    "TOKIWA", "KONGO", "GanganF73", "Shengli3",
    "Zheyou73", "ZS11", "Zephyr", "Tower",
    "Westar", "Golden", "Tanka", "Nugget", "Oro",
    "Regent", "Regina II", "Altex", "Andor", "Reston",
    "Midas", "Triton", "Bronowski", "ERGLU", "CRESOR",
    "MAJOR", "Ramses"
  )

varieties <-
  c(
    varieties,
    varieties_key_pedigree |> filter(TYPE != "founder") |>
    pull(LINE_NAME)
  )

varieties <-
  c(
    varieties_key[match(varieties, varieties_key$LINE_NAME), "ID"],
    "Liho"
  )

extract_window_vcf(
  bcf = "data/pedigree/out/mergeliho_basnpsmaf_a8.bcf.gz",
  chrom = ea_chrom,
  win_start = 1L,
  win_end = 30146539L,
  out_vcf = ea_vcf_path,
  pad_fraction = 0,
  full_chrom = TRUE,
  varieties = varieties
)

vcf <- read.vcfR(ea_vcf_path, verbose = FALSE)
gt <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
hap1 <-
  apply(
    gt, 2,
    function(x) {
      sapply(
        x, function(g) {
          if (is.na(g) || g %in% c(".", "./.", ".|.")) return(NA_integer_)
          alleles <- as.integer(strsplit(g, "[|/]")[[1]])
          sum(alleles, na.rm = TRUE)
        }
      )
    }
  )

mat <- t(hap1)
pos <- as.integer(vcf@fix[, "POS"])

vline_start <- which.min(abs(pos - ea_start))
vline_end   <- which.min(abs(pos - ea_end))
midline <- round((vline_start + vline_end) / 2)

varieties <- intersect(varieties, rownames(mat))
mat <- mat[varieties, ]

colnames(mat) <- seq_len(ncol(mat))
keyfile <- varieties_key |> filter(ID %in% varieties)
line_levels <- keyfile$LINE_NAME[match(varieties, keyfile$ID)]
plot_df <- melt(
  mat,
  varnames = c("HAPLOTYPE", "SITE"),
  value.name = "ALLELE"
) |>
  left_join(keyfile, by = c("HAPLOTYPE" = "ID")) |>
  filter(!is.na(LINE_NAME)) |>
  mutate(
    LINE_NAME = factor(
      LINE_NAME, levels = line_levels
    ),
    ALLELE = factor(ALLELE),
    POS = rep(pos, each = nrow(mat))
  )

plot_df <- plot_haplotypes(
  ea_vcf_path,
  title = sprintf(
    "Erucic acid locus  %s:%d\u2013%d", ea_chrom, ea_start, ea_end
  ),
  keyfile = varieties_key_pedigree,
  varieties = varieties,
  win_start = ea_start,
  win_end = ea_end
)

p_ea <- ggplot(
  plot_df, aes(
      x = SITE, y = LINE_NAME,
      fill = ALLELE
    )
  ) +
    geom_tile() +
    scale_fill_manual(
      values = c("0" = "goldenrod", "2" = "steelblue"),
      na.value = "grey80"
    ) +
    scale_x_continuous(
      breaks = seq(0, 54761, 2500),
      labels = unique(plot_df$POS)[
        seq(1, length(unique(plot_df$POS)), 2500)
      ] / 1000000
    ) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_classic() +
    theme(
      axis.text.y  = element_text(size = 7, margin = margin(r = 0)),
      axis.ticks.y = element_blank(),
      axis.line.y  = element_blank(),
      axis.text.x = element_text(angle = 90),
      legend.position = "none"
    ) +
  scale_fill_manual(
    values = c("0" = "azure2", "2" = "darkgreen"),
    na.value = "grey80"
  ) +
  geom_vline(aes(xintercept = midline), linetype = "dashed", colour = "red", linewidth = 1.5)

png("R/plots/ea_locus_haplotypes.png", width = 700, height = 1500, res = 150)
print(p_ea)
dev.off()

extract_window_vcf(
  bcf = "data/pedigree/out/mergeliho_basnpsmaf_c3.bcf.gz",
  chrom = ec_chrom,
  win_start = 1L,
  win_end = 20000000L,
  out_vcf = ec_vcf_path,
  pad_fraction = 0,
  full_chrom = TRUE,
  varieties = varieties
)

vcf <- read.vcfR(ec_vcf_path, verbose = FALSE)
gt <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
hap1 <-
  apply(
    gt, 2,
    function(x) {
      sapply(
        x, function(g) {
          if (is.na(g) || g %in% c(".", "./.", ".|.")) return(NA_integer_)
          alleles <- as.integer(strsplit(g, "[|/]")[[1]])
          sum(alleles, na.rm = TRUE)
        }
      )
    }
  )

mat <- t(hap1)
pos <- as.integer(vcf@fix[, "POS"])

vline_start <- which.min(abs(pos - ec_start))
vline_end   <- which.min(abs(pos - ec_end))
midline <- round((vline_start + vline_end) / 2)

varieties <- intersect(varieties, rownames(mat))
mat <- mat[varieties, ]

colnames(mat) <- seq_len(ncol(mat))
keyfile <- varieties_key |> filter(ID %in% varieties)
line_levels <- keyfile$LINE_NAME[match(varieties, keyfile$ID)]
plot_df <- melt(
  mat,
  varnames = c("HAPLOTYPE", "SITE"),
  value.name = "ALLELE"
) |>
  left_join(keyfile, by = c("HAPLOTYPE" = "ID")) |>
  filter(!is.na(LINE_NAME)) |>
  mutate(
    LINE_NAME = factor(
      LINE_NAME, levels = line_levels
    ),
    ALLELE = factor(ALLELE),
    POS = rep(pos, each = nrow(mat))
  )

p_ec <- ggplot(
  plot_df, aes(
      x = SITE, y = LINE_NAME,
      fill = ALLELE
    )
  ) +
    geom_tile() +
    scale_fill_manual(
      values = c("0" = "goldenrod", "2" = "steelblue"),
      na.value = "grey80"
    ) +
    scale_x_continuous(
      breaks = seq(1, length(unique(plot_df$POS)), 2500),
      labels = unique(plot_df$POS)[
        seq(1, length(unique(plot_df$POS)), 2500)
      ] / 1000000
    ) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_classic() +
    theme(
      axis.text.y  = element_text(size = 7, margin = margin(r = 0)),
      axis.ticks.y = element_blank(),
      axis.line.y  = element_blank(),
      axis.text.x = element_text(angle = 90),
      legend.position = "none"
    ) +
  scale_fill_manual(
    values = c("0" = "azure2", "2" = "darkgreen"),
    na.value = "grey80"
  ) +
  geom_vline(aes(xintercept = midline), colour = "red", linewidth = 1)

png("R/plots/ec_locus_haplotypes.png", width = 700, height = 1500, res = 150)
print(p_ec)
dev.off()

png("R/plots/test.png", width = 700, height = 400, res = 150)
print(p_ea)
dev.off()
eu_unique_af <- eu_unique_af |>
  inner_join(varieties_key_pedigree, by = c("FOUNDER" = "ID"))

# Then plot frequency of unique alleles across the locus.
p_eu_af <-
  ggplot(
    eu_unique_af, aes(x = POS / 1e6, y = AF, fill = LINE_NAME)
  ) +
  geom_point(shape = 21, size = 2, alpha = 0.8) +
  facet_wrap(~ CHROM, ncol = 1, scale = "free") +
  theme_classic() +
  labs(x = "Position (Mb)", y = "Frequency of unique allele") +
  theme(axis.text = element_text(size = 12), text = element_text(size = 15)) +
  scale_fill_manual(
    values = c(
      "cadetblue2", "firebrick3", "chartreuse3",
      "magenta", "goldenrod", "darkcyan"
    ),
    name = ""
  )

svg("R/plots/eu_locus_unique_af.svg", width = 15, height = 4)
p_eu_af
dev.off()

# 3. Plot Nei's H for cultivars and founders in the candidate regions on C07.

selscan_nei_c7 <- selscan_nei[CHROM == "C07" & POS >= 20e6]

selscan_nei_plot2 <-
  selscan_nei_c7 |>
  select(CHROM, POS, DIV_CULTIVAR, DIV_FOUNDER) |>
  rename(
    DIV_CULTIVAR = "Cultivar",
    DIV_FOUNDER   = "Founder"
  ) |>
  pivot_longer(
    cols = c(Cultivar, Founder),
    names_to = "GROUP", values_to = "DIV"
  )

selscan_nei_shade2 <-
  selection_candidates_merged |>
  filter(CHROM == "C07") |>
  mutate(
    GROUP = "Cultivar",
    DIV = 0.5
  ) |>
  select(CHROM, WIN_START, WIN_END, GROUP, DIV)

selscan_nei_shade2 <-
  selscan_nei_shade2 |>
  bind_rows(
    selscan_nei_shade2 |>
      mutate(GROUP = "Founder")
  )

# Make the plot.
p_c7 <-
  ggplot(
    selscan_nei_plot2, aes(x = POS / 1e6, y = DIV, colour = GROUP)
  ) +
  geom_rect(
    data = selscan_nei_shade2,
    aes(xmin = WIN_START / 1e6, xmax = WIN_END / 1e6, ymin = -Inf, ymax = Inf),
    fill = "darkolivegreen1", alpha = 0.8, inherit.aes = FALSE
  ) +
  geom_point(show.legend = FALSE) +
  facet_wrap(
    ~GROUP, scales = "free_x", ncol = 1
  ) +
  scale_color_manual(
    values = c(
      "Cultivar" = "coral3", "Founder" = "steelblue"
    )
  ) +
  labs(x = "Position (Mb)", y = "Nei's H") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12),
    text = element_text(size = 15),
    strip.background = element_blank()
  ) +
  geom_hline(
    aes(yintercept = 0.5),
    size = 1,
    linetype = "dashed",
    colour = "grey50"
  )

png(
  "R/plots/nei_c7.png",
  width = 6, height = 3,
  units = "in",
  res = 300
)
p_c7
dev.off()
