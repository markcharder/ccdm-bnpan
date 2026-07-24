source("envs/R/genomodel-env.R")

lines_interest <-
  c(
    "Sapphire", "SIREN",
    "YICKADEE", "CHISAYA", "HAYAC",
    "CHIKUZEN", "Beacon", "CLANCY")

data_fta24 <-
  data_fta24 |>
  mutate(
    COLOUR = ifelse(
      LINE_NAME %in% lines_interest,
      LINE_NAME, "Other"
    )
  )

data_fta24_ag <-
  data_fta24 |>
  group_by(LINE_NAME, COLOUR) |>
  summarise(
    SCORE = mean(
      DPI15, na.rm = TRUE
    ),
    SE = sd(DPI15) / sqrt(n()),
    .groups = "drop"
  ) |>
  mutate(
    COLOUR = factor(
      COLOUR,
      levels = c(
        lines_interest[lines_interest %in% LINE_NAME],
        "Other"
      )
    )
  ) |>
  select(LINE_NAME, SCORE, SE, COLOUR)

fta_ln_levels <-
  data_fta24_ag[order(data_fta24_ag$SCORE), ]$LINE_NAME

data_fta24_ag$LINE_NAME <-
  factor(data_fta24_ag$LINE_NAME, levels = fta_ln_levels)

data_glasshouse19_filter <-
  data_glasshouse19_filter |>
  mutate(
    COLOUR = ifelse(
        LINE_NAME %in% lines_interest,
        LINE_NAME, "Other"
    )
  )

data_glasshouse19_ag <-
  data_glasshouse19_filter |>
  group_by(
    LINE_NAME, COLOUR
  ) |>
  summarise(
    SCORE = mean(DPI14, na.rm = TRUE),
    SE = sd(DPI14) / sqrt(n()),
    .groups = "drop"
  ) |>
  mutate(
    COLOUR = factor(
      COLOUR,
      levels = c(
        lines_interest[lines_interest %in% LINE_NAME],
        "Other"
      )
    )
  ) |>
  select(
    LINE_NAME, SCORE, SE, COLOUR
  )

ln_levels <-
  data_glasshouse19_ag[order(data_glasshouse19_ag$SCORE), ]$LINE_NAME

data_glasshouse19_ag$LINE_NAME <-
  factor(
    data_glasshouse19_ag$LINE_NAME, levels = ln_levels
  )

svg("R/plots/3p1_meeting_fta.svg", width = 10, height = 4)
ggplot(data_fta24_ag, aes(x = LINE_NAME, y = SCORE, fill = COLOUR)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(
    aes(
      ymin = SCORE - SE, ymax = SCORE + SE
    ),
    width = 0.2,
    position = position_dodge(0.9)
  ) +
  theme_classic() +
  labs(x = "", y = "Lesion length at 2 weeks") +
  scale_fill_manual(
    values = c(
      "darkred", "darkblue",
      "darkgreen", "darkorange", "purple",
      "chartreuse", "darkcyan", "lightgrey"
    ),
    name = ""
  ) +
  theme(axis.text.x = element_blank())
dev.off()

svg("R/plots/3p1_meeting_gh.svg", width = 10, height = 4)
ggplot(data_glasshouse19_ag, aes(x = LINE_NAME, y = SCORE, fill = COLOUR)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(
    aes(
      ymin = SCORE - SE, ymax = SCORE + SE
    ),
    width = 0.2,
    position = position_dodge(0.9)
  ) +
  theme_classic() +
  labs(x = "", y = "Lesion length at 2 weeks") +
  scale_fill_manual(
    values = c(
      "darkred", "darkblue",
      "darkgreen", "darkorange", "magenta",
      "purple", "chartreuse", "darkcyan",
      "lightgrey"
    ),
    name = ""
  ) +
  theme(axis.text.x = element_blank())
dev.off()

data_conviron2020 <-
  read.csv(
    "data/gwas/inputs/conviron-2020.csv",
    header = TRUE,
    stringsAsFactors = FALSE
  )

data_conviron2020_treatids <-
  read.csv("data/gwas/inputs/conviron-2020_treat-ids.csv") |>
  rename(TREATMENT = treatment, CURTIN_ID = ID ) |>
  dplyr::select(CURTIN_ID, TREATMENT)

data_conviron2020_ids <-
  read.csv("data/gwas/inputs/conviron-2020_ids.csv") |>
  rbind(
    data.frame(
        CURTIN_ID = c("Bonito", "Mako", "Stingray"),
        LINE_NAME = c("Bonito", "ATR Mako", "Stingray")
    )
  )

data_conviron2020 <-
  data_conviron2020 |>
  left_join(
    data_conviron2020_treatids,
    by = "TREATMENT"
    ) |>
  left_join(data_conviron2020_ids, by = "CURTIN_ID") |>
  dplyr::select(CURTIN_ID, LINE_NAME, BENCH, ISOLATE, DPI14)

lines_compare <-
  c(
    "Bonito", "SIREN", "YICKADEE", "CHISAYA",
    "HAYAC"
  )

plot_df <-
  rbind(
    data_glasshouse19 |>
    filter(LINE_NAME %in% lines_compare) |>
    dplyr::select(LINE_NAME, DPI14) |>
    mutate(EXPERIMENT = "Glasshouse 2019"),
    data_conviron2020 |>
    filter(LINE_NAME %in% lines_compare) |>
    dplyr::select(LINE_NAME, DPI14) |>
    mutate(EXPERIMENT = "Conviron 2020")) |>
    rbind(
      data_fta24 |>
      filter(LINE_NAME %in% lines_compare) |>
      mutate(DPI14 = DPI15) |>
      dplyr::select(LINE_NAME, DPI14) |>
      mutate(EXPERIMENT = "FTA 2024")
    )

plot_df.ag <-
  with(plot_df, aggregate(DPI14 ~ LINE_NAME + EXPERIMENT, FUN = mean))
plot_df.ag$SD <- with(plot_df, aggregate(DPI14 ~ LINE_NAME + EXPERIMENT, FUN = sd))$DPI14
plot_df.ag$N <- with(plot_df, aggregate(DPI14 ~ LINE_NAME + EXPERIMENT, FUN = length))$DPI14
plot_df.ag$SE <- plot_df.ag$SD / sqrt(plot_df.ag$N) 
plot_df.ag1 <-
    with(plot_df, aggregate(DPI14 ~ LINE_NAME, FUN = mean))

plot_df.ag1 <-
  plot_df.ag1[order(plot_df.ag1$DPI14), ]

plot_df.ag$LINE_NAME <- factor(plot_df.ag$LINE_NAME, levels = plot_df.ag1$LINE_NAME)

library(ggplot2)

ggplot(plot_df.ag, aes(x = LINE_NAME, y = DPI14, fill = EXPERIMENT)) +
  facet_wrap(~ EXPERIMENT, ncol = 1) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = DPI14 - SE, ymax = DPI14 + SE), width = 0.2, position = position_dodge(0.9)) +
  theme_minimal() +
  labs(title = "DPI14 by Line and Experiment", x = "Line Name", y = "DPI14") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_df2 <- data_glasshouse19 |>
  mutate(
    LINE_FILL = ifelse(LINE_NAME %in% lines_compare, "Compare", "Other"),
  ) |>
  dplyr::select(
    LINE_NAME,
    LINE_FILL,
    DPI14
  )

plot_df2.ag <- with(plot_df2, aggregate(DPI14 ~ LINE_NAME, FUN = mean))
plot_df2.ag$SD <- with(plot_df2, aggregate(DPI14 ~ LINE_NAME, FUN = sd))$DPI14
plot_df2.ag$N <- with(plot_df2, aggregate(DPI14 ~ LINE_NAME, FUN = length))$DPI14
plot_df2.ag$SE <- plot_df2.ag$SD / sqrt(plot_df2.ag$N)

plot_df2.ag$LINE_FILL <- plot_df2.ag$LINE_NAME %in% lines_compare

plot_df2.ag1 <- with(plot_df2, aggregate(DPI14 ~ LINE_NAME, FUN = mean)) |> filter(LINE_NAME %in% shared_lines | LINE_NAME %in% lines_compare)
plot_df2.ag1 <- plot_df2.ag1[order(plot_df2.ag1$DPI14), ]
plot_df2.ag$LINE_NAME <- factor(plot_df2.ag$LINE_NAME, levels = plot_df2.ag1$LINE_NAME)

svg("cpw-talk.svg", width = 15, height = 5)
ggplot(plot_df2.ag, aes(x = LINE_NAME, y = DPI14, fill = LINE_FILL)) +
  geom_bar(stat = "identity", position = position_dodge(), show.legend=FALSE) +
  geom_errorbar(aes(ymin = DPI14 - SE, ymax = DPI14 + SE), width = 0.2, position = position_dodge(0.9)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("lightgrey", "darkred")) +
  xlab("") +
  theme(axis.text=element_text(size=15), axis.text.x=element_blank())
dev.off()

ancestry_estimates <-
  read.table("data/haplotypes/out/ancestry_estimates.txt", h = F, sep = "\t")

chromosomes_map <-
  read.table("data/haplotypes/out/all_chromosomes_map.txt", h=F, sep = "\t")

colnames(ancestry_estimates) <- c("LINE_NAME", "CHROM", "STARTIDX", "ENDIDX", "Bronowski", "Chisaya", "Haya", "Zephyr")
colnames(chromosomes_map) <- c("CHROM", "BP", "IDX")

ancestry_estimates <- ancestry_estimates |>
  left_join(chromosomes_map |> mutate(STARTIDX = IDX, START = BP), by = c("CHROM", "STARTIDX")) |>
  left_join(chromosomes_map |> mutate(ENDIDX = IDX, END = BP), by = c("CHROM", "ENDIDX")) |>
  dplyr::select(LINE_NAME, CHROM, START, END, Bronowski, Chisaya, Haya, Zephyr)

# Reshape data for plotting (wide to long format)
library(tidyr)
ancestry_long <- ancestry_estimates |>
  pivot_longer(
    cols = c(Bronowski, Chisaya, Haya, Zephyr),
    names_to = "Variety",
    values_to = "Ancestry"
  )

# Create the horizontal heat map plot with stacked colors
library(ggplot2)
p <- ggplot(ancestry_long, aes(xmin = START, xmax = END, 
                                ymin = 0, ymax = 1)) +
  geom_rect(aes(fill = Variety, alpha = Ancestry), color = NA) +
  facet_grid(CHROM ~ LINE_NAME, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c("Bronowski" = "#E41A1C", 
                                "Chisaya" = "#377EB8", 
                                "Haya" = "#4DAF4A", 
                                "Zephyr" = "#984EA3"),
                     name = "Variety") +
  scale_alpha_continuous(range = c(0, 1), name = "Ancestry") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text = element_text(size = 8),
    panel.spacing = unit(0.1, "lines") 
  ) +
  labs(x = "Position (bp)", y = "", 
       title = "Ancestry Estimates Across Chromosomes")

# Display the plot
print(p)

# Save the plot
ggsave("ancestry_heatmap.pdf", p, width = 16, height = 10)
ggsave("ancestry_heatmap.png", p, width = 16, height = 10, dpi = 300)

