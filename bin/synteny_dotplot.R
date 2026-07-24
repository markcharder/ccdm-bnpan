#!/usr/bin/env Rscript

library(ggplot2)

# Supply the minimap2 paf and GAP start+ends.
args <- commandArgs(trailingOnly = TRUE)
paf_file <- args[1]
variety <- args[2]

### Process the PAF file ###

# Read the paf.
paf <- 
    read.csv(
        paf_file,
        sep = "\t",
        header = FALSE,
        )

# Drop unused columns and rename.
paf <- paf[, 1:12]
colnames(paf) <- 
    c(
        "qname",
        "qlen",
        "qstart",
        "qend",
        "strand",
        "tname",
        "tlen",
        "tstart",
        "tend",
        "matches",
        "alnlen",
        "mapq"
    )

# Filter paf to remove matches below 10 Kb.
min_length <- 10000
paf_filtered <-
    paf[paf$alnlen >= min_length, ]

# Clean up sequence names - strip PanSN prefix (sample#hap#name) if present,
# otherwise pass through unchanged.
strip_pansn <- function(x) {
    ifelse(grepl("#", x), sub("^[^#]+#[^#]+#", "", x), x)
}
paf_filtered$qname <- strip_pansn(paf_filtered$qname)
paf_filtered$tname <- strip_pansn(paf_filtered$tname)

# Sort sequence names naturally if they match [ABC][0-9]+ (chromosome IDs),
# otherwise fall back to plain alphabetical sort.
natural_sort <- function(x) {
    if (all(grepl('^[A-C][0-9]+$', x))) {
        suffixes <- as.numeric(sub('[A-C]', '', x))
        prefixes <- sub('[0-9]+', '', x)
        x[order(prefixes, suffixes)]
    } else {
        sort(x)
    }
}
target_order <- natural_sort(unique(paf_filtered$tname))

# For the sorted reference contigs,
# get maximum lengths. All lengths for the
# same contig should be identical.
target_lengths <-
    sapply(
        target_order,
        function(x) max(
            paf_filtered$tlen[paf_filtered$tname == x]
        )
    )

# Cumulative sum of target lengths.
# Starts at zero, adds c1 length, then c2 length,
# up to the last contig, but does not add the last contig.
# The 0 index will match the first contig in target_lengths.
target_cumsum <- 
    c(0, cumsum(target_lengths)[-length(target_lengths)])
names(target_cumsum) <- target_order

# Order query contigs by their median alignment start on the reference,
# so they appear in positional order on the y-axis.
query_order <- names(
    sort(
        tapply(
            paf_filtered$tstart,
            paf_filtered$qname,
            median
        )
    )
)
query_lengths <- 
    sapply(
        query_order,
        function(x) max(paf_filtered$qlen[paf_filtered$qname == x])
    )
query_cumsum <- 
    c(0, cumsum(query_lengths)[-length(query_lengths)])
names(query_cumsum) <- query_order

# Get the alignment boundaries per query contig for label placement.
query_label_pos <- sapply(query_order, function(x) {
    rows <- paf_filtered[paf_filtered$qname == x, ]
    query_cumsum[x] + max(rows$qend)
})

# Add the alignment starts and ends to the cumulative coordinates.
paf_filtered$target_cum_start <- 
    target_cumsum[paf_filtered$tname] + paf_filtered$tstart
paf_filtered$target_cum_end <- 
    target_cumsum[paf_filtered$tname] + paf_filtered$tend
paf_filtered$query_cum_start <-
    query_cumsum[paf_filtered$qname] + paf_filtered$qstart
paf_filtered$query_cum_end <- 
    query_cumsum[paf_filtered$qname] + paf_filtered$qend

### Make the plot ###

output_file <- paste0(variety, "-synteny.svg")
svg(output_file, width = 10, height = 10)

# Build label data frame for query contig names.
query_labels <- data.frame(
    x     = max(paf_filtered$target_cum_end),
    y     = query_label_pos,
    label = query_order
)

ggplot(paf_filtered) +
    geom_segment(
        aes(
            x = target_cum_start,
            xend = target_cum_end,
            y = query_cum_start,
            yend = query_cum_end,
            color = strand
            ),
        alpha = 0.6,
        linewidth = 0.3,
        show.legend = FALSE
        ) +
    geom_text(
        data = query_labels,
        aes(x = x, y = y, label = label),
        hjust = 0,
        size = 2.5,
        nudge_x = max(paf_filtered$target_cum_end) * 0.005
    ) +
    scale_color_manual(
        values = c(
            "+" = "#2E86AB", "-" = "#A23B72"
        )
    ) +
    scale_x_continuous(
        breaks = unlist(lapply(target_order, function(chr) {
            start <- target_cumsum[chr]
            len   <- target_lengths[which(target_order == chr)]
            mb_breaks <- pretty(c(0, len), n = 4)
            mb_breaks <- mb_breaks[mb_breaks >= 0 & mb_breaks <= len]
            start + mb_breaks
        })),
        labels = unlist(lapply(target_order, function(chr) {
            len <- target_lengths[which(target_order == chr)]
            mb_breaks <- pretty(c(0, len), n = 4)
            mb_breaks <- mb_breaks[mb_breaks >= 0 & mb_breaks <= len]
            paste0(round(mb_breaks / 1e6, 1), " Mb")
        })),
        sec.axis = sec_axis(
            ~ .,
            breaks = sapply(target_order, function(chr) {
                target_cumsum[chr] + target_lengths[which(target_order == chr)] / 2
            }),
            labels = target_order
        )
    ) +
    scale_y_continuous(
        breaks = unlist(lapply(query_order, function(ctg) {
            start <- query_cumsum[ctg]
            len   <- query_lengths[which(query_order == ctg)]
            mb_breaks <- pretty(c(0, len), n = 4)
            mb_breaks <- mb_breaks[mb_breaks >= 0 & mb_breaks <= len]
            start + mb_breaks
        })),
        labels = unlist(lapply(query_order, function(ctg) {
            len <- query_lengths[which(query_order == ctg)]
            mb_breaks <- pretty(c(0, len), n = 4)
            mb_breaks <- mb_breaks[mb_breaks >= 0 & mb_breaks <= len]
            paste0(round(mb_breaks / 1e6, 1), " Mb")
        }))
    ) +
    theme_bw() +
    theme(
        axis.text.x       = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y       = element_text(size = 6, colour = "grey50"),
        axis.text.x.top   = element_text(angle = 45, hjust = 0, size = 6, colour = "grey50"),
        plot.margin       = margin(5, 80, 5, 5)
    ) +
    coord_fixed(ratio = 1, clip = "off") +
    labs(x = "Reference", y = "Query (Mb)") +
    ggtitle(variety)

dev.off()
