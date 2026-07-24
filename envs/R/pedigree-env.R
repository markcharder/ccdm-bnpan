require(dplyr)
require(tidyr)
require(ggplot2)
require(ggrepel)
require(pheatmap)
require(cowplot)
require(data.table)
require(GenomicRanges)
require(vcfR)
require(reshape2)

# Nei's haplotype diversity.

nei_hap_div <- function(haps) {
  haps <- haps[!is.na(haps)]
  n    <- length(haps)
  if (n < 2L) return(NA_real_)
  freqs <- table(haps) / n
  (n / (n - 1L)) * (1 - sum(freqs^2))
}

make_haps <- function(mat, idx) {
  sub <- mat[idx, , drop = FALSE]
  hap <- apply(sub, 2, paste, collapse = "")
  hap[colSums(is.na(sub)) > 0L] <- NA_character_
  hap
}

# For dotplots.
generate_colors <- function(n){

  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("Install RColorBrewer.")

  }

  if (n <= 12) {
    RColorBrewer::brewer.pal(n, "Set3")[1:n]

  } else {
    colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n)

  }

}

# This is to read paf files and make a plot.
plot_synteny <-
  function(
    input_file,
    identity_lower = 0,
    min_query_aln = 400000,
    min_align = 10000,
    min_ref_len = 1000000,
    plot_size = 15
  ) {

  # Read alignments from PAF.
  alignments = read.table(
    input_file,
    stringsAsFactors = FALSE,
    row.names = NULL,
    fill = TRUE, header = FALSE
  )[, c(1:12)]

  # Set to double.
  alignments[, c(2:4, 7:12)] <- apply(
    alignments[, c(2:4, 7:12)], 2, as.double
  )

  # Column names.
  colnames(
    alignments
  )[1:12] <- c(
    "queryID", "queryLen", "queryStart",
    "queryEnd", "strand", "refID",
    "refLen", "refStart", "refEnd",
    "numResidueMatches", "lenAln", "mapQ"
  )

  # Get average percent identity across an alignment.
  alignments$percentID <- alignments$numResidueMatches / alignments$lenAln

  if (identity_lower) {
    alignments$percentID[
      which(
        alignments$percentID < identity_lower
      )
    ] <- identity_lower
  }

  queryStartTemp <- alignments$queryStart

  # Flip starts, ends for negative strand alignments
  alignments$queryStart[
    which(alignments$strand == "-")
  ] <- alignments$queryEnd[
    which(alignments$strand == "-")
  ]
  alignments$queryEnd[
    which(alignments$strand == "-")
  ] <- queryStartTemp[
    which(alignments$strand == "-")
  ]
  rm(queryStartTemp)

  cat(paste0("\nNumber of alignments: ", nrow(alignments), "\n"))
  cat(paste0("Number of query sequences: ", length(unique(alignments$queryID)), "\n"))

  # Sort by reference chromosome sizes.
  chromMax <- tapply(alignments$refLen, alignments$refID, max)
  refIDsToKeepOrdered <- names(sort(chromMax, decreasing = TRUE))

  # Filter queries and alignments by length.
  queryLenAgg <- tapply(alignments$lenAln, alignments$queryID, sum)
  alignments <- alignments[
    which(alignments$queryID %in% names(queryLenAgg)[which(queryLenAgg > min_query_aln)]),
  ]

  alignments <- alignments[which(alignments$lenAln > min_align),]
  alignments <- alignments[which(alignments$refLen > min_ref_len),]
  queryLenAgg <- tapply(alignments$lenAln, alignments$queryID, sum)
  alignments <- alignments[
    which(alignments$queryID %in% names(queryLenAgg)[which(queryLenAgg > min_query_aln)]),
  ]

  cat(paste0("\nAfter filtering... Number of alignments: ", nrow(alignments),"\n"))
  cat(paste0("After filtering... Number of query sequences: ", length(unique(alignments$queryID)),"\n\n"))

  if (nrow(alignments) == 0) {
    cat("Error: All alignments were discarded after filtering, no valid alignments remain.\n\n.")
    quit()
  }

  # Filter ref IDs.
  refIDsToKeepOrdered <- refIDsToKeepOrdered[
    which(refIDsToKeepOrdered %in% alignments$refID)
  ]

  # Sort dataframe on ref.
  alignments$refID <- factor(
    alignments$refID, levels = refIDsToKeepOrdered
  )
  alignments <- alignments[
    with(alignments,order(refID,refStart)),
  ]
  chromMax <- tapply(
    alignments$refLen, alignments$refID, max
  )

  # Make new dot plot coords.
  alignments$refStart2 <-
    alignments$refStart + sapply(
      as.character(alignments$refID),
      function(x) ifelse(
          x == names(chromMax)[1], 0, cumsum(as.numeric(chromMax))[match(x, names(chromMax)) - 1]
      )
    )

  alignments$refEnd2 <-
    alignments$refEnd + sapply(
      as.character(alignments$refID), function(x) ifelse(
        x == names(chromMax)[1], 0, cumsum(as.numeric(chromMax))[match(x, names(chromMax)) - 1]
      )
    )

  # Sort query ID based on longest alignment.
  alignments$queryID <- factor(
    alignments$queryID,
    levels=unique(as.character(alignments$queryID))
  )
  queryMaxAlnIndex <- tapply(
    alignments$lenAln, alignments$queryID, which.max, simplify = FALSE
  )
  alignments$queryID <- factor(
    alignments$queryID, levels = unique(
      as.character(alignments$queryID)
    )[order(mapply(
      function(x, i)
      alignments$refStart2[which(i == alignments$queryID)][x],
      queryMaxAlnIndex,
       names(queryMaxAlnIndex)
    ))]
  )

  queryLenAggPerRef <-
    sapply(
      (
        levels(alignments$queryID)), function(x) tapply(
          alignments$lenAln[which(alignments$queryID == x)], alignments$refID[which(alignments$queryID == x)],
          sum
      )
    )
  if(length(levels(alignments$refID)) > 1){
    queryID_Ref <- apply(queryLenAggPerRef, 2, function(x) rownames(queryLenAggPerRef)[which.max(x)])
  
  } else {
    queryID_Ref <- sapply(queryLenAggPerRef, function(x) names(queryLenAggPerRef)[which.max(x)])
  
  }
  
  alignments$queryID <- factor(alignments$queryID, levels = (levels(alignments$queryID))[order(match(queryID_Ref, levels(alignments$refID)))])
  queryMax <- tapply(alignments$queryLen, alignments$queryID, max)

  # Flip query coords if necesary.
  queryRevComp <-
    tapply(
      alignments$queryEnd - alignments$queryStart, alignments$queryID, function(x) sum(x)
    ) < 0
  queryRevComp <-
    names(queryRevComp)[which(queryRevComp)]
  alignments$queryStart[which(alignments$queryID %in% queryRevComp)] <-
    queryMax[match(as.character(alignments$queryID[which(alignments$queryID %in% queryRevComp)]), names(queryMax))] -
    alignments$queryStart[which(alignments$queryID %in% queryRevComp)] + 1
  alignments$queryEnd[which(alignments$queryID %in% queryRevComp)] <-
    queryMax[match(as.character(alignments$queryID[which(alignments$queryID %in% queryRevComp)]), names(queryMax))] -
    alignments$queryEnd[which(alignments$queryID %in% queryRevComp)] + 1

  alignments$queryStart2 <-
    alignments$queryStart +
    sapply(as.character(alignments$queryID), function(x) ifelse(x == names(queryMax)[1], 0, cumsum(queryMax)[match(x, names(queryMax)) - 1]) )
  alignments$queryEnd2 <-
    alignments$queryEnd +
    sapply(as.character(alignments$queryID), function(x) ifelse(x == names(queryMax)[1], 0, cumsum(queryMax)[match(x, names(queryMax)) - 1]) )

  options(warn = -1) # turn off warnings
  gp <- ggplot(alignments) +
    theme_bw() +
    theme(
      text = element_text(size = 12),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.y = element_text(angle = 15),
      axis.text.x = element_text(hjust = 1, angle = 45)
    ) +
    scale_color_distiller(
      palette = "Spectral",
      limits  = c(0.90, 1.00),
      oob     = scales::squish,
      breaks  = c(0.90, 0.95, 0.97, 0.99, 1.00),
      labels  = function(x) paste0(x * 100, "%")
    ) +
    labs(color = "Identity")

  if (length(unique(alignments$refID)) == 1){
    reflen <- unique(alignments$refLen)
    xbreaks <- seq(0, reflen, reflen/10)

    if (reflen/10 > 1e9){
      xlables <- paste(round(xbreaks/1e9), "GB")

    } else if (reflen/10 > 1e6) {
      xlables <- paste(round(xbreaks/1e6), "MB")

    } else if(reflen/10 > 1e3){
      xlables <- paste(round(xbreaks/1e3), "KB")

    } else{
      xlables <- paste(round(xbreaks), "bp")
    }

    gp <- gp + scale_x_continuous(
      expand = c(0, 0), limits = c(0, reflen + 0.1), breaks = xbreaks, labels = xlables
    ) +
      xlab(unique(alignments$refID))

  } else {
    gp <- gp + 
      theme(panel.grid.major.x=element_blank()) +
      geom_vline(xintercept = cumsum(as.numeric(chromMax)), col="#ebebeb") + 
      scale_x_continuous(expand = c(0, 0), limits = c(0, sum(as.numeric(chromMax)) + 0.1), 
        breaks = cumsum(as.numeric(chromMax)) - chromMax/2, labels = substr(levels(alignments$refID), start = 1, stop = 20)) + 
      xlab("Reference")
  }

  if (length(unique(alignments$queryID)) == 1){
    queryLen <- unique(alignments$queryLen)
    ybreaks <- seq(0, queryLen, queryLen/10)

    if (queryLen/10 > 1e9){
      ylables <- paste(round(ybreaks/1e9), "GB")

    } else if (queryLen/10 > 1e6) {
      ylables <- paste(round(ybreaks/1e6), "MB")

    } else if(queryLen/10 > 1e3){
      ylables <- paste(round(ybreaks/1e3), "KB")

    } else{
      ylables <- paste(round(ybreaks), "bp")

    }

    gp <- gp + scale_y_continuous(expand = c(0, 0), limits = c(0, queryLen + 0.1), breaks = ybreaks, labels = ylables) +
      ylab(unique(alignments$queryID))

  } else{
    gp <- gp +
      theme(
        panel.grid.major.y=element_blank()
      ) +
      geom_hline(
        yintercept = cumsum(as.numeric(queryMax)), col="#ebebeb"
      ) +
      scale_y_continuous(
        expand = c(0, 0), limits = c(
          0, sum(as.numeric(queryMax)) + 0.1
        ),
        breaks = cumsum(
          as.numeric(queryMax)
        ) -
        queryMax / 2,
        labels = substr(levels(alignments$queryID), start = 1, stop = 20)
      ) +
      ylab("Query")
  }

  gp <-
    gp +
    geom_segment(
      aes(
        x = refStart2,
        xend = refEnd2,
        y = queryStart2,
        yend = queryEnd2,
        color = percentID
      )
    )

  return(gp)

}

plot_synteny_panel <- function(
  paf_files, titles,
  ncol = NULL, byrow = TRUE, ...
  ) {
  stopifnot(length(paf_files) == length(titles))
  plots <- mapply(
    function(f, ttl) {
      p <- plot_synteny(f, ...)
      p +
        ggtitle(ttl) +
        theme(
          axis.text.x     = element_blank(),
          axis.ticks.x    = element_blank(),
          axis.text.y     = element_blank(),
          axis.ticks.y    = element_blank(),
          axis.title.x    = element_blank(),
          axis.title.y    = element_blank(),
          legend.position = "none"
        )
    },
    paf_files, titles,
    SIMPLIFY = FALSE
  )
  if (is.null(ncol)) ncol <- ceiling(sqrt(length(plots)))
  plot_grid(plotlist = plots, ncol = ncol, byrow = byrow)
}

# Functions for plotting haplotypes.
extract_window_vcf <-
  function(
    bcf, chrom, win_start,
    win_end, out_vcf, varieties,
    pad_fraction = 0.2,
    full_chrom = FALSE
  ) {

    # Get the size in base pairs of the padding to extract.
    pad <- as.integer((win_end - win_start) * pad_fraction)

    # If we want just the region.
    if (!full_chrom) {

      # Extract the region on the chromosome, start / end + pad.
      region <-
        sprintf(
          "%s:%d-%d", chrom, max(1L, win_start - pad),
          win_end + pad
        )

      cmd <- sprintf(
        paste0(
          "bcftools view --samples %s --regions %s ",
          "--force-samples --output-type z --output-file ",
          "%s %s && bcftools index %s"
        ),
        paste(varieties, collapse = ","), region, out_vcf, bcf, out_vcf
      )

    } else {

      # If we want the whole chromosome.
      cmd <- sprintf(
        paste0(
          "bcftools view --samples %s --regions %s ",
          "--force-samples --output-type z --output-file ",
          "%s %s && bcftools index %s"
        ),
        paste(varieties, collapse = ","), chrom, out_vcf, bcf, out_vcf
      )

    }

    # Run the command and return the VCF.
    system(cmd)
    invisible(out_vcf)

  }

plot_haplotypes <- function(
  vcf_path, title = "",
  keyfile, varieties,
  win_start, win_end
) {

  # Read the VCF.
  print("Reading the VCF...")
  vcf <- read.vcfR(vcf_path, verbose = FALSE)
  print("VCF read.")

  # Extract phased GT matrix (sites x samples)
  print("Extracting GT matrix...")
  gt <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
  print("GT matrix extracted.")

  # Sums the alleles for each genotype.
  print("Summing alleles across haplotypes...")
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
  print("Alleles summed.")

  # Transpose to samples x sites.
  print("Transposing haplotype matrix...")
  mat <- t(hap1)
  print("Haplotype matrix transposed.")

  # Map window boundaries in base pairs to nearest site.
  print("Mapping window boundaries to nearest sites...")
  pos <- as.integer(vcf@fix[, "POS"])
  print("Window boundaries mapped.")

  if (length(pos) == 0L) {
    warning(
      sprintf(
        "No variants in VCF %s — returning NULL", vcf_path
      )
    )
    return(NULL)
  }

  print("Preparing to plot haplotypes...")
  # These are the positions to draw the vertical lines -
  # start and end of the gene.
  vline_start <- which.min(abs(pos - win_start))
  vline_end   <- which.min(abs(pos - win_end))
  midline <- round((vline_start + vline_end) / 2)

  # The user specifies the order of the varieties, and 
  # which to keep, in the varieties vector.
  varieties <- intersect(varieties, rownames(mat))
  mat <- mat[varieties, ]

  # Use a plain numeric index as the x axis.
  colnames(mat) <- seq_len(ncol(mat))

  # Get line levels from the key file.
  keyfile <- keyfile |> filter(ID %in% varieties)
  line_levels <- keyfile$LINE_NAME[match(varieties, keyfile$ID)]

  # Long format for ggplot.
  df <- melt(
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
  print("Data prepared for plotting.")

  print("Returning plot...")
  return(df)

}

# Get the modified VanRaden GRM.
pedig_ldfilt_plink <- "data/pedigree/inputs/pedig_ldfilt_plink"
snpmat_plink <-
  fread(paste0(pedig_ldfilt_plink, ".raw"), data.table = FALSE)

snpmat_extract <- as.matrix(snpmat_plink[, 7:ncol(snpmat_plink)])
rownames(snpmat_extract) <- snpmat_plink$IID

# Flip: plink2 counts REF allele, we need ALT allele count.
snpmat_alt <- 2 - snpmat_extract

# Read founder base allele frequencies.
variant_ids_file <- "data/pedigree/inputs/variant_ids.txt"
variant_ids <-
  read.table(
    variant_ids_file, col.names = c("CHROM", "POS", "ID")
  )

pedig_base_af_file <- "data/pedigree/inputs/pedig_base_af.txt"
pedig_base_af <-
  read.table(
    pedig_base_af_file, col.names = c("CHROM", "POS", "AF")
  ) |>
  left_join(variant_ids, by = c("CHROM", "POS"))

# plink2 names columns as CHROM:POS_REF_ALT_counted-allele.
snp_keys <- sub("_.*", "", colnames(snpmat_alt))
colnames(snpmat_alt) <- snp_keys

# Align AF vector to SNP matrix column order.
af_matches <-
  pedig_base_af$AF[
    match(snp_keys, pedig_base_af$ID)
  ]

# Drop any SNPs with no AF match.
keep <- !is.na(af_matches)
snpmat_alt <- snpmat_alt[, keep]
pvect <- af_matches[keep]

# Canadian lines:
# STELLAR DH
# CRESOR
# TARGET
# Tribute
# Tower
# Q2
# Altex
# Andor
# Golden
# Midas
# Regent
# Triton
# Reston

# Compute VanRaden GRM manually.
zmat <- sweep(snpmat_alt, 2, 2 * pvect, "-")
denom <- 2 * sum(pvect * (1 - pvect))
kinmat_founder <- tcrossprod(zmat) / denom

# The key to get variety names.
varieties_key_pedigree <-
  read.table(
    "doc/info/varieties_key_pedigree.txt",
    sep = "\t",
    header = TRUE
  )

# Variety IDs - maps the bcf / kinship IDs to the line names in the AGG.
varieties_key <-
  read.csv(
    "doc/info/varieties_key.txt",
    header = TRUE, sep = "\t"
  )

varieties_key <-
  varieties_key |> bind_rows(
    data.frame(
      ID = c("Liho", "ganganF73_ref", "no2127_ref", "quintaA_ref", "shengli3_ref", "zheyou73_ref", "zs11_ref"),
      LINE_NAME = c("Liho", "GanganF73", "No2127", "QuintaA", "Shengli3", "Zheyou73", "ZS11")
    )
  )

# Calculate the sum of the off-diagonal elements of the GRM.
pedigcultiv_ids <-
  varieties_key_pedigree |>
  filter(TYPE == "cultivar") |>
  pull(ID)

kinmat_cultiv <-
  kinmat_founder[pedigcultiv_ids, pedigcultiv_ids]

total_sum <- sum(kinmat_cultiv)
diag_sum <- sum(diag(kinmat_cultiv))
sumg_offdiag <- total_sum - diag_sum
print(sumg_offdiag)

# Subset for just the lines in the pedigree.
pedigree_key <-
  varieties_key_pedigree |>
  filter(ID %in% rownames(kinmat_founder)) |>
  mutate(LINE_NAME = ifelse(LINE_NAME == "HAYAC", "HAYA", LINE_NAME))

# Get the outputs from selscan.
selscan_out_files <-
  list.files(
    "data/pedigree/inputs/",
    pattern = "ihs_cultiv_.*\\.windows",
    full.names = TRUE
  )

for (i in seq_along(selscan_out_files)) {

  selscan_out_tmp <-
    read.table(
      selscan_out_files[[i]],
      header = TRUE
    )

  if (i == 1) {
    selscan_out <- selscan_out_tmp

  } else {
    selscan_out <- rbind(selscan_out, selscan_out_tmp)

  }

}

selscan_ihs <- selscan_out

selscan_af <-
  fread(
    "data/pedigree/out/cultivpedig_ba_selscan_af.txt",
  )

window_pi <-
  fread(
    "data/pedigree/out/cultivall_basnpsmaf_bcf_pi.txt.windowed.pi"
  )

window_tajima <-
  fread(
    "data/pedigree/inputs/cultivall_ba_tajima.Tajima.D"
  )

eu_unique_af <-
  read.table(
    "data/pedigree/out/cultivars_founders_eu_af.txt",
    header = FALSE
  )

colnames(eu_unique_af) <-c(
  "CHROM", "POS", "FOUNDER", "N_CULTIVARS", "AF", "ID_LIST"
  )

eu_a8_haplotypes <-
  fread(
    "data/pedigree/out/cultivars_founders_eu_compare-a8.txt",
    head = FALSE,
    fill = TRUE
  )

eu_c3_haplotypes <-
  fread(
    "data/pedigree/out/cultivars_founders_eu_compare-c3.txt",
    head = FALSE,
    fill = TRUE
  )

theta_w_global <-
  read.table(
    "data/pedigree/out/theta_w_pedigree.txt",
    header = TRUE
  )

alleles_retained <-
  read.table(
    "data/pedigree/out/alleles_retained.txt",
    header = TRUE
  )

theta_w_global_sv <-
  read.table(
    "data/pedigree/out/theta_w_pedigree_sv.txt",
    header = TRUE
  )

alleles_retained_sv <-
  read.table(
    "data/pedigree/out/alleles_retained_sv.txt",
    header = TRUE
  )