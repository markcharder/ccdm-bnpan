set.seed(123)

# Source the environment.
source("envs/R/genomodel-env.R")

chroms_rename <-
  read.table(
    "doc/info/chromosomes_rename_key.txt",
    header = TRUE
  ) |>
  mutate(chr = seq(1, 19))

gwas_blink_fta_npca <-
  fread("data/gwas/out/fta_npca_gwas.txt") |>
  left_join(chroms_rename, by = "chr") |>
  mutate(
    padj = p.adjust(P.value, method = "BH"),
    ENVIRONMENT = "FTA"
 )

gwas_blink_glasshouse_npca <-
  fread("data/gwas/out/glasshouse_npca_gwas.txt") |>
  left_join(chroms_rename, by = "chr") |>
  mutate(
    padj = p.adjust(P.value, method = "BH"),
    ENVIRONMENT = "Glasshouse"
  )

gwas_blink_ftanw_npca <-
  fread("data/gwas/out/ftanw_npca_gwas.txt") |>
  left_join(chroms_rename, by = "chr") |>
  mutate(
    padj = p.adjust(P.value, method = "BH"),
    ENVIRONMENT = "FTANW"
  )

gwas_full <-
  bind_rows(
    gwas_blink_fta_npca,
    gwas_blink_glasshouse_npca,
    gwas_blink_ftanw_npca
  )

gwas_full_sig <-
  gwas_full |>
  filter(padj < 0.05)

gwas_full_menv <-
  gwas_full |>
  filter(rs %in% gwas_full_sig$rs)