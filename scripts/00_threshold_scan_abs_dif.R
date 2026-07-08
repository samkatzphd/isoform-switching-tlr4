#!/usr/bin/env Rscript
# Explore how switching gene/isoform counts respond to min |dIF|
# after q-value filtering. Does not change config defaults.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
})
source("utils/bootstrap.R")
cfg <- load_yaml_config("config/config.yml")
iq <- as.numeric(cfg$significance$isoform_q %||% 0.05)
gq <- as.numeric(cfg$significance$gene_q %||% 0.05)

score <- function(path, label) {
  z <- as_tibble(readRDS(path))
  z$q_i <- as.numeric(z$isoform_switch_q_value)
  z$q_g <- if ("gene_switch_q_value" %in% names(z)) as.numeric(z$gene_switch_q_value) else NA_real_
  z$dIF_n <- as.numeric(z$dIF)
  z$abs_dIF <- abs(z$dIF_n)
  z$sig_q <- is.finite(z$q_i) & z$q_i < iq &
    (is.na(z$q_g) | (is.finite(z$q_g) & z$q_g < gq))
  z$dataset <- label
  z
}

files <- list.files("data/processed", pattern = "^isoformFeatures_.*\\.rds$", full.names = TRUE)
all <- bind_rows(lapply(files, function(f) {
  lab <- sub("^isoformFeatures_(.+)\\.rds$", "\\1", basename(f))
  score(f, lab)
}))
sig <- all |> filter(.data$sig_q)

message("Significant isoforms by dataset:")
print(count(sig, .data$dataset))

message("|dIF| quantiles among q-significant isoforms:")
print(
  sig |>
    group_by(.data$dataset) |>
    summarize(
      n = n(),
      q10 = quantile(.data$abs_dIF, 0.1, na.rm = TRUE),
      q25 = quantile(.data$abs_dIF, 0.25, na.rm = TRUE),
      q50 = quantile(.data$abs_dIF, 0.5, na.rm = TRUE),
      q75 = quantile(.data$abs_dIF, 0.75, na.rm = TRUE),
      q90 = quantile(.data$abs_dIF, 0.9, na.rm = TRUE),
      .groups = "drop"
    )
)

thresholds <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3)
rows <- list()
for (th in thresholds) {
  for (lab in unique(sig$dataset)) {
    s <- sig |> filter(.data$dataset == lab, .data$abs_dIF >= th)
    rows[[length(rows) + 1L]] <- tibble(
      dataset = lab,
      min_abs_dif = th,
      n_sig_isoforms = nrow(s),
      n_sig_genes = dplyr::n_distinct(s$gene_id)
    )
  }
}
scan <- bind_rows(rows)
print(scan, n = 100)
ensure_dir("results/tables")
utils::write.csv(scan, "results/tables/min_abs_dif_threshold_scan.csv", row.names = FALSE)

fig_dir <- ensure_dir("results/figures/threshold_scan")
p1 <- ggplot(sig, aes(x = .data$abs_dIF)) +
  geom_histogram(bins = 40, fill = "#6a51a3", color = "white", linewidth = 0.1) +
  geom_vline(xintercept = c(0.1, 0.2), linetype = 2, color = "grey30") +
  facet_wrap(~ .data$dataset, scales = "free_y") +
  theme_bw(base_size = 11) +
  labs(
    title = "|dIF| among q-significant isoforms",
    subtitle = "Dashed lines mark |dIF| = 0.10 and 0.20",
    x = "|dIF|", y = "Count"
  )
ggsave(file.path(fig_dir, "fig_sig_abs_dif_hist.png"), p1, width = 9, height = 6, dpi = 180)

p2 <- ggplot(scan, aes(x = .data$min_abs_dif, y = .data$n_sig_genes, color = .data$dataset)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  theme_bw(base_size = 11) +
  labs(
    title = "Switching gene count vs min |dIF| (after q filter)",
    x = "min |dIF| threshold",
    y = "Genes with >= 1 switching isoform"
  )
ggsave(file.path(fig_dir, "fig_gene_count_vs_min_abs_dif.png"), p2, width = 8, height = 4.8, dpi = 180)

p3 <- ggplot(scan, aes(x = .data$min_abs_dif, y = .data$n_sig_isoforms, color = .data$dataset)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  theme_bw(base_size = 11) +
  labs(
    title = "Switching isoform count vs min |dIF| (after q filter)",
    x = "min |dIF| threshold",
    y = "Significant isoforms"
  )
ggsave(file.path(fig_dir, "fig_isoform_count_vs_min_abs_dif.png"), p3, width = 8, height = 4.8, dpi = 180)

# UT overlap stability
is_real <- function(x) {
  x <- as.character(x)
  !is.na(x) & nzchar(x) & !grepl("^(XLOC_|ENS[GTFP][0-9])", x)
}
pick_symbol <- function(names) {
  names <- unique(as.character(names))
  names <- names[!is.na(names) & nzchar(names)]
  if (!length(names)) return(NA_character_)
  good <- names[is_real(names)]
  if (length(good)) good[[1]] else names[[1]]
}

t <- score("data/processed/isoformFeatures_T_UT.rds", "T_UT")
u <- score("data/processed/isoformFeatures_U_UT.rds", "U_UT")
t$gene_name <- ave(as.character(t$gene_name), t$gene_id, FUN = pick_symbol)
u$gene_name <- ave(as.character(u$gene_name), u$gene_id, FUN = pick_symbol)

ov <- list()
for (th in thresholds) {
  tg <- unique(t$gene_name[t$sig_q & t$abs_dIF >= th & is_real(t$gene_name)])
  ug <- unique(u$gene_name[u$sig_q & u$abs_dIF >= th & is_real(u$gene_name)])
  shared <- length(intersect(tg, ug))
  t_only <- length(setdiff(tg, ug))
  u_only <- length(setdiff(ug, tg))
  ov[[length(ov) + 1L]] <- tibble(
    min_abs_dif = th,
    shared = shared,
    T_only = t_only,
    U_only = u_only,
    n_T = shared + t_only,
    n_U = shared + u_only,
    jaccard = if ((shared + t_only + u_only) > 0L) shared / (shared + t_only + u_only) else NA_real_
  )
}
ov <- bind_rows(ov)
print(ov)
utils::write.csv(ov, "results/tables/ut_overlap_vs_min_abs_dif.csv", row.names = FALSE)

ov_long <- bind_rows(
  ov |> transmute(min_abs_dif, class = "shared", n = .data$shared),
  ov |> transmute(min_abs_dif, class = "T_only", n = .data$T_only),
  ov |> transmute(min_abs_dif, class = "U_only", n = .data$U_only)
)
p4 <- ggplot(ov_long, aes(x = .data$min_abs_dif, y = .data$n, color = .data$class)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  theme_bw(base_size = 11) +
  labs(
    title = "UT gene overlap classes vs min |dIF|",
    x = "min |dIF|", y = "Genes"
  )
ggsave(file.path(fig_dir, "fig_ut_overlap_vs_min_abs_dif.png"), p4, width = 7.5, height = 4.5, dpi = 180)

message("Wrote threshold scan tables/figures under results/")
