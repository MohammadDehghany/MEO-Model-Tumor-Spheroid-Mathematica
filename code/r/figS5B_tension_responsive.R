#!/usr/bin/env Rscript

# ------------------------------------------------------------
# Figure S5B reproduction script (published diffusion Smart-seq3D spheroids)
# Tension-responsive gene program (quadratic NB GLM vs null)
#
# NOTE
# - This script reanalyzes a published dataset. The input data files
#   (e.g., allCounts.rds, FACS_bc.csv) are not intended to be redistributed here.
#
# Method summary
# - Score: fractional abundance = sum(UMIs in set) / total UMIs per cell
# - Xi correlation (XICOR)
# - NB GLM (quadratic vs null): y_set ~ r + r^2 + offset(log(total_UMIs))
# - Plot: background cell-level points + binned mean + error bars + NB fit
#
# Outputs (default):
#   results_geneset/Fig_Tension_responsive_quad.pdf
#   results_geneset/Tension_responsive_quad_summary.csv
#   results_geneset/Tension_responsive_quad_genes_used.csv
#   results_geneset/Tension_responsive_quad_sessionInfo.txt
#
# How to run (from repo root):
#   Rscript transcriptomics/figS5B_tension_responsive_quad.R
#   Rscript transcriptomics/figS5B_tension_responsive_quad.R --counts=allCounts.rds --facs=FACS_bc.csv --experiment=MDA_SS_3D
# ------------------------------------------------------------

# -------------------- SIMPLE ARG PARSER ----------------------
parse_kv_args <- function(args) {
  # accepts --key=value
  out <- list()
  for (a in args) {
    if (!grepl("^--[^=]+=.+", a)) next
    key <- sub("^--([^=]+)=.+$", "\\1", a)
    val <- sub("^--[^=]+=(.+)$", "\\1", a)
    out[[key]] <- val
  }
  out
}

args <- commandArgs(trailingOnly = TRUE)
kv <- parse_kv_args(args)

# -------------------- USER SETTINGS (DEFAULTS) ---------------
COUNTS_RDS <- if (!is.null(kv$counts)) kv$counts else "allCounts.rds"
FACS_CSV   <- if (!is.null(kv$facs)) kv$facs else "FACS_bc.csv"
EXPERIMENT <- if (!is.null(kv$experiment)) kv$experiment else "MDA_SS_3D"

OUT_DIR <- if (!is.null(kv$out_dir)) kv$out_dir else "results_geneset"
OUT_TAG <- if (!is.null(kv$out_tag)) kv$out_tag else "Tension_responsive_quad"

N_BINS     <- if (!is.null(kv$bins)) as.integer(kv$bins) else 10
MIN_BIN_N  <- if (!is.null(kv$min_bin_n)) as.integer(kv$min_bin_n) else 30
VIS_SUBSAMPLE_N <- if (!is.null(kv$vis_n)) as.integer(kv$vis_n) else 20000
SEED <- if (!is.null(kv$seed)) as.integer(kv$seed) else 1

FIG_PDF <- file.path(OUT_DIR, paste0("Fig_", OUT_TAG, ".pdf"))
SUMMARY_CSV <- file.path(OUT_DIR, paste0(OUT_TAG, "_summary.csv"))
GENES_USED_CSV <- file.path(OUT_DIR, paste0(OUT_TAG, "_genes_used.csv"))
SESSIONINFO_TXT <- file.path(OUT_DIR, paste0(OUT_TAG, "_sessionInfo.txt"))

# --- Background dots (visualization only) ---
BG_POINT_COLOR  <- "grey35"
BG_POINT_ALPHA  <- 0.22
BG_POINT_SIZE   <- 0.55

# Binned points
BIN_POINT_SIZE   <- 4.2
BIN_POINT_FILL   <- "#E31A1C"
BIN_POINT_EDGE   <- "white"
BIN_POINT_STROKE <- 1.1

# Error bars
ERRBAR_COLOR     <- "black"
ERRBAR_WIDTH     <- 0.018
ERRBAR_LINEWIDTH <- 0.9

# Fit line
FIT_LINE_WIDTH <- 1.3

# Tensile/tension-responsive gene set (as provided)
my_geneset <- c(
  "AKR1B1","AKR1C1",
  "AIM1","CYR61","DLC1","DUSP1","ECT2","GADD45B","GLS","SCHIP1","SERPINE1","SHCBP1","TNS1",
  "CD274","CDH4","CLDN1","NEGR1","PVR"
)

# -------------------- PRE-FLIGHT CHECKS -----------------------
if (!file.exists(COUNTS_RDS)) stop("Counts file not found: ", COUNTS_RDS, call. = FALSE)
if (!file.exists(FACS_CSV))   stop("FACS file not found: ", FACS_CSV, call. = FALSE)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -------------------- PACKAGE CHECK --------------------------
required_pkgs <- c("dplyr","tibble","readr","ggplot2","MASS","XICOR","scales")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing R packages: ", paste(missing_pkgs, collapse = ", "),
    "\nInstall with: install.packages(c(", paste(sprintf('"%s"', missing_pkgs), collapse = ", "), "))",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(MASS)
  library(XICOR)
  library(scales)
})

writeLines(capture.output(sessionInfo()), SESSIONINFO_TXT)

# -------------------- THEME ----------------------------------
theme_nature <- function(base_size = 11, base_family = "sans") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(size = base_size + 2, face = "bold", hjust = 0),
      plot.subtitle   = ggplot2::element_text(size = base_size, hjust = 0, margin = ggplot2::margin(b = 6)),
      axis.title      = ggplot2::element_text(size = base_size + 1),
      axis.text       = ggplot2::element_text(size = base_size),
      axis.line       = ggplot2::element_line(linewidth = 0.9, colour = "black"),
      axis.ticks      = ggplot2::element_line(linewidth = 0.8, colour = "black"),
      axis.ticks.length = grid::unit(2.2, "mm"),
      panel.grid      = ggplot2::element_blank(),
      plot.margin     = ggplot2::margin(8, 10, 8, 8, "pt")
    )
}

# -------------------- HELPERS --------------------------------
load_counts_gene_by_cell <- function(counts_rds_path) {
  allCounts <- readr::read_rds(counts_rds_path)
  if (!("GeneID" %in% names(allCounts))) stop("Expected GeneID column in allCounts.rds", call. = FALSE)
  
  geneNames <- as.character(allCounts$GeneID)
  geneNames <- make.unique(geneNames)  # guard duplicates
  
  mat <- allCounts %>%
    dplyr::select(-GeneID) %>%
    as.matrix()
  
  storage.mode(mat) <- "numeric"
  mat <- t(mat)
  colnames(mat) <- geneNames
  mat
}

load_facs_with_radius <- function(facs_csv_path, experiment) {
  facs_raw <- readr::read_csv(facs_csv_path, show_col_types = FALSE)
  
  facs <- facs_raw %>%
    dplyr::filter(Experiment == experiment) %>%
    dplyr::mutate(
      PE_val   = readr::parse_number(as.character(`PE.Cy5..YG..A`)),
      FITC_val = readr::parse_number(as.character(`FITC.A`))
    ) %>%
    dplyr::filter(is.finite(PE_val), is.finite(FITC_val), PE_val > 0, FITC_val > 0) %>%
    dplyr::mutate(
      calcein = log10(FITC_val)
    )
  
  rng <- range(facs$calcein, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || (rng[2] - rng[1]) <= 0) {
    stop("Cannot compute radialPosition: Calcein range is zero or non-finite.", call. = FALSE)
  }
  
  facs <- facs %>%
    dplyr::mutate(
      radialPosition = (calcein - rng[1]) / (rng[2] - rng[1])
    ) %>%
    dplyr::select(barcode, radialPosition)
  
  facs$radialPosition <- pmin(1, pmax(0, facs$radialPosition))
  facs
}

fit_nb_quadratic <- function(df_nb) {
  m0 <- MASS::glm.nb(y_set ~ offset(log(N_total)), data = df_nb)
  m2 <- MASS::glm.nb(y_set ~ r + I(r^2) + offset(log(N_total)), data = df_nb)
  
  ll0 <- as.numeric(stats::logLik(m0))
  ll2 <- as.numeric(stats::logLik(m2))
  LR  <- 2 * (ll2 - ll0)
  df  <- attr(stats::logLik(m2), "df") - attr(stats::logLik(m0), "df")
  p_lrt <- stats::pchisq(LR, df = df, lower.tail = FALSE)
  
  co <- summary(m2)$coefficients
  list(
    p_lrt  = p_lrt,
    beta0  = unname(co["(Intercept)", "Estimate"]),
    beta_r = unname(co["r", "Estimate"]),
    beta_r2 = unname(co["I(r^2)", "Estimate"])
  )
}

classify_trend <- function(beta_r, beta_r2) {
  b <- beta_r; c <- beta_r2
  r_star <- if (is.finite(b) && is.finite(c) && c != 0) -b / (2 * c) else NA_real_
  in_domain <- is.finite(r_star) && r_star > 0 && r_star < 1
  
  trend_class <- dplyr::case_when(
    in_domain && c < 0 ~ "intermediate peak",
    in_domain && c > 0 ~ "U-shape (minimum)",
    is.finite(b) && is.finite(c) && b > 0 && (b + 2*c) > 0 ~ "periphery-increasing (monotone)",
    is.finite(b) && is.finite(c) && b < 0 && (b + 2*c) < 0 ~ "core-increasing (monotone decrease)",
    TRUE ~ "other/weak"
  )
  
  peak_location <- dplyr::case_when(
    in_domain && r_star < 1/3 ~ "peak closer to core",
    in_domain && r_star > 2/3 ~ "peak closer to periphery",
    in_domain ~ "peak intermediate",
    TRUE ~ NA_character_
  )
  
  list(r_star = r_star, trend_class = trend_class, peak_location = peak_location)
}

binned_summary <- function(scores_df, bins = 10) {
  brks <- seq(0, 1, length.out = bins + 1)
  scores_df %>%
    dplyr::mutate(r_bin = cut(r, breaks = brks, include.lowest = TRUE)) %>%
    dplyr::group_by(r_bin) %>%
    dplyr::summarise(
      r_center  = mean(r),
      mean_frac = mean(frac),
      se        = stats::sd(frac) / sqrt(dplyr::n()),
      ci_low    = pmax(mean_frac - 1.96 * se, 0),
      ci_high   = mean_frac + 1.96 * se,
      n         = dplyr::n(),
      .groups   = "drop"
    ) %>%
    dplyr::arrange(r_center)
}

# -------------------- MAIN -----------------------------------
counts_mat <- load_counts_gene_by_cell(COUNTS_RDS)
facs <- load_facs_with_radius(FACS_CSV, EXPERIMENT)

counts_df <- as.data.frame(counts_mat) %>%
  tibble::rownames_to_column("barcode")

genes_in_matrix <- setdiff(colnames(counts_df), "barcode")
counts_df$N_total <- rowSums(as.matrix(counts_df[, genes_in_matrix, drop = FALSE]))

cf <- facs %>%
  dplyr::inner_join(counts_df, by = "barcode") %>%
  dplyr::filter(is.finite(radialPosition), is.finite(N_total), N_total > 0)

genes_requested <- unique(my_geneset)
genes_present <- intersect(genes_requested, colnames(cf))
genes_missing <- setdiff(genes_requested, genes_present)

if (length(genes_present) == 0) stop("None of the genes in the set are present in the matrix.", call. = FALSE)
if (length(genes_missing) > 0) message("Missing genes: ", paste(genes_missing, collapse = ", "))

readr::write_csv(
  tibble::tibble(
    gene_requested = genes_requested,
    present_in_matrix = genes_requested %in% colnames(cf)
  ),
  GENES_USED_CSV
)

y_set <- rowSums(as.matrix(cf[, genes_present, drop = FALSE]))
frac  <- y_set / cf$N_total

scores_df <- tibble::tibble(
  barcode = cf$barcode,
  r       = cf$radialPosition,
  y_set   = as.numeric(y_set),
  N_total = as.numeric(cf$N_total),
  frac    = as.numeric(frac)
) %>%
  dplyr::filter(is.finite(r), is.finite(frac), is.finite(y_set), is.finite(N_total), N_total > 0)

f_min_plot <- 1 / sum(scores_df$N_total)
scores_df  <- scores_df %>% dplyr::mutate(frac_plot = frac + f_min_plot)

set.seed(SEED)
vis_df <- if (is.finite(VIS_SUBSAMPLE_N) && nrow(scores_df) > VIS_SUBSAMPLE_N) {
  scores_df %>% dplyr::slice_sample(n = VIS_SUBSAMPLE_N)
} else {
  scores_df
}

# XICOR is on CRAN, install.packages("XICOR") works.
xi_p <- XICOR::xicor(scores_df$frac, scores_df$r, pvalue = TRUE)$pval  # returns list incl pval

nb_fit <- fit_nb_quadratic(scores_df %>% dplyr::select(y_set, N_total, r))
trend  <- classify_trend(nb_fit$beta_r, nb_fit$beta_r2)

fc_peri_core     <- exp(nb_fit$beta_r + nb_fit$beta_r2)
log2fc_peri_core <- (nb_fit$beta_r + nb_fit$beta_r2) / log(2)

bin_df <- binned_summary(scores_df, bins = N_BINS) %>%
  dplyr::mutate(
    mean_plot    = mean_frac + f_min_plot,
    ci_low_plot  = ifelse(n >= MIN_BIN_N, ci_low + f_min_plot, NA_real_),
    ci_high_plot = ifelse(n >= MIN_BIN_N, ci_high + f_min_plot, NA_real_)
  )

w <- bin_df$mean_frac / sum(bin_df$mean_frac)
w[!is.finite(w)] <- 0
avgR <- sum(bin_df$r_center * w)
entropy_H <- -sum(ifelse(w > 0, w * log(w), 0))

rr <- seq(0, 1, length.out = 300)
fit_curve <- tibble::tibble(
  r = rr,
  frac_fit      = exp(nb_fit$beta0 + nb_fit$beta_r * r + nb_fit$beta_r2 * r^2),
  frac_fit_plot = frac_fit + f_min_plot
)

y_low  <- max(f_min_plot, as.numeric(stats::quantile(scores_df$frac_plot, 0.005, na.rm = TRUE)))
y_high <- as.numeric(stats::quantile(scores_df$frac_plot, 0.995, na.rm = TRUE))
if (!is.finite(y_low) || !is.finite(y_high) || y_low >= y_high) {
  y_low <- max(f_min_plot, min(scores_df$frac_plot, na.rm = TRUE))
  y_high <- max(scores_df$frac_plot, na.rm = TRUE)
}

subtitle_txt <- paste0(
  OUT_TAG,
  " | Xi p=", signif(xi_p, 3),
  " | NB LRT p=", signif(nb_fit$p_lrt, 3),
  " | FC(periphery/core)=", signif(fc_peri_core, 3),
  " | avgR=", signif(avgR, 3),
  if (!is.na(trend$peak_location)) paste0(" | ", trend$peak_location) else "",
  if (is.finite(VIS_SUBSAMPLE_N) && nrow(scores_df) > VIS_SUBSAMPLE_N)
    paste0(" | background=", VIS_SUBSAMPLE_N, " cells (random)")
  else ""
)

p <- ggplot() +
  geom_point(
    data = vis_df,
    aes(x = r, y = frac_plot),
    color = BG_POINT_COLOR,
    alpha = BG_POINT_ALPHA,
    size  = BG_POINT_SIZE
  ) +
  geom_errorbar(
    data = bin_df,
    aes(x = r_center, ymin = ci_low_plot, ymax = ci_high_plot),
    inherit.aes = FALSE,
    width = ERRBAR_WIDTH,
    linewidth = ERRBAR_LINEWIDTH,
    color = ERRBAR_COLOR
  ) +
  geom_point(
    data = bin_df,
    aes(x = r_center, y = mean_plot),
    inherit.aes = FALSE,
    shape  = 21,
    fill   = BIN_POINT_FILL,
    color  = BIN_POINT_EDGE,
    stroke = BIN_POINT_STROKE,
    size   = BIN_POINT_SIZE
  ) +
  geom_line(
    data = fit_curve,
    aes(x = r, y = frac_fit_plot),
    inherit.aes = FALSE,
    color = "black",
    linewidth = FIT_LINE_WIDTH
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    trans  = "log10",
    limits = c(y_low, y_high),
    oob    = scales::oob_squish,
    breaks = scales::breaks_log(n = 4),
    labels = scales::label_scientific(digits = 1)
  ) +
  labs(
    title    = trend$trend_class,
    subtitle = subtitle_txt,
    x = "Radial position (0 = core, 1 = periphery)",
    y = "Gene-set fractional abundance (sum UMIs in set / total UMIs)"
  ) +
  theme_nature()

print(p)

# Cairo fallback
use_cairo <- isTRUE(capabilities("cairo"))
device_fun <- if (use_cairo) grDevices::cairo_pdf else grDevices::pdf

ggsave(
  filename = FIG_PDF,
  plot = p,
  width = 5.8, height = 4.3, units = "in",
  device = device_fun,
  dpi = 600
)

cat("\n==================== Gene-set radial trend ====================\n")
cat("Tag:                           ", OUT_TAG, "\n", sep = "")
cat("Counts:                         ", COUNTS_RDS, "\n", sep = "")
cat("FACS:                           ", FACS_CSV, "\n", sep = "")
cat("Experiment:                     ", EXPERIMENT, "\n", sep = "")
cat("Cells used (analysis):          ", nrow(scores_df), "\n", sep = "")
if (nrow(vis_df) != nrow(scores_df)) cat("Cells plotted (background):     ", nrow(vis_df), " (random subset)\n", sep = "")
cat("Genes present:                  ", length(genes_present), "\n", sep = "")
cat("Genes missing:                  ", length(genes_missing), "\n", sep = "")
cat("Xi correlation p-value:         ", formatC(xi_p, digits = 4, format = "g"), "\n", sep = "")
cat("NB LRT p-value (quad vs null):  ", formatC(nb_fit$p_lrt, digits = 4, format = "g"), "\n", sep = "")
cat("FC(periphery/core):             ", formatC(fc_peri_core, digits = 3, format = "f"),
    " (log2FC=", formatC(log2fc_peri_core, digits = 3, format = "f"), ")\n", sep = "")
cat("avgR (0 core -> 1 periphery):   ", formatC(avgR, digits = 4, format = "f"), "\n", sep = "")
cat("Entropy H:                      ", formatC(entropy_H, digits = 4, format = "f"), "\n", sep = "")
cat("Session info saved:             ", normalizePath(SESSIONINFO_TXT), "\n", sep = "")
cat("Figure saved:                   ", normalizePath(FIG_PDF), "\n", sep = "")
cat("Summary saved:                  ", normalizePath(SUMMARY_CSV), "\n", sep = "")
cat("Genes used saved:               ", normalizePath(GENES_USED_CSV), "\n", sep = "")
cat("===============================================================\n\n")

summary_tbl <- tibble::tibble(
  tag = OUT_TAG,
  genes_requested = paste(genes_requested, collapse = ", "),
  genes_present = paste(genes_present, collapse = ", "),
  genes_missing = paste(genes_missing, collapse = ", "),
  n_cells = nrow(scores_df),
  xi_p = xi_p,
  nb_lrt_p = nb_fit$p_lrt,
  beta_r = nb_fit$beta_r,
  beta_r2 = nb_fit$beta_r2,
  r_star = trend$r_star,
  trend_class = trend$trend_class,
  peak_location = trend$peak_location,
  fc_periphery_core = fc_peri_core,
  log2fc_periphery_core = log2fc_peri_core,
  avgR = avgR,
  entropy_H = entropy_H,
  min_bin_n = MIN_BIN_N
)
readr::write_csv(summary_tbl, SUMMARY_CSV)
