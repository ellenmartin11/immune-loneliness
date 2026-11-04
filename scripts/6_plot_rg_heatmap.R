#!/usr/bin/env Rscript
# plot_rg_heatmap.R -- LDSC genetic-correlation heatmap for one trait group
#
# Reads an rg matrix (from run_ldsr_rg.R's rg_matrix.tsv), subsets to the
# traits in a group file, and draws an annotated tile heatmap with
# FDR-corrected significance asterisks.
#
# Usage:
#   Rscript plot_rg_heatmap.R <rg_matrix.tsv> <group_traits.txt> <out_prefix> [rg_long.tsv] [title]
#
#   Rscript plot_rg_heatmap.R ldsr_rg/rg_matrix.tsv group1_psych_traits.txt fig_rg_psych
#   Rscript plot_rg_heatmap.R ldsr_rg/rg_matrix.tsv group1_psych_traits.txt fig_rg_psych ldsr_rg/rg_long.tsv "Psychiatric"
#
# If rg_long.tsv is omitted, the script looks for rg_long.tsv in the same
# directory as rg_matrix.tsv. FDR (Benjamini-Hochberg) is computed across the
# group's off-diagonal pairs. Asterisks: * q<0.05, ** q<0.01, *** q<0.001.
#
# Output: <out_prefix>.svg and <out_prefix>.png
# Env overrides: RG_WIDTH / RG_HEIGHT (inches, defaults 9 / 7.5)

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))

# SVG device: prefer svglite; fall back to base svg() (cairo); else PNG only
has_svglite <- requireNamespace("svglite", quietly = TRUE)
has_cairo   <- isTRUE(capabilities()[["cairo"]])

# matplotlib "coolwarm" colormap, 256 colors (matches seaborn cmap="coolwarm")
COOLWARM <- c(
  "#3A4CC0", "#3B4DC1", "#3C4FC3", "#3E51C4", "#3F53C6", "#4054C7",
  "#4156C9", "#4258CA", "#435ACC", "#455BCD", "#465DCF", "#475FD0",
  "#4860D1", "#4962D3", "#4B64D4", "#4C66D6", "#4D67D7", "#4E69D8",
  "#506BDA", "#516CDB", "#526EDC", "#5370DD", "#5571DE", "#5673E0",
  "#5775E1", "#5876E2", "#5A78E3", "#5B79E4", "#5C7BE5", "#5D7DE6",
  "#5F7EE7", "#6080E8", "#6182EA", "#6383EA", "#6485EB", "#6586EC",
  "#6788ED", "#6889EE", "#698BEF", "#6B8DF0", "#6C8EF1", "#6D90F1",
  "#6F91F2", "#7093F3", "#7194F4", "#7395F4", "#7497F5", "#7598F6",
  "#779AF6", "#789BF7", "#7A9DF8", "#7B9EF8", "#7CA0F9", "#7EA1F9",
  "#7FA2FA", "#80A4FA", "#82A5FB", "#83A6FB", "#85A8FB", "#86A9FC",
  "#87AAFC", "#89ACFC", "#8AADFD", "#8BAEFD", "#8DAFFD", "#8EB1FD",
  "#90B2FE", "#91B3FE", "#92B4FE", "#94B5FE", "#95B7FE", "#97B8FE",
  "#98B9FE", "#99BAFE", "#9BBBFE", "#9CBCFE", "#9DBDFE", "#9FBEFE",
  "#A0BFFE", "#A2C0FE", "#A3C1FE", "#A4C2FE", "#A6C3FD", "#A7C4FD",
  "#A8C5FD", "#AAC6FD", "#ABC7FC", "#ACC8FC", "#AEC9FC", "#AFCAFB",
  "#B0CBFB", "#B2CBFB", "#B3CCFA", "#B4CDFA", "#B6CEF9", "#B7CFF9",
  "#B8CFF8", "#B9D0F8", "#BBD1F7", "#BCD1F6", "#BDD2F6", "#BED3F5",
  "#C0D3F5", "#C1D4F4", "#C2D4F3", "#C3D5F2", "#C5D5F2", "#C6D6F1",
  "#C7D6F0", "#C8D7EF", "#C9D7EE", "#CAD8EE", "#CCD8ED", "#CDD9EC",
  "#CED9EB", "#CFD9EA", "#D0DAE9", "#D1DAE8", "#D2DAE7", "#D3DBE6",
  "#D5DBE5", "#D6DBE4", "#D7DBE2", "#D8DBE1", "#D9DCE0", "#DADCDF",
  "#DBDCDE", "#DCDCDD", "#DDDCDB", "#DEDBDA", "#DFDBD9", "#E0DAD7",
  "#E1DAD6", "#E2D9D4", "#E3D9D3", "#E4D8D1", "#E5D8D0", "#E6D7CF",
  "#E7D6CD", "#E7D6CC", "#E8D5CA", "#E9D4C9", "#EAD3C7", "#EBD3C6",
  "#ECD2C4", "#ECD1C3", "#EDD0C1", "#EDCFC0", "#EECFBE", "#EFCEBC",
  "#EFCDBB", "#F0CCB9", "#F1CBB8", "#F1CAB6", "#F2C9B5", "#F2C8B3",
  "#F2C7B2", "#F3C6B0", "#F3C5AF", "#F4C4AD", "#F4C3AB", "#F4C2AA",
  "#F5C1A8", "#F5C0A7", "#F5BFA5", "#F6BDA4", "#F6BCA2", "#F6BBA0",
  "#F6BA9F", "#F6B99D", "#F6B79C", "#F6B69A", "#F7B598", "#F7B397",
  "#F7B295", "#F7B194", "#F7B092", "#F7AE91", "#F7AD8F", "#F6AB8D",
  "#F6AA8C", "#F6A98A", "#F6A789", "#F6A687", "#F6A486", "#F6A384",
  "#F5A182", "#F5A081", "#F59E7F", "#F49D7E", "#F49B7C", "#F49A7B",
  "#F39879", "#F39678", "#F39576", "#F29375", "#F29173", "#F19072",
  "#F18E70", "#F08D6F", "#F08B6D", "#EF896C", "#EE876A", "#EE8669",
  "#ED8467", "#EC8266", "#EC8064", "#EB7F63", "#EA7D61", "#EA7B60",
  "#E9795E", "#E8775D", "#E7755C", "#E6745A", "#E67259", "#E57057",
  "#E46E56", "#E36C54", "#E26A53", "#E16852", "#E06650", "#DF644F",
  "#DE624E", "#DD604C", "#DC5E4B", "#DB5C4A", "#DA5A48", "#D95847",
  "#D85646", "#D75444", "#D65243", "#D44F42", "#D34D40", "#D24B3F",
  "#D1493E", "#CF463D", "#CE443C", "#CD423A", "#CC3F39", "#CA3D38",
  "#C93B37", "#C83835", "#C63534", "#C53233", "#C43032", "#C22D31",
  "#C12A30", "#BF282E", "#BE232D", "#BC1F2C", "#BB1A2B", "#B9162A",
  "#B81129", "#B60D28", "#B50827", "#B30326"
)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript plot_rg_heatmap.R <rg_matrix.tsv> <group_traits.txt> <out_prefix> [rg_long.tsv] [title]")
}
mat_file <- args[1]; group_file <- args[2]; out_prefix <- args[3]
long_file <- if (length(args) >= 4 && nzchar(args[4])) args[4] else
  file.path(dirname(mat_file), "rg_long.tsv")
plot_title <- if (length(args) >= 5 && nzchar(args[5])) args[5] else
  tools::toTitleCase(sub("_traits\\.txt$", "", basename(group_file)))

m  <- fread(mat_file)                    # column 1 = row trait names
gn <- names(m)[-1]
m  <- data.frame(m, check.names = FALSE, stringsAsFactors = FALSE)
rownames(m) <- m[[1]]; m[[1]] <- NULL
mat <- as.matrix(m)
colnames(mat) <- gn

traits <- trimws(readLines(group_file, warn = FALSE))
traits <- traits[nzchar(traits) & traits != "phenotype"]

missing <- setdiff(traits, rownames(mat))
if (length(missing)) {
  message("NOTE: not in the rg matrix (skipped): ", paste(missing, collapse = ", "))
  traits <- intersect(traits, rownames(mat))
}
if (length(traits) < 2) stop("Fewer than 2 group traits found in the matrix.")

sub <- mat[traits, traits, drop = FALSE]

# order by hierarchical clustering of the rg values
dmat <- as.dist(1 - abs(sub))
if (any(is.na(dmat))) {
  napairs <- which(is.na(sub), arr.ind = TRUE)
  message("NOTE: NA rg values - clustering skipped, file order used. NA pairs: ",
          paste(unique(paste(rownames(sub)[napairs[,1]], rownames(sub)[napairs[,2]], sep="-")),
                collapse = ", "))
  ord <- seq_len(nrow(sub))
} else {
  ord <- hclust(dmat)$order
}
sub <- sub[ord, ord, drop = FALSE]

# ---- FDR-corrected significance (BH over the group's off-diagonal pairs) ----
star_of <- function(q) ifelse(is.na(q), "", ifelse(q < 0.001, "***",
                                                   ifelse(q < 0.01, "**", ifelse(q < 0.05, "*", ""))))
qmat <- matrix(NA_real_, nrow(sub), ncol(sub), dimnames = dimnames(sub))
if (file.exists(long_file)) {
  long <- fread(long_file)
  setnames(long, tolower(names(long)))
  pr <- long[trait1 %in% traits & trait2 %in% traits & trait1 != trait2,
             .(a = pmin(trait1, trait2), b = pmax(trait1, trait2), p)]
  pr <- pr[!duplicated(pr, by = c("a", "b"))]
  qmat <- matrix(NA_real_, nrow(sub), ncol(sub), dimnames = dimnames(sub))
  if (nrow(pr)) {
    pr[, q := p.adjust(p, method = "BH")]
    for (i in seq_len(nrow(pr))) {
      if (pr$a[i] %in% rownames(qmat) && pr$b[i] %in% rownames(qmat)) {
        qmat[pr$a[i], pr$b[i]] <- qmat[pr$b[i], pr$a[i]] <- pr$q[i]
      }
    }
  }
} else {
  message("NOTE: rg_long.tsv not found - drawing without significance asterisks.")
}

plot_dt <- data.table(
  row = rep(rownames(sub), times = ncol(sub)),
  col = rep(colnames(sub), each = nrow(sub)),
  rg  = as.vector(sub),
  q   = as.vector(qmat)
)
plot_dt[, label := ifelse(is.na(rg), "", sprintf("%.2f", rg))]
plot_dt[row != col, label := paste0(label, star_of(q))]
# white text on dark tiles, black on light (seaborn-style)
plot_dt[, frac := pmin(pmax((rg + 1) / 2, 0), 1)]
plot_dt[, idx := pmin(as.integer(frac * 255) + 1, 256)]
plot_dt[, lum := sapply(idx, function(i) {
  cc <- col2rgb(COOLWARM[i]) / 255
  lum <- 0.299 * cc[1] + 0.587 * cc[2] + 0.114 * cc[3]
  ifelse(is.na(lum) || lum > 0.5, "#000000", "#FFFFFF")
})]

short <- function(x) sub(" isolation PanUKB|_PanUKB| PanUKB", "", x)
plot_dt[, row_s := short(as.character(row))]
plot_dt[, col_s := short(as.character(col))]
plot_dt[, row := factor(row_s, levels = short(rev(rownames(sub))))]
plot_dt[, col := factor(col_s, levels = short(colnames(sub)))]

w <- as.numeric(Sys.getenv("RG_WIDTH",  unset = 9))
h <- as.numeric(Sys.getenv("RG_HEIGHT", unset = 7.5))

p <- ggplot(plot_dt, aes(col, row, fill = rg)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = label, color = lum), size = as.numeric(Sys.getenv("RG_TEXT_SIZE", unset = 2))) +
  scale_color_identity() +
  scale_fill_gradientn(colours = COOLWARM, na.value = "grey85",
                       name = "rg", limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  scale_x_discrete(position = "top") +
  coord_fixed() +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0.2),
        panel.grid = element_blank(),
        plot.margin = margin(45, 20, 20, 20),
        plot.title = element_text(face = "bold", hjust = 0.5,
                                  margin = margin(b = 45)),
        text = element_text(family = "Liberation Sans")) +
  labs(x = NULL, y = NULL, title = plot_title,
       caption = "* q(FDR) < 0.05   ** q < 0.01   *** q < 0.001 (BH, within-panel pairs)")

if (has_svglite) {
  ggsave(paste0(out_prefix, ".svg"), p, width = w, height = h)
} else if (has_cairo) {
  ggsave(paste0(out_prefix, ".svg"), p, width = w, height = h, device = svg)
} else {
  message("NOTE: no svglite and no cairo - skipping SVG, writing PNG only.")
}
ggsave(paste0(out_prefix, ".png"), p, width = w, height = h, dpi = 300)
cat("Written:", paste0(out_prefix, ".svg"), "and", paste0(out_prefix, ".png"), "\n")

