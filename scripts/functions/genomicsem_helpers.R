
# --------------------------------------------------
# Calculate pairwise rg, SE and p-value from S and V
# --------------------------------------------------
S <- covstruc$S
V <- covstruc$V

# Assign trait names to rows
if (is.null(rownames(S))) {
  rownames(S) <- colnames(S)
}

# Validate
stopifnot(
  identical(rownames(S), colnames(S)),
  isSymmetric(S),
  nrow(V) == nrow(S) * (nrow(S) + 1) / 2
)

# Index lower triangle of S, column by column
# Assumes GenomicSEM's standard vech ordering
idx <- which(
  lower.tri(S, diag = TRUE),
  arr.ind = TRUE
)

k <- nrow(S)

vech_index <- matrix(
  NA_integer_,
  nrow = k,
  ncol = k,
  dimnames = dimnames(S)
)

for (a in seq_len(nrow(idx))) {
  
  i <- idx[a, 1]
  j <- idx[a, 2]
  
  vech_index[i, j] <- a
  vech_index[j, i] <- a
}

get_rg <- function(trait1, trait2) {
  
  # Genetic variances and covariance
  a <- S[trait1, trait1]
  b <- S[trait2, trait2]
  c <- S[trait1, trait2]
  
  if (!all(is.finite(c(a, b, c))) || a <= 0 || b <= 0) {
    return(data.frame(
      rg = NA_real_,
      se = NA_real_,
      p = NA_real_
    ))
  }
  
  # Genetic correlation
  rg <- c / sqrt(a * b)
  
  # Corresponding positions in V
  ids <- c(
    vech_index[trait1, trait2],
    vech_index[trait1, trait1],
    vech_index[trait2, trait2]
  )
  
  # Delta-method gradient
  gradient <- c(
    1 / sqrt(a * b),
    -rg / (2 * a),
    -rg / (2 * b)
  )
  
  # Sampling variance of rg
  variance <- as.numeric(
    t(gradient) %*%
      V[ids, ids, drop = FALSE] %*%
      gradient
  )
  
  se <- if (is.finite(variance) && variance > 0) {
    sqrt(variance)
  } else {
    NA_real_
  }
  
  # Two-sided Wald p-value
  p <- if (is.finite(se) && se > 0) {
    2 * pnorm(-abs(rg / se))
  } else {
    NA_real_
  }
  
  data.frame(
    rg = rg,
    se = se,
    p = p
  )
}



# --------------------------------------------------
# Generate pairwise rg results for a trait family
# --------------------------------------------------

make_rg_data <- function(traits) {
  
  traits <- unique(traits)
  
  # Check that all traits exist
  missing <- setdiff(traits, colnames(S))
  
  if (length(missing) > 0) {
    stop(
      "Missing traits: ",
      paste(missing, collapse = ", ")
    )
  }
  
  # All unique trait pairs (excluding diagonal)
  pairs <- combn(traits, 2, simplify = FALSE)
  
  # Calculate rg for every pair
  results <- purrr::map_dfr(
    pairs,
    function(pair) {
      
      estimates <- get_rg(pair[1], pair[2])
      
      dplyr::bind_cols(
        data.frame(
          trait1 = pair[1],
          trait2 = pair[2]
        ),
        estimates
      )
    }
  )
  
  # BH FDR correction within this family
  results$p_fdr <- p.adjust(
    results$p,
    method = "BH"
  )
  
  results
}

#####################
### PLOT HEATMAPS ###
#####################

# Palette
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

# Use the supplied Matplotlib coolwarm palette
coolwarm <- COOLWARM

library(ComplexHeatmap)
library(circlize)
library(grid)

plot_family_heatmap <- function(
    family,
    traits,
    results,
    S,
    output_dir = "results/genomic_sem/heatmaps",
    cell_size_mm = 13
) {
  
  # ---------------------------------------
  # 1. Genetic correlation matrix
  # ---------------------------------------
  
  if (is.null(rownames(S))) {
    rownames(S) <- colnames(S)
  }
  
  stopifnot(all(traits %in% colnames(S)))
  
  rg_mat <- cov2cor(
    S[traits, traits, drop = FALSE]
  )
  
  if (any(!is.finite(rg_mat)) ||
      any(abs(rg_mat) > 1 + 1e-8)) {
    stop("Invalid or out-of-range genetic correlations.")
  }
  
  # ---------------------------------------
  # 2. FDR significance annotations
  # ---------------------------------------
  
  stars <- matrix(
    "",
    nrow = length(traits),
    ncol = length(traits),
    dimnames = dimnames(rg_mat)
  )
  
  for (i in seq_len(nrow(results))) {
    
    a <- results$trait1[i]
    b <- results$trait2[i]
    q <- results$p_fdr[i]
    
    sig <- if (is.na(q)) {
      ""
    } else if (q < 0.001) {
      "***"
    } else if (q < 0.01) {
      "**"
    } else if (q < 0.05) {
      "*"
    } else {
      ""
    }
    
    stars[a, b] <- sig
    stars[b, a] <- sig
  }
  
  labels <- matrix(
    paste0(sprintf("%.2f", rg_mat), stars),
    nrow = nrow(rg_mat),
    dimnames = dimnames(rg_mat)
  )
  
  # ---------------------------------------
  # 3. Exact supplied coolwarm palette
  # ---------------------------------------
  
  color_fun <- circlize::colorRamp2(
    seq(-1, 1, length.out = length(COOLWARM)),
    COOLWARM,
    space = "RGB"
  )
  
  # ---------------------------------------
  # 4. Hierarchical clustering
  # ---------------------------------------
  
  distance <- as.dist(1 - rg_mat)
  
  hc <- hclust(
    distance,
    method = "average"
  )
  
  # ---------------------------------------
  # 5. Heatmap dimensions
  # ---------------------------------------
  
  n <- length(traits)
  
  grid_size <- grid::unit(
    n * cell_size_mm,
    "mm"
  )
  
  # ---------------------------------------
  # 6. Build heatmap
  # ---------------------------------------
  
  ht <- ComplexHeatmap::Heatmap(
    rg_mat,
    
    name = "Genetic correlation",
    
    col = color_fun,
    
    # Cluster but hide dendrograms
    cluster_rows = hc,
    cluster_columns = hc,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    
    # Trait names
    row_names_side = "left",
    column_names_side = "bottom",
    
    row_names_gp = grid::gpar(fontsize = 10),
    column_names_gp = grid::gpar(fontsize = 10),
    
    column_names_rot = 45,
    
    # Grid appearance
    rect_gp = grid::gpar(
      col = "white",
      lwd = 0.8
    ),
    
    width = grid_size,
    height = grid_size,
    
    # Display rg and significance
    cell_fun = function(j, i, x, y, width, height, fill) {
      
      grid::grid.text(
        labels[i, j],
        x = x,
        y = y,
        gp = grid::gpar(
          fontsize = 8,
          col = "black"
        )
      )
      
    },
    
    # Legend
    heatmap_legend_param = list(
      title = expression(r[g]),
      at = seq(-1, 1, by = 0.25),
      labels = sprintf(
        "%.2f",
        seq(-1, 1, by = 0.25)
      ),
      legend_height = grid_size,
      grid_width = grid::unit(7, "mm"),
      title_gp = grid::gpar(
        fontsize = 11,
        fontface = "bold"
      ),
      labels_gp = grid::gpar(fontsize = 9)
    ),
    
    show_heatmap_legend = TRUE
  )
  
  # ---------------------------------------
  # 7. Save figures
  # ---------------------------------------
  
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  file_base <- file.path(
    output_dir,
    paste0(family, "_rg_heatmap")
  )
  
  # Allow room for names and legend
  fig_width <- n * cell_size_mm / 25.4 + 4
  fig_height <- n * cell_size_mm / 25.4 + 3
  
  # PNG
  png(
    paste0(file_base, ".png"),
    width = fig_width,
    height = fig_height,
    units = "in",
    res = 300
  )
  
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "right"
  )
  
  dev.off()
  
  # Vector PDF
  pdf(
    paste0(file_base, ".pdf"),
    width = fig_width,
    height = fig_height
  )
  
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "right"
  )
  
  dev.off()
  
  invisible(ht)
}
