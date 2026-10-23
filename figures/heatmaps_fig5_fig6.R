# Fig 6B (23 matrisome hits x 12 proteomics donors) and Fig 5F (10 CAF markers x 12 lines) heatmaps,
# drawn at final size with one set of pheatmap settings; columns run from most to least aligned matrix.
# In:  results/proteomics/run_old_keyence/, results/thesis/fig5/fig5F_marker_input.csv (figures/fig5_panels.R)
# Out: results/thesis/fig6/fig6B_matrisome_heatmap.pdf, results/thesis/fig5/fig5F_marker_heatmap.pdf, z-score csvs

suppressPackageStartupMessages({ library(pheatmap); library(grid); library(RColorBrewer) })

# USER INPUTS
ROOT     <- "."
RUN_DIR  <- file.path(ROOT, "results/proteomics/run_old_keyence")        # Fig 6 run (design.rds, matrisome_intersection.csv)
F5_INPUT <- file.path(ROOT, "results/thesis/fig5/fig5F_marker_input.csv") # from figures/fig5_panels.R
OUT6     <- file.path(ROOT, "results/thesis/fig6")
OUT5     <- file.path(ROOT, "results/thesis/fig5")

# sizes in PDF pt at final size
CELL    <- 8       # heatmap cell, square
TREE_H  <- 12      # row dendrogram width
LEG_W   <- 6       # color bar width
LEG_H   <- 56      # color bar height (pheatmap's own bar is a fixed 150 pt and overflows the 80-pt Fig 5F)
PT_TICK <- 0.36    # color bar frame and tick line width
PT_LAB  <- 6.5     # row / column labels, legend text, annotation names (= tick-label size elsewhere)
Z_CLIP  <- 2       # color scale runs -Z_CLIP .. +Z_CLIP
PSEUDO  <- 1       # added to marker intensities before log2

DIV_COLORS   <- c("Core matrisome" = "#5E4B8B", "Matrisome-associated" = "#D9A441")   # as proteomics/04_figures.R
BATCH_COLORS <- c(PC1514 = "#BDBDBD", PC1556 = "#6E6E6E")
GROUP_COLORS <- c(proCAF = "#E74095", restCAF = "#288D92")
SD_COLORS    <- brewer.pal(9, "Greys")[2:7]
HEAT_COLORS  <- colorRampPalette(rev(brewer.pal(7, "RdBu")))(100)
FONT <- "Arial"

# de-identified aliases, the same 17 as proteomics/04_figures.R
ALIAS <- c(P201015T1 = "P015T1", P220804N1 = "P804N1", P220804T1 = "P804T1",
           P220815T1 = "P815T1", P220908N1 = "P908N1", P231207T1 = "P207T1",
           P240229T1 = "P229T1", P140807T1 = "P807T1", P210923T1 = "P923T1",
           P220201N1 = "P201N1", P220212T1 = "P212T1", P220519T1 = "P519T1",
           P221108T1 = "P108T1", P231207NA = "P207NA", P240711N1 = "P711N1",
           P220118N1 = "P118N1", P240509N1 = "P509N1")
stopifnot(all(ALIAS == paste0("P", substring(names(ALIAS), 5))))

dir.create(OUT6, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT5, showWarnings = FALSE, recursive = TRUE)
FS <- 72.27 / 72                                  # PDF pt -> grid "points" (1/72.27 in)

# STEP 1. Fig 6B data as in proteomics/04_figures.R STEP 3: donor means of log2 intensity, z-scored per protein
d       <- readRDS(file.path(RUN_DIR, "design.rds"))
samples <- d$samples
quant   <- d$mat_quant_use
ov      <- read.csv(file.path(RUN_DIR, "matrisome_intersection.csv"))
stopifnot(nrow(ov) == 23, identical(colnames(quant$E), samples$sample), all(ov$UNIPROT %in% rownames(quant$E)))

donor_ids <- levels(samples$donor)
E_sig   <- quant$E[ov$UNIPROT, , drop = FALSE]
E_donor <- sapply(donor_ids, function(dd) rowMeans(E_sig[, samples$donor == dd, drop = FALSE]))
info <- unique(samples[, c("donor", "org_sd_deg", "batch", "group")])
info <- info[match(donor_ids, info$donor), ]
ord6 <- order(info$org_sd_deg)
hm6  <- t(scale(t(E_donor[, ord6, drop = FALSE])))
rownames(hm6) <- ifelse(is.na(ov$Genes) | ov$Genes == "", ov$UNIPROT, ov$Genes)
stopifnot(all(colnames(hm6) %in% names(ALIAS)), ncol(hm6) == 12, nrow(hm6) == 23,
          all(abs(rowMeans(hm6)) < 1e-9))
colnames(hm6) <- ALIAS[colnames(hm6)]
col_annot6 <- data.frame(`CAF type` = info$group[ord6], Batch = info$batch[ord6],
                         `SD around 0 (deg)` = info$org_sd_deg[ord6],
                         row.names = colnames(hm6), check.names = FALSE)
row_annot6 <- data.frame(Division = ov$Division, row.names = rownames(hm6))
cat(sprintf("STEP 1  Fig 6B: %d proteins x %d donors, SD %.1f -> %.1f deg left to right\n",
            nrow(hm6), ncol(hm6), min(info$org_sd_deg), max(info$org_sd_deg)))

# STEP 2. Fig 5F data: log2(intensity + 1), z-scored per marker
# the +1 keeps zeros (not detected: IL6 in 8 of 12 lines, FBLN5 in 3) as the lowest value of that marker
f5 <- read.csv(F5_INPUT, stringsAsFactors = FALSE)
MARKERS <- setdiff(names(f5), c("donor", "alias", "cafType", "sd_around_0_deg", "perc_within_15deg"))
stopifnot(nrow(f5) == 12, length(MARKERS) == 10, !is.unsorted(f5$sd_around_0_deg),
          all(f5$donor %in% names(ALIAS)), all(f5$alias == ALIAS[f5$donor]),
          all(as.matrix(f5[, MARKERS]) >= 0))
M   <- t(as.matrix(f5[, MARKERS]))                      # rows = markers, columns = lines
colnames(M) <- f5$alias
hm5 <- t(scale(t(log2(M + PSEUDO))))
stopifnot(all(abs(rowMeans(hm5)) < 1e-9), !anyNA(hm5))
rownames(hm5)[rownames(hm5) == "aSMA"] <- "αSMA"
col_annot5 <- data.frame(`CAF type` = f5$cafType, `SD around 0 (deg)` = f5$sd_around_0_deg,
                         row.names = colnames(hm5), check.names = FALSE)
cat(sprintf("STEP 2  Fig 5F: %d markers x %d lines, SD %.1f -> %.1f deg; zeros (not detected): %d cells\n",
            nrow(hm5), ncol(hm5), min(f5$sd_around_0_deg), max(f5$sd_around_0_deg), sum(M == 0)))

# STEP 3. The two heatmaps, one set of settings
# pheatmap measures text while it lays out the plot, so a cairo device that knows the font
# must be open (writing to nowhere); otherwise R opens pdf() and leaves a stray Rplots.pdf
ANN_COLORS <- list(Division = DIV_COLORS, `CAF type` = GROUP_COLORS, Batch = BATCH_COLORS,
                   `SD around 0 (deg)` = SD_COLORS)
BREAKS <- seq(-Z_CLIP, Z_CLIP, length.out = 101)

cairo_pdf(nullfile(), family = FONT)
ph6 <- pheatmap(hm6, annotation_col = col_annot6, annotation_row = row_annot6, annotation_colors = ANN_COLORS,
                cluster_cols = FALSE, cluster_rows = TRUE, clustering_distance_rows = "euclidean",
                clustering_method = "ward.D2", treeheight_row = TREE_H,
                color = HEAT_COLORS, breaks = BREAKS, legend_breaks = c(-2, -1, 0, 1, 2),
                border_color = NA, cellwidth = CELL, cellheight = CELL,
                fontsize = PT_LAB * FS, fontsize_row = PT_LAB * FS, fontsize_col = PT_LAB * FS,
                angle_col = "90", silent = TRUE)
ph5 <- pheatmap(hm5, annotation_col = col_annot5, annotation_colors = ANN_COLORS,
                cluster_cols = FALSE, cluster_rows = TRUE, clustering_distance_rows = "euclidean",
                clustering_method = "ward.D2", treeheight_row = TREE_H,
                color = HEAT_COLORS, breaks = BREAKS, legend_breaks = c(-2, -1, 0, 1, 2),
                border_color = NA, cellwidth = CELL, cellheight = CELL,
                fontsize = PT_LAB * FS, fontsize_row = PT_LAB * FS, fontsize_col = PT_LAB * FS,
                angle_col = "90", silent = TRUE)
invisible(dev.off())

# STEP 4. One shared color bar
# pheatmap's bar is a fixed 150 pt tall; replace it with a LEG_W x LEG_H bar, same 100 colors, ticks at -2..2
LEG_BREAKS <- c(-2, -1, 0, 1, 2)
tick_at <- (LEG_BREAKS + Z_CLIP) / (2 * Z_CLIP)
LWD     <- PT_TICK * 96 / 72                      # 1 R lwd = 1/96 in
bar   <- rectGrob(x = 0, y = seq(0, 1, length.out = 101)[-101], width = 1, height = 1 / 100,
                  just = c("left", "bottom"), gp = gpar(fill = HEAT_COLORS, col = NA))
frame <- rectGrob(gp = gpar(fill = NA, col = "black", lwd = LWD))
ticks <- segmentsGrob(x0 = unit(1, "npc"), x1 = unit(1, "npc") + unit(1.5, "bigpts"),
                      y0 = tick_at, y1 = tick_at, gp = gpar(col = "black", lwd = LWD))
labs  <- textGrob(LEG_BREAKS, x = unit(1, "npc") + unit(3, "bigpts"), y = tick_at, just = "left",
                  gp = gpar(fontsize = PT_LAB * FS))
title <- textGrob("z-score", x = 0, y = unit(1, "npc") + unit(3, "bigpts"), just = c("left", "bottom"),
                  gp = gpar(fontsize = PT_LAB * FS, fontface = "bold"))
colorbar <- gTree(children = gList(bar, frame, ticks, labs, title),
                  vp = viewport(x = unit(2, "bigpts"), y = unit(1, "npc") - unit(8, "bigpts"), just = c("left", "top"),
                                width = unit(LEG_W, "bigpts"), height = unit(LEG_H, "bigpts")))
gt6 <- ph6$gtable; gt5 <- ph5$gtable
i6 <- which(gt6$layout$name == "legend"); i5 <- which(gt5$layout$name == "legend")
stopifnot(length(i6) == 1, length(i5) == 1)
gt6$grobs[[i6]] <- colorbar
gt5$grobs[[i5]] <- colorbar
# the legend column is sized for pheatmap's own bar; give the "z-score" title room before the annotation legend
gt6$widths[gt6$layout$l[i6]] <- unit(LEG_W + 30, "bigpts")
gt5$widths[gt5$layout$l[i5]] <- unit(LEG_W + 30, "bigpts")
stopifnot(LEG_H + 8 <= CELL * nrow(hm5) + 2 * CELL)   # the bar must end above the bottom of the smaller heatmap
cairo_pdf(nullfile(), family = FONT)
w6 <- convertWidth(sum(gt6$widths), "in", valueOnly = TRUE); h6 <- convertHeight(sum(gt6$heights), "in", valueOnly = TRUE)
w5 <- convertWidth(sum(gt5$widths), "in", valueOnly = TRUE); h5 <- convertHeight(sum(gt5$heights), "in", valueOnly = TRUE)
invisible(dev.off())
stopifnot(is.finite(w6), is.finite(h6), is.finite(w5), is.finite(h5), w6 > 0, h6 > 0, w5 > 0, h5 > 0)

# STEP 5. Write each at its natural size (cells are CELL pt, nothing is rescaled)
out6 <- file.path(OUT6, "fig6B_matrisome_heatmap.pdf")
cairo_pdf(out6, width = w6, height = h6, family = FONT); grid.draw(gt6); invisible(dev.off())
out5 <- file.path(OUT5, "fig5F_marker_heatmap.pdf")
cairo_pdf(out5, width = w5, height = h5, family = FONT); grid.draw(gt5); invisible(dev.off())
write.csv(round(hm6, 3), file.path(OUT6, "fig6B_zscores.csv"))
write.csv(round(hm5, 3), file.path(OUT5, "fig5F_zscores.csv"))
cat(sprintf("Fig 6B  heatmap body %3d x %3d pt | page %5.1f x %5.1f pt  %s\n", CELL * ncol(hm6), CELL * nrow(hm6), w6 * 72, h6 * 72, out6))
cat(sprintf("Fig 5F  heatmap body %3d x %3d pt | page %5.1f x %5.1f pt  %s\n", CELL * ncol(hm5), CELL * nrow(hm5), w5 * 72, h5 * 72, out5))
cat("row order 6B:", paste(rownames(hm6)[ph6$tree_row$order], collapse = " "), "\n")
cat("row order 5F:", paste(rownames(hm5)[ph5$tree_row$order], collapse = " "), "\n")
