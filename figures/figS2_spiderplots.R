# Supplementary Fig S2: one spider plot per Fig 4 field (54 fields, 5,606 tracks detected in >= 40 of
# 48 frames), every track shifted to start at (0, 0), drawn as Fig 4B; proCAF block on top, restCAF below.
# Drawn at final size; axes in px (0.64 um/px), same range in every panel; dy on image rows, not flipped.
# Out: results/thesis/suppfig2/suppfig2_spiderplots.pdf (+ 200-dpi PNG preview), suppfig2_fields.csv

suppressPackageStartupMessages({ library(ggplot2); library(grid); library(gtable); library(readxl) })

# USER INPUTS
ROOT      <- "."
TRAJ_DIR  <- file.path(ROOT, "results/guidance/traj")                          # from guidance/02_spots_to_trajectories.py
META_XLSX <- file.path(ROOT, "data/data_tracking/analysis/tracking_metadata.xlsx")
GPF_CSV   <- file.path(ROOT, "results/guidance/guidance_per_field.csv")        # Fig 4 per-field n_tracks
OUT_DIR   <- file.path(ROOT, "results/thesis/suppfig2")

MIN_FRAMES <- 40                       # a track must be detected in >= 40 of the 48 frames (Fig 4)
BAD_FIELDS <- "P170713_field_13"       # shares 27 of 48 raw frames with field_12; excluded everywhere
AX_LIM     <- 2200                     # px, half-width of every panel; must exceed the largest excursion
AX_BREAKS  <- c(-1000, 0, 1000)        # px; 0 gets a tick but no label (the labels would collide)
AX_LABELS  <- c("-1000", "", "1000")
TRACK_ALPHA <- 0.7                     # as Fig 4B
TUMOR_ORDER <- c("P201015", "PancT6", "P110616", "P170713")   # basal-like, basal-like, classical, classical

PANEL      <- 47                       # pt, side of every square panel; largest that fits 8 columns + labels in 6.5 in
GAP_BLOCK  <- 10                       # pt between the two blocks
# text in PDF pt at final size: ticks, strips, track counts 6.5; axis and block titles 7
PT_TICKLAB <- 6.5;  PT_TITLE <- 7.0;  PT_STRIP <- 6.5;  PT_N <- 6.5
PT_AXIS    <- 0.59; PT_TICK <- 0.36;  PT_TICK_LEN <- 0.95;  PT_BORDER <- 0.36
PT_TRACK   <- 0.25                     # line width of one track

PRO  <- "#E74095";  REST <- "#288D92"
FONT <- "Arial"                        # cairo substitutes Liberation Sans (metric-identical) on the cluster

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# unit conversions (exact)
LW <- 96 / (72 * .pt)                  # PDF pt of line width -> ggplot linewidth
FS <- 72.27 / 72                       # PDF pt of font size  -> ggplot size

# STEP 1. Every field's trajectories -> kept tracks, each shifted to start at (0, 0)
files <- list.files(TRAJ_DIR, pattern = "^trajectories_.*\\.csv$", recursive = TRUE, full.names = TRUE)
stopifnot(length(files) == 55)

pieces <- vector("list", length(files))
for (i in seq_along(files)) {
  field    <- sub("^trajectories_(.*)\\.csv$", "\\1", basename(files[i]))
  category <- basename(dirname(files[i]))                              # pro_basal, rest_classical, ...
  if (field %in% BAD_FIELDS) next
  d <- read.csv(files[i])
  stopifnot(all(c("particle", "frame", "x", "y") %in% names(d)), !anyNA(d[, c("particle", "frame", "x", "y")]))
  d <- d[order(d$particle, d$frame), ]
  stopifnot(!any(duplicated(d[, c("particle", "frame")])))              # one position per track per frame
  n_det <- ave(d$frame, d$particle, FUN = length)                       # frames in which the track was detected
  d <- d[n_det >= MIN_FRAMES, ]
  first <- !duplicated(d$particle)                                      # rows are sorted: first row = first frame
  x0 <- d$x[first][match(d$particle, d$particle[first])]
  y0 <- d$y[first][match(d$particle, d$particle[first])]
  pieces[[i]] <- data.frame(field = field, category = category, track = paste(field, d$particle),
                            frame = d$frame,
                            particle = d$particle,
                            dx = d$x - x0,                                # px, + = right
                            dy = d$y - y0,                                # px, image rows (y down), not flipped, as Fig 4B
                            stringsAsFactors = FALSE)
}
tr <- do.call(rbind, pieces)

# the Fig 4 track set, field by field
gpf   <- read.csv(GPF_CSV, stringsAsFactors = FALSE)
n_trk <- table(tr$field[!duplicated(tr$track)])                       # tracks per field
stopifnot(length(n_trk) == 54, setequal(names(n_trk), gpf$field),
          all(as.integer(n_trk[gpf$field]) == gpf$n_tracks), sum(n_trk) == 5606)
# nothing may fall outside the common axis range
stopifnot(max(abs(tr$dx)) < AX_LIM, max(abs(tr$dy)) < AX_LIM)
cat(sprintf("STEP 1  %s tracks in %d fields = the Fig 4 set (per-field counts match guidance_per_field.csv)\n",
            format(sum(n_trk), big.mark = ","), length(n_trk)))
cat(sprintf("        largest excursion from the start: %.0f px (%.2f mm); axis range +/- %d px\n",
            max(abs(c(tr$dx, tr$dy))), max(abs(c(tr$dx, tr$dy))) * 0.64 / 1000, AX_LIM))

# STEP 2. Field -> cfECM donor, CAF subtype, tumor line (join on the donor string, not the field number)
meta <- as.data.frame(read_excel(META_XLSX, sheet = "Sheet1"))
meta$field <- sprintf("%s_field_%02d", meta$Tumor, as.integer(meta$Field))
stopifnot(!any(duplicated(meta$field)), all(unique(tr$field) %in% meta$field))
fld <- data.frame(field = sort(unique(tr$field)), stringsAsFactors = FALSE)
fld$tumor     <- meta$Tumor[match(fld$field, meta$field)]
fld$donor     <- meta$cfECM[match(fld$field, meta$field)]
fld$cafType   <- meta$cafType[match(fld$field, meta$field)]
fld$tumorType <- meta$tumorType[match(fld$field, meta$field)]
fld$category  <- tr$category[match(fld$field, tr$field)]
fld$n_tracks  <- as.integer(n_trk[fld$field])
# the folder a field sits in must agree with the metadata, for subtype and tumor type
stopifnot(!anyNA(fld$donor), all(sub("_.*", "", fld$category) == sub("CAF$", "", fld$cafType)),
          all(sub(".*_", "", fld$category) == fld$tumorType))
stopifnot(setequal(fld$tumor, TUMOR_ORDER),
          all(fld$tumorType[fld$tumor %in% TUMOR_ORDER[1:2]] == "basal"),
          all(fld$tumorType[fld$tumor %in% TUMOR_ORDER[3:4]] == "classical"))
# de-identified alias used on every figure: P201015T1 -> P015T1, P220804N1 -> P804N1
fld$alias <- paste0("P", substring(fld$donor, 5))
stopifnot(all(nchar(fld$alias) == 6), length(unique(fld$donor)) == 14, length(unique(fld$alias)) == 14)
# one field per donor x tumor line
stopifnot(!any(duplicated(fld[, c("donor", "tumor")])))
ext <- tapply(pmax(abs(tr$dx), abs(tr$dy)), tr$field, max)
fld$max_excursion_px <- round(as.numeric(ext[fld$field]), 1)                 # largest |dx| or |dy| of any track

tr <- merge(tr, fld[, c("field", "tumor", "alias", "cafType")], by = "field")
cat("STEP 2  donors per subtype:", paste(names(table(fld$cafType[!duplicated(fld$donor)])),
                                        table(fld$cafType[!duplicated(fld$donor)]), collapse = ", "), "\n")

# every donor x tumor-line cell of the grid that has no field: say why
grid_all <- merge(unique(fld[, c("alias", "donor", "cafType")]), data.frame(tumor = TUMOR_ORDER), by = NULL)
grid_all$key <- paste(grid_all$donor, grid_all$tumor)
missing <- grid_all[!grid_all$key %in% paste(fld$donor, fld$tumor), ]
meta_key <- paste(meta$cfECM, meta$Tumor)
missing$why <- ifelse(paste(missing$donor, missing$tumor) %in% meta_key[meta$field %in% BAD_FIELDS],
                      "excluded", "no field")
cat("        empty grid cells:\n"); print(missing[, c("alias", "tumor", "cafType", "why")], row.names = FALSE)

# STEP 3. Facet order, labels, one color per track
TUMOR_LAB <- c(P201015 = "P201015\nbasal-like", PancT6 = "PancT6\nbasal-like",
               P110616 = "P110616\nclassical",   P170713 = "P170713\nclassical")
tr$tumor       <- factor(tr$tumor, levels = TUMOR_ORDER)
missing$tumor  <- factor(missing$tumor, levels = TUMOR_ORDER)
fld$tumor_f    <- factor(fld$tumor, levels = TUMOR_ORDER)

# colors as ggplot draws colour = factor(track) for one field (as Fig 4B): the default hue
# palette spread over that field's tracks, in track-ID order
tr$col <- NA_character_
for (fd in unique(tr$field)) {
  ids <- sort(unique(tr$particle[tr$field == fd]))
  pal <- scales::hue_pal()(length(ids))
  idx <- tr$field == fd
  tr$col[idx] <- pal[match(tr$particle[idx], ids)]
}
stopifnot(!anyNA(tr$col))

# STEP 4. Theme (absolute sizes) and one plot per CAF subtype
THEME <- theme_classic(base_size = PT_TITLE * FS, base_family = FONT) +
  theme(axis.line         = element_line(linewidth = PT_AXIS * LW, colour = "black", lineend = "square"),
        axis.ticks        = element_line(linewidth = PT_TICK * LW, colour = "grey20"),
        axis.ticks.length = unit(PT_TICK_LEN, "bigpts"),
        axis.text         = element_text(size = PT_TICKLAB * FS, colour = "black"),
        axis.text.x       = element_text(margin = margin(t = 1.5)),
        axis.text.y       = element_text(margin = margin(r = 1.5)),
        axis.title        = element_text(size = PT_TITLE * FS, colour = "black"),
        axis.title.x      = element_text(margin = margin(t = 2.5)),
        axis.title.y      = element_text(margin = margin(r = 2.5)),
        plot.title        = element_text(size = PT_TITLE * FS, face = "bold", hjust = 0, margin = margin(b = 3)),
        plot.title.position = "panel",
        strip.background  = element_blank(),
        strip.text.x      = element_text(size = PT_STRIP * FS, colour = "black", margin = margin(b = 2)),
        strip.text.y      = element_text(size = PT_STRIP * FS, colour = "black", angle = 0, hjust = 0,
                                         lineheight = 0.95, margin = margin(l = 3)),
        panel.border      = element_rect(fill = NA, colour = "grey55", linewidth = PT_BORDER * LW),
        panel.grid.major  = element_line(colour = "grey92", linewidth = PT_BORDER * LW),
        panel.spacing     = unit(3, "bigpts"),
        legend.position   = "none",
        plot.margin       = margin(2, 2, 2, 2),
        plot.background   = element_blank(),
        panel.background  = element_blank())

blocks <- list()
for (sub_type in c("proCAF", "restCAF")) {
  d  <- tr[tr$cafType == sub_type, ]
  f  <- fld[fld$cafType == sub_type, ]
  f$tumor <- factor(f$tumor, levels = TUMOR_ORDER)
  mi <- missing[missing$cafType == sub_type, ]
  n_don <- length(unique(f$donor))
  p <- ggplot(d, aes(dx, dy, group = track, colour = col)) +
    geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed", linewidth = PT_BORDER * LW) +
    geom_vline(xintercept = 0, colour = "grey50", linetype = "dashed", linewidth = PT_BORDER * LW) +
    geom_path(linewidth = PT_TRACK * LW, alpha = TRACK_ALPHA, lineend = "round") +   # joins consecutive detections (TrackMate closes gaps of up to 4 frames)
    # number of tracks in the field, top-left corner
    geom_text(data = f, aes(x = -AX_LIM * 0.93, y = AX_LIM * 0.93, label = n_tracks), inherit.aes = FALSE,
              hjust = 0, vjust = 1, size = PT_N * FS / .pt, family = FONT, colour = "black") +
    scale_colour_identity() +
    scale_x_continuous(limits = c(-AX_LIM, AX_LIM), breaks = AX_BREAKS, labels = AX_LABELS, expand = c(0, 0)) +
    scale_y_continuous(limits = c(-AX_LIM, AX_LIM), breaks = AX_BREAKS, labels = AX_LABELS, expand = c(0, 0)) +
    facet_grid(tumor ~ alias, labeller = labeller(tumor = TUMOR_LAB), drop = FALSE) +
    labs(x = "Δx (px)", y = "Δy (px)",
         title = sprintf("%s matrix (%d cfECM donors)", sub_type, n_don)) +
    THEME + theme(plot.title = element_text(colour = ifelse(sub_type == "proCAF", PRO, REST)))
  if (nrow(mi) > 0)
    p <- p + geom_label(data = mi, aes(x = 0, y = 0, label = why), inherit.aes = FALSE,     # white box hides the dashed lines
                        size = PT_N * FS / .pt, family = FONT, colour = "grey45", fill = "white",
                        label.size = 0, label.padding = unit(1, "bigpts"))
  blocks[[sub_type]] <- p
}

# STEP 5. Force every panel to PANEL x PANEL pt, stack the two blocks, write
# text is measured with a device open, so open a cairo device writing to nowhere before
# ggplotGrob(): the default pdf() device does not know the font and leaves a stray Rplots.pdf
cairo_pdf(nullfile(), family = FONT)
gts <- list(); wh <- list()
for (k in names(blocks)) {
  g   <- ggplotGrob(blocks[[k]])
  loc <- g$layout[grepl("^panel", g$layout$name), ]
  for (l in unique(loc$l)) g$widths[l]  <- unit(PANEL, "bigpts")
  for (t in unique(loc$t)) g$heights[t] <- unit(PANEL, "bigpts")
  stopifnot(length(unique(loc$t)) == 4)                                  # four tumor lines in every block
  gts[[k]] <- g
  wh[[k]]  <- c(w = sum(convertWidth(g$widths, "bigpts", valueOnly = TRUE)),
                h = sum(convertHeight(g$heights, "bigpts", valueOnly = TRUE)))
}
invisible(dev.off())
page_w <- max(wh$proCAF["w"], wh$restCAF["w"])                           # bigpts
page_h <- wh$proCAF["h"] + GAP_BLOCK + wh$restCAF["h"]
cat(sprintf("STEP 5  block widths: proCAF %.1f pt, restCAF %.1f pt; limit %.0f pt\n", wh$proCAF["w"], wh$restCAF["w"], 6.5 * 72))
stopifnot(page_w <= 6.5 * 72)                                            # fits the 1-in frame of a letter page

out_pdf <- file.path(OUT_DIR, "suppfig2_spiderplots.pdf")
for (dev_kind in c("pdf", "png")) {
  if (dev_kind == "pdf") cairo_pdf(out_pdf, width = page_w / 72, height = page_h / 72, family = FONT)
  if (dev_kind == "png") png(sub("\\.pdf$", "_200dpi.png", out_pdf), width = page_w / 72, height = page_h / 72,
                             units = "in", res = 200, type = "cairo", family = FONT, bg = "white")
  grid.newpage()
  # both blocks are left-aligned; their y-axis columns have the same width, so the panels line up
  pushViewport(viewport(x = 0, y = unit(1, "npc"), width = unit(wh$proCAF["w"], "bigpts"),
                        height = unit(wh$proCAF["h"], "bigpts"), just = c("left", "top")))
  grid.draw(gts$proCAF); popViewport()
  pushViewport(viewport(x = 0, y = 0, width = unit(wh$restCAF["w"], "bigpts"),
                        height = unit(wh$restCAF["h"], "bigpts"), just = c("left", "bottom")))
  grid.draw(gts$restCAF); popViewport()
  invisible(dev.off())
}
cat(sprintf("STEP 5  page %.1f x %.1f pt (%.2f x %.2f in), panels %d x %d pt -> %s\n",
            page_w, page_h, page_w / 72, page_h / 72, PANEL, PANEL, out_pdf))

# STEP 6. One row per panel, for the legend
fld <- fld[order(fld$cafType, fld$alias, fld$tumor_f), ]
write.csv(fld[, c("field", "cafType", "alias", "donor", "tumor", "tumorType", "n_tracks", "max_excursion_px")],
          file.path(OUT_DIR, "suppfig2_fields.csv"), row.names = FALSE)
cat(sprintf("STEP 6  tracks per field %d-%d (median %d); max excursion per field %.0f-%.0f px\n",
            min(fld$n_tracks), max(fld$n_tracks), as.integer(median(fld$n_tracks)),
            min(fld$max_excursion_px), max(fld$max_excursion_px)))
cat("written to", OUT_DIR, "\n")
