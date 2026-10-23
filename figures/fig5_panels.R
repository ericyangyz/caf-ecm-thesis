# Fig 5B-E panels at final size, the Fig 5F heatmap input, and the Supplementary Fig S3 panels.
# B: DeCAF score vs ECM fiber spread; C: Young's modulus (median, mean) vs spread; D/E: linearity of
# forward progression and directional change rate by CAF subtype (cfECM donor means); S3: same per track.
# Out: results/thesis/fig5/ (panels, fig5F_marker_input.csv, fig5_stats.csv), results/thesis/suppfig3/

suppressPackageStartupMessages({ library(ggplot2); library(grid); library(readxl) })

# USER INPUTS
ROOT       <- "."
TM_DIR     <- file.path(ROOT, "data/data_tracking/data/by_CAF")            # TrackMate export, Pro/ and Rest/
META_XLSX  <- file.path(ROOT, "data/data_tracking/analysis/tracking_metadata.xlsx")
ECM_CSV    <- "./results/old_keyence_alignment_per_sample.csv"
DECAF_TSV  <- file.path(ROOT, "data/data_decaf/decaf_scores_20260922.tsv")
STIFF_XLSX <- file.path(ROOT, "data/data_stiffness/items.xlsx")
OUT_DIR    <- file.path(ROOT, "results/thesis/fig5")
SUPP3_DIR  <- file.path(ROOT, "results/thesis/suppfig3")                   # Supp Fig S3 (track-level) panels

MIN_SPOTS    <- 40                     # keep tracks detected in >= 40 of 48 frames (the Fig 4 filter)
MIN_DURATION <- 40                     # frames; the labmate's own filter (TRACK_DURATION > 40), read-check only
GPF_CSV      <- file.path(ROOT, "results/guidance/guidance_per_field.csv")   # Fig 4 per-field n_tracks
BAD_FIELDS   <- "P170713_field_13"     # shares 27 of 48 raw frames with field_12; excluded everywhere
RANDOM_DCR   <- pi / 2                 # rad/frame, mean turning angle of an unbiased random walk
CP_LINES_IN  <- FALSE                  # the 3 chronic-pancreatitis lines have no ECM data anyway

# panels are drawn at final size (PDF pt); lines and points as Fig 4E-H, text 6.5 pt, titles 7 pt
PANEL_H     <- 95      # y axis length of every box / scatter panel
PANEL_W_BOX <- 98      # x axis length, two-group boxplots (= Fig 4E)
PANEL_W_SCT <- 99      # x axis length, scatter panels     (= Fig 4F)
PT_AXIS     <- 0.59;  PT_TICK <- 0.36;  PT_TICK_LEN <- 0.95
PT_POINT    <- 2.05;  PT_FIT <- 0.59;   PT_BOX <- 0.345;   PT_MEDIAN <- 0.78;  PT_DASH <- 0.37
PT_TICKLAB  <- 6.5;   PT_TITLE <- 7.0;  PT_LEGEND <- 6.5;  PT_STAT <- 7.0;  PT_STRIP <- 7.0
PT_LEG_ROW  <- 8.4     # legend row pitch

PRO  <- "#E74095";  REST <- "#288D92"
CAF_COLS   <- c(proCAF = PRO, restCAF = REST)
CAF_SHAPES <- c(proCAF = 16, restCAF = 15)
FONT <- "Arial"                        # cairo substitutes Liberation Sans (metric-identical) on the cluster

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(SUPP3_DIR, showWarnings = FALSE, recursive = TRUE)

# unit conversions (exact)
LW <- 96 / (72 * .pt)                  # PDF pt of line width -> ggplot linewidth
FS <- 72.27 / 72                       # PDF pt of font size  -> ggplot size
PS <- FS / (0.75 * .pt)                # PDF pt of point diameter -> ggplot size (stroke = 0)

# STEP 1. Labmate's TrackMate tracks: 55 files -> one table
# linearity = mean straight-line speed / mean speed (dimensionless); dcr in rad per frame (30 min)
files <- c(list.files(file.path(TM_DIR, "Pro"),  pattern = "\\.csv$", full.names = TRUE),
           list.files(file.path(TM_DIR, "Rest"), pattern = "\\.csv$", full.names = TRUE))
stopifnot(length(files) == 55)

pieces <- vector("list", length(files))
for (i in seq_along(files)) {
  nm <- basename(files[i])
  # names are inconsistent: P110616_field_2_stitched..., P170713_field_4_stitch..., PancT6_field6_stitched...
  m <- regmatches(nm, regexec("^(P[0-9]+|PancT6)_field_?([0-9]+)_stitch(ed)?_mask_tracks\\.csv$", nm))[[1]]
  stopifnot(length(m) == 4)
  d <- read.csv(files[i], check.names = FALSE, stringsAsFactors = FALSE)
  d <- d[, !duplicated(names(d))]                        # TRACK_ID appears twice in the export
  stopifnot(all(c("TRACK_DURATION", "NUMBER_SPOTS", "LINEARITY_OF_FORWARD_PROGRESSION",
                  "MEAN_DIRECTIONAL_CHANGE_RATE") %in% names(d)))
  pieces[[i]] <- data.frame(
    field     = sprintf("%s_field_%02d", m[2], as.integer(m[3])),
    tumor     = m[2],
    cafType   = ifelse(grepl("/Pro/", files[i], fixed = TRUE), "proCAF", "restCAF"),
    duration  = as.numeric(d$TRACK_DURATION),                       # frames
    n_spots   = as.numeric(d$NUMBER_SPOTS),                         # frames in which the track was detected
    linearity = as.numeric(d$LINEARITY_OF_FORWARD_PROGRESSION),     # dimensionless
    dcr       = as.numeric(d$MEAN_DIRECTIONAL_CHANGE_RATE),         # rad / frame
    stringsAsFactors = FALSE)
}
tm <- do.call(rbind, pieces)
stopifnot(nrow(tm) == 22777, !anyNA(tm$duration), !anyNA(tm$n_spots), length(unique(tm$field)) == 55)
# units from the data: durations in frames (max 47 of 48), at most 48 spots, rates in rad/frame
stopifnot(max(tm$duration) == 47, max(tm$n_spots) == 48, max(tm$dcr, na.rm = TRUE) < pi, min(tm$dcr, na.rm = TRUE) >= 0)

# reproduce the labmate's plot (55 fields, TRACK_DURATION > 40, one-sided Pro > Rest):
# his figure says t p = 0.927, Wilcoxon p = 0.684
chk <- tm[tm$duration > MIN_DURATION, ]
chk$cafType <- factor(chk$cafType, levels = c("proCAF", "restCAF"))
p_chk_t <- t.test(linearity ~ cafType, data = chk, alternative = "greater")$p.value
p_chk_w <- wilcox.test(linearity ~ cafType, data = chk, alternative = "greater", exact = FALSE)$p.value
stopifnot(abs(p_chk_t - 0.927) < 0.001, abs(p_chk_w - 0.684) < 0.001)
cat(sprintf("STEP 1  read-check vs labmate's plot: t p = %.3f (his 0.927), w p = %.3f (his 0.684)  OK\n",
            p_chk_t, p_chk_w))

# thesis track set (>= 40 of 48 frames, 54 fields) must equal the Fig 4 set field by field
kept <- tm[tm$n_spots >= MIN_SPOTS & !tm$field %in% BAD_FIELDS, ]
gpf  <- read.csv(GPF_CSV, stringsAsFactors = FALSE)
n_kept <- table(kept$field)
stopifnot(length(n_kept) == 54, setequal(names(n_kept), gpf$field),
          all(as.integer(n_kept[gpf$field]) == gpf$n_tracks), nrow(kept) == 5606)
cat("        track filter NUMBER_SPOTS >= 40 reproduces Fig 4 n_tracks in 54/54 fields (5,606 tracks)\n")

meta <- as.data.frame(read_excel(META_XLSX, sheet = "Sheet1"))
meta$field <- sprintf("%s_field_%02d", meta$Tumor, as.integer(meta$Field))
stopifnot(!any(duplicated(meta$field)), all(kept$field %in% meta$field))
kept$donor     <- meta$cfECM[match(kept$field, meta$field)]
kept$tumorType <- meta$tumorType[match(kept$field, meta$field)]
# the Pro/Rest folder must agree with the metadata subtype for every field
stopifnot(all(kept$cafType == meta$cafType[match(kept$field, meta$field)]))
stopifnot(length(unique(kept$donor)) == 14, !anyNA(kept$tumorType))
cat(sprintf("        %s tracks in 54 fields (%d proCAF / %d restCAF fields), 14 donors\n",
            format(nrow(kept), big.mark = ","),
            length(unique(kept$field[kept$cafType == "proCAF"])),
            length(unique(kept$field[kept$cafType == "restCAF"]))))

# STEP 2. Field means -> donor means, and the tests at each level
# unit of replication = cfECM donor (mean over its fields), 6 proCAF vs 8 restCAF, Welch t as in Fig 4E
fld <- aggregate(cbind(linearity, dcr) ~ field + donor + cafType + tumorType, data = kept, FUN = mean)
fld$n_tracks <- as.integer(table(kept$field)[fld$field])
stopifnot(nrow(fld) == 54, all(fld$n_tracks >= 39))

dn <- aggregate(cbind(linearity, dcr) ~ donor + cafType, data = fld, FUN = mean)
dn$n_fields <- as.integer(table(fld$donor)[dn$donor])
dn$cafType  <- factor(dn$cafType, levels = c("proCAF", "restCAF"))
stopifnot(nrow(dn) == 14, all(table(dn$cafType) == c(6, 8)), all(dn$n_fields %in% 3:4))

tt_lin_dn <- t.test(linearity ~ cafType, data = dn)          # Welch, two-sided, donors
tt_dcr_dn <- t.test(dcr       ~ cafType, data = dn)
kept$cafType <- factor(kept$cafType, levels = c("proCAF", "restCAF"))
tt_lin_tr <- t.test(linearity ~ cafType, data = kept)        # Welch, two-sided, tracks (Supp Fig S3, comparison only)
tt_dcr_tr <- t.test(dcr       ~ cafType, data = kept)
n_tr      <- table(kept$cafType)                             # tracks per subtype

# the same, split by tumor subtype (2 basal-like lines, 2 classical lines per donor)
dn_t <- aggregate(cbind(linearity, dcr) ~ donor + cafType + tumorType, data = fld, FUN = mean)
dn_t$cafType   <- factor(dn_t$cafType, levels = c("proCAF", "restCAF"))
dn_t$tumorType <- factor(dn_t$tumorType, levels = c("basal", "classical"), labels = c("basal-like", "classical"))
stopifnot(nrow(dn_t) == 28)
tt_lin_bas <- t.test(linearity ~ cafType, data = dn_t[dn_t$tumorType == "basal-like", ])
tt_lin_cla <- t.test(linearity ~ cafType, data = dn_t[dn_t$tumorType == "classical",  ])
tt_dcr_bas <- t.test(dcr ~ cafType, data = dn_t[dn_t$tumorType == "basal-like", ])
tt_dcr_cla <- t.test(dcr ~ cafType, data = dn_t[dn_t$tumorType == "classical",  ])

cat(sprintf("STEP 2  linearity  donors: proCAF %.3f (%.3f-%.3f) vs restCAF %.3f (%.3f-%.3f), Welch p = %.3f | tracks p = %.3g\n",
            mean(dn$linearity[dn$cafType == "proCAF"]),  min(dn$linearity[dn$cafType == "proCAF"]),  max(dn$linearity[dn$cafType == "proCAF"]),
            mean(dn$linearity[dn$cafType == "restCAF"]), min(dn$linearity[dn$cafType == "restCAF"]), max(dn$linearity[dn$cafType == "restCAF"]),
            tt_lin_dn$p.value, tt_lin_tr$p.value))
cat(sprintf("        dir. change donors: proCAF %.3f (%.3f-%.3f) vs restCAF %.3f (%.3f-%.3f) rad/frame, Welch p = %.3f | tracks p = %.3g | random walk = %.4f\n",
            mean(dn$dcr[dn$cafType == "proCAF"]),  min(dn$dcr[dn$cafType == "proCAF"]),  max(dn$dcr[dn$cafType == "proCAF"]),
            mean(dn$dcr[dn$cafType == "restCAF"]), min(dn$dcr[dn$cafType == "restCAF"]), max(dn$dcr[dn$cafType == "restCAF"]),
            tt_dcr_dn$p.value, tt_dcr_tr$p.value, RANDOM_DCR))
cat(sprintf("        by tumor subtype (donor level): linearity basal p = %.3f, classical p = %.3f; dir. change basal p = %.3f, classical p = %.3f\n",
            tt_lin_bas$p.value, tt_lin_cla$p.value, tt_dcr_bas$p.value, tt_dcr_cla$p.value))

# STEP 3. ECM fiber spread per donor (OrientationPy on the single-plane confocal images)
# sd_around_0_deg: energy-weighted SD of fiber orientation around the dominant direction, deg (lower = more aligned)
ecm <- read.csv(ECM_CSV, stringsAsFactors = FALSE, colClasses = c(sample = "character"))
stopifnot(nrow(ecm) == 15, all(c("sd_around_0_deg", "perc_within_15deg", "group", "donor") %in% names(ecm)))
# two donors were never tracked, so their donor string is empty in this table; fill from the sample code
fill <- c("509" = "P240509N1", "0118N" = "P220118N1")
idx  <- is.na(ecm$donor) | ecm$donor == ""
stopifnot(sum(idx) == 2, all(ecm$sample[idx] %in% names(fill)))
ecm$donor[idx] <- fill[ecm$sample[idx]]
stopifnot(!any(duplicated(ecm$donor)), all(table(ecm$group) == c(7, 8)))
ecm <- ecm[, c("donor", "group", "sd_around_0_deg", "perc_within_15deg")]
names(ecm)[2] <- "cafType"
# de-identified alias used on every figure: P201015T1 -> P015T1, P220804N1 -> P804N1, P231207NA -> P207NA
ecm$alias <- paste0("P", substring(ecm$donor, 5))
stopifnot(all(nchar(ecm$alias) == 6), !any(duplicated(ecm$alias)))

# STEP 4. DeCAF scores, stiffness and marker intensities; join to the ECM table
decaf <- read.delim(DECAF_TSV, stringsAsFactors = FALSE, na.strings = "")
stopifnot(nrow(decaf) == 19, all(decaf$DeCAF == ifelse(decaf$DeCAF_score >= 0.5, "proCAF", "restCAF")),
          !any(duplicated(decaf$Line)))
if (!CP_LINES_IN) decaf <- decaf[decaf$Disease == "PDAC", ]
stopifnot(nrow(decaf) == 16)

# YM_MEAN / YM_MEDIAN: Young's modulus in Pa; marker columns (aSMA..PDGFRab) in arbitrary units.
# The sheet's own LINEARITY column (0.13-1.84, definition unknown) is not used.
stiff <- as.data.frame(read_excel(STIFF_XLSX, sheet = "Sheet1"))
stopifnot(nrow(stiff) == 14, identical(names(stiff)[1:4], c("Line", "LINEARITY", "YM_MEAN", "YM_MEDIAN")),
          !any(duplicated(stiff$Line)), all(stiff$YM_MEAN > 0), all(stiff$YM_MEDIAN > 0))
MARKERS <- names(stiff)[5:ncol(stiff)]
stopifnot(length(MARKERS) == 10, MARKERS[1] == "aSMA")
# every stiffness line must have a DeCAF call; the subtype comes from the DeCAF table
stopifnot(all(stiff$Line %in% decaf$Line))
stiff$cafType <- factor(decaf$DeCAF[match(stiff$Line, decaf$Line)], levels = c("proCAF", "restCAF"))
stiff$YM_MEAN_kPa   <- stiff$YM_MEAN   / 1000
stiff$YM_MEDIAN_kPa <- stiff$YM_MEDIAN / 1000

# 5b: ECM spread x DeCAF score; 13 of the 15 imaged donors have a score
b <- merge(ecm, decaf[, c("Line", "DeCAF_score", "DeCAF")], by.x = "donor", by.y = "Line")
stopifnot(nrow(b) == 13, all(b$cafType == b$DeCAF),
          setequal(setdiff(ecm$donor, b$donor), c("P220815T1", "P220118N1")))
b$cafType <- factor(b$cafType, levels = c("proCAF", "restCAF"))

# 5c: ECM spread x stiffness; 13 of the 14 stiffness lines were imaged
c3 <- merge(ecm, stiff[, c("Line", "YM_MEAN_kPa", "YM_MEDIAN_kPa")], by.x = "donor", by.y = "Line")
stopifnot(nrow(c3) == 13, setdiff(stiff$Line, c3$donor) == "P240815N1")
c3$cafType <- factor(c3$cafType, levels = c("proCAF", "restCAF"))

# 5f: markers x ECM spread; 13 lines have markers, 12 of them were imaged
mk <- stiff[complete.cases(stiff[, MARKERS]), c("Line", "cafType", MARKERS)]
stopifnot(nrow(mk) == 13, setdiff(stiff$Line, mk$Line) == "P220804T1")
f6 <- merge(ecm, mk[, c("Line", MARKERS)], by.x = "donor", by.y = "Line")
stopifnot(nrow(f6) == 12)

# STEP 5. Statistics for the scatter panels and the stiffness comparison
# exact = FALSE for Spearman: at n = 13 the exact p can underflow or fail on ties
ct_b_sd   <- cor.test(b$sd_around_0_deg,   b$DeCAF_score, method = "pearson")
sp_b_sd   <- cor.test(b$sd_around_0_deg,   b$DeCAF_score, method = "spearman", exact = FALSE)
ct_b_p15  <- cor.test(b$perc_within_15deg, b$DeCAF_score, method = "pearson")
sp_b_p15  <- cor.test(b$perc_within_15deg, b$DeCAF_score, method = "spearman", exact = FALSE)
r2_b_sd   <- summary(lm(DeCAF_score ~ sd_around_0_deg,   data = b))$r.squared
r2_b_p15  <- summary(lm(DeCAF_score ~ perc_within_15deg, data = b))$r.squared

ct_c_sd_med  <- cor.test(c3$sd_around_0_deg,   c3$YM_MEDIAN_kPa, method = "pearson")
sp_c_sd_med  <- cor.test(c3$sd_around_0_deg,   c3$YM_MEDIAN_kPa, method = "spearman", exact = FALSE)
ct_c_sd_mean <- cor.test(c3$sd_around_0_deg,   c3$YM_MEAN_kPa,   method = "pearson")
sp_c_sd_mean <- cor.test(c3$sd_around_0_deg,   c3$YM_MEAN_kPa,   method = "spearman", exact = FALSE)
ct_c_p15_med <- cor.test(c3$perc_within_15deg, c3$YM_MEDIAN_kPa, method = "pearson")
sp_c_p15_med <- cor.test(c3$perc_within_15deg, c3$YM_MEDIAN_kPa, method = "spearman", exact = FALSE)

tt_ym_med  <- t.test(YM_MEDIAN_kPa ~ cafType, data = stiff)   # 7 proCAF vs 7 restCAF lines
tt_ym_mean <- t.test(YM_MEAN_kPa   ~ cafType, data = stiff)

cat(sprintf("STEP 5  5b  DeCAF ~ SD(deg):  Pearson r = %.2f p = %.3g | Spearman rho = %.2f p = %.3g | R2 = %.2f  (n = 13)\n",
            ct_b_sd$estimate, ct_b_sd$p.value, sp_b_sd$estimate, sp_b_sd$p.value, r2_b_sd))
cat(sprintf("        5b  DeCAF ~ %%<15deg:  Pearson r = %.2f p = %.3g | Spearman rho = %.2f p = %.3g | R2 = %.2f\n",
            ct_b_p15$estimate, ct_b_p15$p.value, sp_b_p15$estimate, sp_b_p15$p.value, r2_b_p15))
cat(sprintf("        5c  YM median ~ SD(deg): Pearson r = %.2f p = %.3g | Spearman rho = %.2f p = %.3g (n = 13); YM mean ~ SD: r = %.2f p = %.3g\n",
            ct_c_sd_med$estimate, ct_c_sd_med$p.value, sp_c_sd_med$estimate, sp_c_sd_med$p.value,
            ct_c_sd_mean$estimate, ct_c_sd_mean$p.value))
cat(sprintf("        stiffness by subtype (kPa): median YM proCAF %.2f (%.2f-%.2f) vs restCAF %.2f (%.2f-%.2f), Welch p = %.3f; mean YM p = %.3f (7 vs 7 lines)\n",
            mean(stiff$YM_MEDIAN_kPa[stiff$cafType == "proCAF"]), min(stiff$YM_MEDIAN_kPa[stiff$cafType == "proCAF"]), max(stiff$YM_MEDIAN_kPa[stiff$cafType == "proCAF"]),
            mean(stiff$YM_MEDIAN_kPa[stiff$cafType == "restCAF"]), min(stiff$YM_MEDIAN_kPa[stiff$cafType == "restCAF"]), max(stiff$YM_MEDIAN_kPa[stiff$cafType == "restCAF"]),
            tt_ym_med$p.value, tt_ym_mean$p.value))

# p written as 0.95 at or above 0.1, 0.038 / 0.004 down to 0.001, 4.4e-07 below
p_all <- c(D = tt_lin_dn$p.value, E = tt_dcr_dn$p.value, Dt = tt_lin_tr$p.value, Et = tt_dcr_tr$p.value,
           B = ct_b_sd$p.value, B15 = ct_b_p15$p.value,
           C = ct_c_sd_med$p.value, Cm = ct_c_sd_mean$p.value, C15 = ct_c_p15_med$p.value,
           YM = tt_ym_med$p.value, YMm = tt_ym_mean$p.value,
           Db = tt_lin_bas$p.value, Dc = tt_lin_cla$p.value, Eb = tt_dcr_bas$p.value, Ec = tt_dcr_cla$p.value)
p_txt <- ifelse(p_all >= 0.1, sprintf("%.2f", p_all),
                ifelse(p_all >= 0.001, sprintf("%.3f", p_all), sprintf("%.1e", p_all)))
names(p_txt) <- names(p_all)

lab_D  <- sprintf("p = %s", p_txt["D"]);   lab_E  <- sprintf("p = %s", p_txt["E"])
lab_B  <- sprintf("r = %.2f, p = %s, R² = %.2f", ct_b_sd$estimate, p_txt["B"], r2_b_sd)
lab_C  <- sprintf("r = %.2f, p = %s", ct_c_sd_med$estimate, p_txt["C"])
lab_Cm <- sprintf("r = %.2f, p = %s", ct_c_sd_mean$estimate, p_txt["Cm"])
lab_Dt <- sprintf("p = %s", p_txt["Dt"]);  lab_Et <- sprintf("p = %s", p_txt["Et"])

# STEP 6. Theme, every size absolute (as Fig 4)
THEME <- theme_classic(base_size = PT_TITLE * FS, base_family = FONT) +
  theme(axis.line         = element_line(linewidth = PT_AXIS * LW, colour = "black", lineend = "square"),
        axis.ticks        = element_line(linewidth = PT_TICK * LW, colour = "grey20"),
        axis.ticks.length = unit(PT_TICK_LEN, "bigpts"),
        axis.text         = element_text(size = PT_TICKLAB * FS, colour = "black"),
        axis.text.x       = element_text(margin = margin(t = 1.5), lineheight = 0.95),
        axis.text.y       = element_text(margin = margin(r = 1.5)),
        axis.title        = element_text(size = PT_TITLE * FS, colour = "black"),
        axis.title.x      = element_text(margin = margin(t = 2.5)),
        axis.title.y      = element_text(margin = margin(r = 2.5)),
        plot.title        = element_text(size = PT_STAT * FS, hjust = 0.5, lineheight = 1.15,
                                         margin = margin(b = 2.5)),
        strip.background  = element_blank(),
        strip.text        = element_text(size = PT_STRIP * FS, colour = "black", margin = margin(b = 2)),
        panel.spacing     = unit(6, "bigpts"),
        legend.text       = element_text(size = PT_LEGEND * FS, margin = margin(l = 1.5)),
        legend.title      = element_text(size = PT_LEGEND * FS),
        legend.key.height = unit(PT_LEG_ROW, "bigpts"),
        legend.key.width  = unit(7, "bigpts"),
        legend.key        = element_blank(),
        legend.background = element_blank(),
        legend.margin     = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.spacing    = unit(0, "pt"),
        legend.position   = "inside",
        legend.position.inside = c(0.04, 0.99),
        legend.justification   = c(0, 1),
        plot.margin       = margin(3, 3, 2, 2),
        plot.background   = element_blank(),
        panel.background  = element_blank())

# STEP 7. Panels
# x title as in Fig 3H / 6E; at 7 pt it runs 1.5 pt past each end of the 99-pt x axis, inside the margin
X_SD  <- "SD around dominant direction (°)"
# y titles more than ~10 % longer than the 95-pt y axis go on two lines
Y_DECAF <- "DeCAF score\n(proCAF probability)"
Y_DCR   <- "Directional change rate\n(rad per frame)"

# 5b. DeCAF score against ECM spread, one point per donor (n = 13)
# dashed line = the 0.5 classification threshold; grey line = least-squares fit
pB <- ggplot(b, aes(sd_around_0_deg, DeCAF_score)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = PT_DASH * LW, colour = "grey35") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey45", linewidth = PT_FIT * LW) +
  geom_point(aes(colour = cafType, shape = cafType), size = PT_POINT * PS, stroke = 0) +
  scale_colour_manual(values = CAF_COLS) + scale_shape_manual(values = CAF_SHAPES) +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(x = X_SD, y = Y_DECAF, colour = NULL, shape = NULL, title = lab_B) +
  THEME + theme(legend.position.inside = c(0.98, 0.02), legend.justification = c(1, 0))
# 5c. Stiffness against ECM spread, one point per line (n = 13)
pC <- ggplot(c3, aes(sd_around_0_deg, YM_MEDIAN_kPa)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey45", linewidth = PT_FIT * LW) +
  geom_point(aes(colour = cafType, shape = cafType), size = PT_POINT * PS, stroke = 0) +
  scale_colour_manual(values = CAF_COLS) + scale_shape_manual(values = CAF_SHAPES) +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50)) +
  scale_y_continuous(limits = c(0, 3.5), breaks = c(0, 1, 2, 3)) +
  labs(x = X_SD, y = "Young's modulus, median (kPa)", colour = NULL, shape = NULL, title = lab_C) +
  THEME + theme(legend.position.inside = c(0.98, 0.98), legend.justification = c(1, 1))
# the same panel with the mean Young's modulus
pCm <- pC %+% c3 + aes(y = YM_MEAN_kPa) + labs(y = "Young's modulus, mean (kPa)", title = lab_Cm)
stopifnot(all(c3$YM_MEDIAN_kPa < 3.5), all(c3$YM_MEAN_kPa < 3.5))   # no point dropped by the 0-3.5 limits
# 5d / 5e. Migration metrics by CAF subtype, one point per cfECM donor (6 vs 8)
n_caf <- table(dn$cafType)
X_LAB <- c(proCAF = sprintf("proCAF\n(n = %d)", n_caf["proCAF"]), restCAF = sprintf("restCAF\n(n = %d)", n_caf["restCAF"]))
pD <- ggplot(dn, aes(cafType, linearity)) +
  geom_boxplot(aes(fill = cafType), colour = "black", alpha = 0.30, width = 0.5,
               linewidth = PT_BOX * LW, fatten = PT_MEDIAN / PT_BOX, staplewidth = 0.4, outlier.shape = NA) +
  geom_point(aes(colour = cafType, shape = cafType), size = PT_POINT * PS, stroke = 0,
             position = position_jitter(width = 0.08, height = 0, seed = 1)) +
  scale_fill_manual(values = CAF_COLS) + scale_colour_manual(values = CAF_COLS) + scale_shape_manual(values = CAF_SHAPES) +
  scale_x_discrete(labels = X_LAB) +
  scale_y_continuous(limits = c(0.15, 0.25), breaks = c(0.15, 0.20, 0.25)) +
  labs(x = NULL, y = "Linearity of forward progression", title = lab_D) +
  THEME + theme(legend.position = "none")
stopifnot(all(dn$linearity > 0.15 & dn$linearity < 0.25))       # nothing may be dropped by the limits

# dashed line = pi/2 rad/frame, the mean turning angle of an unbiased random walk
pE <- ggplot(dn, aes(cafType, dcr)) +
  geom_hline(yintercept = RANDOM_DCR, linetype = "dashed", linewidth = PT_DASH * LW, colour = "grey35") +
  geom_boxplot(aes(fill = cafType), colour = "black", alpha = 0.30, width = 0.5,
               linewidth = PT_BOX * LW, fatten = PT_MEDIAN / PT_BOX, staplewidth = 0.4, outlier.shape = NA) +
  geom_point(aes(colour = cafType, shape = cafType), size = PT_POINT * PS, stroke = 0,
             position = position_jitter(width = 0.08, height = 0, seed = 1)) +
  scale_fill_manual(values = CAF_COLS) + scale_colour_manual(values = CAF_COLS) + scale_shape_manual(values = CAF_SHAPES) +
  scale_x_discrete(labels = X_LAB) +
  scale_y_continuous(limits = c(1.45, 1.70), breaks = c(1.5, 1.6, 1.7)) +
  labs(x = NULL, y = Y_DCR, title = lab_E) +
  THEME + theme(legend.position = "none")
stopifnot(all(dn$dcr > 1.45 & dn$dcr < 1.70))

# Supp Fig S3: the same two metrics with every track as one observation (the labmate's unit), for
# comparison only: a field's tracks share one matrix and imaging session, so this p overstates the
# evidence. Boxes over all 5,606 tracks, outliers as faint dots.
X_TR <- c(proCAF  = sprintf("proCAF\n(%s tracks)",  format(n_tr["proCAF"],  big.mark = ",")),
          restCAF = sprintf("restCAF\n(%s tracks)", format(n_tr["restCAF"], big.mark = ",")))
stopifnot(!anyNA(kept$dcr), !anyNA(kept$linearity))
# the axis limits must not drop a single track (a dropped value would change the box statistics)
stopifnot(min(kept$dcr) >= 0, max(kept$dcr) <= 2.5, min(kept$linearity) >= 0, max(kept$linearity) <= 1.5)
pEt <- ggplot(kept, aes(cafType, dcr)) +
  geom_hline(yintercept = RANDOM_DCR, linetype = "dashed", linewidth = PT_DASH * LW, colour = "grey35") +
  geom_boxplot(aes(fill = cafType), colour = "black", alpha = 0.30, width = 0.5,
               linewidth = PT_BOX * LW, fatten = PT_MEDIAN / PT_BOX, staplewidth = 0.4,
               outlier.size = 0.6 * PS, outlier.stroke = 0, outlier.alpha = 0.25) +
  scale_fill_manual(values = CAF_COLS) +
  scale_x_discrete(labels = X_TR) +
  scale_y_continuous(limits = c(0, 2.5), breaks = c(0, 0.5, 1.0, 1.5, 2.0, 2.5)) +
  labs(x = NULL, y = Y_DCR, title = lab_Et) +
  THEME + theme(legend.position = "none")
pDt <- ggplot(kept, aes(cafType, linearity)) +
  geom_boxplot(aes(fill = cafType), colour = "black", alpha = 0.30, width = 0.5,
               linewidth = PT_BOX * LW, fatten = PT_MEDIAN / PT_BOX, staplewidth = 0.4,
               outlier.size = 0.6 * PS, outlier.stroke = 0, outlier.alpha = 0.25) +
  scale_fill_manual(values = CAF_COLS) +
  scale_x_discrete(labels = X_TR) +
  scale_y_continuous(limits = c(0, 1.5), breaks = c(0, 0.5, 1.0, 1.5)) +
  labs(x = NULL, y = "Linearity of forward progression", title = lab_Dt) +
  THEME + theme(legend.position = "none")

# 5f. The marker heatmap is drawn by figures/heatmaps_fig5_fig6.R (shared with Fig 6B). Its input:
# the 12 lines with markers and a fiber-orientation measurement, most -> least aligned, raw intensities (a.u.)
f6 <- f6[order(f6$sd_around_0_deg), ]
write.csv(f6[, c("donor", "alias", "cafType", "sd_around_0_deg", "perc_within_15deg", MARKERS)],
          file.path(OUT_DIR, "fig5F_marker_input.csv"), row.names = FALSE)

# STEP 8. Fix each plotting area to its slot size and write
panels <- list(B = pB, C = pC, Cm = pCm, D = pD, E = pE, S3A = pEt, S3B = pDt)
files  <- c(B  = "fig5B_decaf_vs_spread.pdf",  C = "fig5C_stiffness_median_vs_spread.pdf",
            Cm = "fig5C_stiffness_mean_vs_spread.pdf",
            D  = "fig5D_linearity_by_caf.pdf", E = "fig5E_dir_change_by_caf.pdf",
            S3A = "suppfig3A_dir_change_tracks.pdf", S3B = "suppfig3B_linearity_tracks.pdf")
dirs    <- c(B = OUT_DIR, C = OUT_DIR, Cm = OUT_DIR, D = OUT_DIR, E = OUT_DIR, S3A = SUPP3_DIR, S3B = SUPP3_DIR)
widths  <- c(B = PANEL_W_SCT, C = PANEL_W_SCT, Cm = PANEL_W_SCT, D = PANEL_W_BOX, E = PANEL_W_BOX,
             S3A = PANEL_W_BOX, S3B = PANEL_W_BOX)
heights <- c(B = PANEL_H, C = PANEL_H, Cm = PANEL_H, D = PANEL_H, E = PANEL_H, S3A = PANEL_H, S3B = PANEL_H)

sizes <- data.frame()
for (k in names(panels)) {
  # text is measured with a device open, so open a cairo device writing to nowhere before
  # ggplotGrob(): the default pdf() device does not know the font and leaves a stray Rplots.pdf
  cairo_pdf(nullfile(), family = FONT)
  g   <- ggplotGrob(panels[[k]])
  loc <- g$layout[grepl("^panel", g$layout$name), ]
  stopifnot(nrow(loc) >= 1)
  for (l in unique(loc$l)) g$widths[l]  <- unit(widths[[k]],  "bigpts")
  for (t in unique(loc$t)) g$heights[t] <- unit(heights[[k]], "bigpts")
  page_w <- sum(convertWidth(g$widths,   "in", valueOnly = TRUE))
  page_h <- sum(convertHeight(g$heights, "in", valueOnly = TRUE))
  invisible(dev.off())

  cairo_pdf(file.path(dirs[[k]], files[[k]]), width = page_w, height = page_h, family = FONT)
  grid.draw(g)
  invisible(dev.off())
  sizes <- rbind(sizes, data.frame(panel = k, file = files[[k]], n_panels = nrow(loc),
                                   plot_w_pt = widths[[k]] * nrow(loc), plot_h_pt = heights[[k]],
                                   page_w_pt = round(page_w * 72, 1), page_h_pt = round(page_h * 72, 1)))
  cat(sprintf("%-4s plotting area %5.1f x %5.1f pt x%d | page %5.1f x %5.1f pt  %s\n",
              k, widths[[k]], heights[[k]], nrow(loc), page_w * 72, page_h * 72, files[[k]]))
}
write.csv(sizes, file.path(OUT_DIR, "fig5_panel_sizes.csv"), row.names = FALSE)

# STEP 9. Every number for the legend, in one table
stats <- rbind(
  data.frame(panel = "5B", test = "Pearson, DeCAF score ~ SD around dominant direction (deg)", n = 13,
             estimate = unname(ct_b_sd$estimate), p = ct_b_sd$p.value, extra = sprintf("R2 = %.3f; Spearman rho = %.2f p = %.3g", r2_b_sd, sp_b_sd$estimate, sp_b_sd$p.value)),
  data.frame(panel = "5B_pct15", test = "Pearson, DeCAF score ~ % fibers within 15 deg", n = 13,
             estimate = unname(ct_b_p15$estimate), p = ct_b_p15$p.value, extra = sprintf("R2 = %.3f; Spearman rho = %.2f p = %.3g", r2_b_p15, sp_b_p15$estimate, sp_b_p15$p.value)),
  data.frame(panel = "5C", test = "Pearson, median Young's modulus (kPa) ~ SD (deg)", n = 13,
             estimate = unname(ct_c_sd_med$estimate), p = ct_c_sd_med$p.value, extra = sprintf("Spearman rho = %.2f p = %.3g", sp_c_sd_med$estimate, sp_c_sd_med$p.value)),
  data.frame(panel = "5C_mean", test = "Pearson, mean Young's modulus (kPa) ~ SD (deg)", n = 13,
             estimate = unname(ct_c_sd_mean$estimate), p = ct_c_sd_mean$p.value, extra = sprintf("Spearman rho = %.2f p = %.3g", sp_c_sd_mean$estimate, sp_c_sd_mean$p.value)),
  data.frame(panel = "5C_pct15", test = "Pearson, median Young's modulus (kPa) ~ % within 15 deg", n = 13,
             estimate = unname(ct_c_p15_med$estimate), p = ct_c_p15_med$p.value, extra = sprintf("Spearman rho = %.2f p = %.3g", sp_c_p15_med$estimate, sp_c_p15_med$p.value)),
  data.frame(panel = "extra_YM", test = "Welch t, median Young's modulus (kPa), proCAF vs restCAF lines", n = 14,
             estimate = diff(rev(tapply(stiff$YM_MEDIAN_kPa, stiff$cafType, mean))), p = tt_ym_med$p.value,
             extra = sprintf("proCAF mean %.2f (range %.2f-%.2f), restCAF mean %.2f (range %.2f-%.2f) kPa",
                             mean(stiff$YM_MEDIAN_kPa[stiff$cafType == "proCAF"]), min(stiff$YM_MEDIAN_kPa[stiff$cafType == "proCAF"]), max(stiff$YM_MEDIAN_kPa[stiff$cafType == "proCAF"]),
                             mean(stiff$YM_MEDIAN_kPa[stiff$cafType == "restCAF"]), min(stiff$YM_MEDIAN_kPa[stiff$cafType == "restCAF"]), max(stiff$YM_MEDIAN_kPa[stiff$cafType == "restCAF"]))),
  data.frame(panel = "extra_YM_mean", test = "Welch t, mean Young's modulus (kPa), proCAF vs restCAF lines", n = 14,
             estimate = diff(rev(tapply(stiff$YM_MEAN_kPa, stiff$cafType, mean))), p = tt_ym_mean$p.value,
             extra = sprintf("proCAF mean %.2f, restCAF mean %.2f kPa", mean(stiff$YM_MEAN_kPa[stiff$cafType == "proCAF"]), mean(stiff$YM_MEAN_kPa[stiff$cafType == "restCAF"]))),
  data.frame(panel = "5D", test = "Welch t, linearity of forward progression, donor means", n = 14,
             estimate = diff(rev(tapply(dn$linearity, dn$cafType, mean))), p = tt_lin_dn$p.value,
             extra = sprintf("proCAF %.3f (%.3f-%.3f, n = 6), restCAF %.3f (%.3f-%.3f, n = 8)",
                             mean(dn$linearity[dn$cafType == "proCAF"]), min(dn$linearity[dn$cafType == "proCAF"]), max(dn$linearity[dn$cafType == "proCAF"]),
                             mean(dn$linearity[dn$cafType == "restCAF"]), min(dn$linearity[dn$cafType == "restCAF"]), max(dn$linearity[dn$cafType == "restCAF"]))),
  data.frame(panel = "5E", test = "Welch t, mean directional change rate (rad/frame), donor means", n = 14,
             estimate = diff(rev(tapply(dn$dcr, dn$cafType, mean))), p = tt_dcr_dn$p.value,
             extra = sprintf("proCAF %.3f (%.3f-%.3f, n = 6), restCAF %.3f (%.3f-%.3f, n = 8); random walk = %.4f",
                             mean(dn$dcr[dn$cafType == "proCAF"]), min(dn$dcr[dn$cafType == "proCAF"]), max(dn$dcr[dn$cafType == "proCAF"]),
                             mean(dn$dcr[dn$cafType == "restCAF"]), min(dn$dcr[dn$cafType == "restCAF"]), max(dn$dcr[dn$cafType == "restCAF"]), RANDOM_DCR)),
  data.frame(panel = "5D_tracks", test = "Welch t, linearity, every track of the Fig 4 set (comparison only)", n = nrow(kept),
             estimate = diff(rev(tapply(kept$linearity, kept$cafType, mean))), p = tt_lin_tr$p.value,
             extra = sprintf("proCAF %.4f (n = %d), restCAF %.4f (n = %d)", mean(kept$linearity[kept$cafType == "proCAF"]), n_tr["proCAF"], mean(kept$linearity[kept$cafType == "restCAF"]), n_tr["restCAF"])),
  data.frame(panel = "5E_tracks", test = "Welch t, directional change rate, every track (comparison only)", n = nrow(kept),
             estimate = diff(rev(tapply(kept$dcr, kept$cafType, mean))), p = tt_dcr_tr$p.value,
             extra = sprintf("proCAF %.4f (n = %d), restCAF %.4f (n = %d)", mean(kept$dcr[kept$cafType == "proCAF"]), n_tr["proCAF"], mean(kept$dcr[kept$cafType == "restCAF"]), n_tr["restCAF"])),
  data.frame(panel = "extra_5D_by_tumor", test = "Welch t, linearity, donor means within tumor subtype", n = 14,
             estimate = NA, p = NA, extra = sprintf("basal-like p = %.3f, classical p = %.3f", tt_lin_bas$p.value, tt_lin_cla$p.value)),
  data.frame(panel = "extra_5E_by_tumor", test = "Welch t, directional change rate, donor means within tumor subtype", n = 14,
             estimate = NA, p = NA, extra = sprintf("basal-like p = %.3f, classical p = %.3f", tt_dcr_bas$p.value, tt_dcr_cla$p.value)),
  data.frame(panel = "5F", test = "heatmap drawn by 03_heatmaps_final_size.R from fig5F_marker_input.csv", n = 12,
             estimate = NA, p = NA, extra = paste("lines left->right:", paste(f6$alias, collapse = " "))))
write.csv(stats, file.path(OUT_DIR, "fig5_stats.csv"), row.names = FALSE)
write.csv(dn,  file.path(OUT_DIR, "fig5DE_donor_means.csv"), row.names = FALSE)
write.csv(fld, file.path(OUT_DIR, "fig5DE_field_means.csv"), row.names = FALSE)
write.csv(b[, c("donor", "alias", "cafType", "sd_around_0_deg", "perc_within_15deg", "DeCAF_score")],
          file.path(OUT_DIR, "fig5B_points.csv"), row.names = FALSE)
write.csv(c3[, c("donor", "alias", "cafType", "sd_around_0_deg", "perc_within_15deg", "YM_MEAN_kPa", "YM_MEDIAN_kPa")],
          file.path(OUT_DIR, "fig5C_points.csv"), row.names = FALSE)
cat("written to", OUT_DIR, "\n")
