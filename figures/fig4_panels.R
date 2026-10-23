# Fig 4E-H at final size: guidance G per cfECM donor by CAF subtype (E), against local fiber
# order R_local (F), on basal-like vs classical lines (G), and against R_local per tumor subtype (H).
# G is the TrackMate G from guidance_per_field.csv; the rho / p text goes above each panel.
# Out: results/guidance/subtype/final_size/fig4{E,F,G,H}_*.pdf

library(ggplot2)
library(grid)

# USER INPUTS
ROOT    <- "."
IN_CSV  <- file.path(ROOT, "results/guidance/alignment_scales.csv")    # matrix: R_local, donor, subtype
IN_GPF  <- file.path(ROOT, "results/guidance/guidance_per_field.csv")  # G from guidance/03_guidance.py (TrackMate)
OUT_DIR <- file.path(ROOT, "results/guidance/subtype/final_size")

# panels are drawn at final size; every size below is in PDF pt (1/72 in)
PANEL_H     <- 95                                   # y axis length, the same for all four
PANEL_W     <- c(E = 98, F = 99, G = 95, H = 102.5) # x axis length; G is square
PT_AXIS     <- 0.59    # axis line width
PT_TICK     <- 0.36    # tick mark width
PT_TICK_LEN <- 0.95    # tick mark length
PT_TICKLAB  <- 5.4     # tick label font size
PT_TITLE    <- 5.8     # axis title font size
PT_LEGEND   <- 4.9     # legend font size
PT_LEG_ROW  <- 6.8     # legend row pitch
PT_STAT     <- 5.6     # rho / p annotation font size
PT_POINT    <- 2.05    # point diameter (circle) / side (square)
PT_FIT      <- 0.59    # regression line width (panel H)
PT_LINK     <- 0.37    # grey basal-classical connector width (panel H)
PT_BOX      <- 0.345   # box and whisker line width (panel E)
PT_MEDIAN   <- 0.78    # median line width (panel E)
PT_DASH     <- 0.37    # equality line width (panel G)

# one common G scale in all four panels
G_LIMITS <- c(0.04, 0.31)     # donor-level G spans 0.053 to 0.292 (dimensionless, -1 to +1)
G_BREAKS <- c(0.1, 0.2, 0.3)

PRO       <- "#E74095"
REST      <- "#288D92"
BASAL     <- "#F89B0E"
CLASSICAL <- "#3C4DA8"
FONT      <- "Arial"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# unit conversions from PDF pt: 1 R lwd = 1/96 in = 0.75 pt; ggplot font sizes are in 1/72.27 in;
# pch 15 / 16 are 0.75 * fontsize across, with fontsize = size * .pt when stroke = 0
LW <- 96 / (72 * .pt)
FS <- 72.27 / 72
PS <- FS / (0.75 * .pt)

# STEP 1. Data
d   <- read.csv(IN_CSV, stringsAsFactors = FALSE)
gpf <- read.csv(IN_GPF, stringsAsFactors = FALSE)

# alignment_scales.csv still carries the old trackpy G; swap in the TrackMate G field by field
stopifnot(!any(duplicated(d$field)), !any(duplicated(gpf$field)),
          setequal(d$field, gpf$field), !"P170713_field_13" %in% gpf$field)
d$G <- gpf$G[match(d$field, gpf$field)]
stopifnot(!anyNA(d$G), nrow(d) == 54, sum(gpf$n_steps) > sum(d$n_steps))

d$tumorType <- factor(d$tumorType, levels = c("basal", "classical"))
d$cafType   <- factor(ifelse(d$cafType == "pro", "proCAF", "restCAF"))

# one value per cfECM donor: the mean over that donor's 3-4 fields
dn <- aggregate(cbind(G, R_local) ~ donor, data = d, FUN = mean)
dn <- merge(dn, unique(d[, c("donor", "cafType")]), by = "donor")

# the same, separately for the basal-like and the classical tumor lines
dn_g <- aggregate(G ~ donor + tumorType, data = d, FUN = mean)
bas  <- dn_g[dn_g$tumorType == "basal",     c("donor", "G")]; names(bas)[2] <- "G_basal"
cla  <- dn_g[dn_g$tumorType == "classical", c("donor", "G")]; names(cla)[2] <- "G_classical"
dn2  <- merge(merge(dn, bas, by = "donor"), cla, by = "donor")
stopifnot(nrow(dn) == 14, nrow(dn2) == 14)

# ggplot would silently drop anything outside the common G scale
stopifnot(all(dn$G > G_LIMITS[1] & dn$G < G_LIMITS[2]),
          all(dn2$G_basal > G_LIMITS[1] & dn2$G_basal < G_LIMITS[2]),
          all(dn2$G_classical > G_LIMITS[1] & dn2$G_classical < G_LIMITS[2]))

# STEP 2. Statistics and the text above each panel
# exact = FALSE: at n = 14 R's exact Spearman p underflows to 0
tt  <- t.test(G ~ cafType, data = dn)                                             # E
ct  <- cor.test(dn$R_local, dn$G, method = "spearman", exact = FALSE)             # F
rbc <- cor.test(dn2$G_basal, dn2$G_classical, method = "spearman", exact = FALSE) # G
ctb <- cor.test(dn2$R_local, dn2$G_basal,     method = "spearman", exact = FALSE) # H
ctc <- cor.test(dn2$R_local, dn2$G_classical, method = "spearman", exact = FALSE) # H

# p written as 0.004 at or above 0.001, as 4.4e-07 below
p_all <- c(E = tt$p.value, F = ct$p.value, G = rbc$p.value, Hb = ctb$p.value, Hc = ctc$p.value)
p_txt <- ifelse(p_all >= 0.001, sprintf("%.3f", p_all), sprintf("%.1e", p_all))
if (p_all["E"] >= 0.01) p_txt["E"] <- sprintf("%.2f", p_all["E"])

RHO <- "ρ"
lab_E <- sprintf("p = %s", p_txt["E"])
lab_F <- sprintf("%s = %.2f, p = %s", RHO, ct$estimate,  p_txt["F"])
lab_G <- sprintf("%s = %.2f, p = %s", RHO, rbc$estimate, p_txt["G"])
lab_H <- sprintf("Basal-like: %s = %.2f, p = %s\nClassical: %s = %.2f, p = %s",
                 RHO, ctb$estimate, p_txt["Hb"], RHO, ctc$estimate, p_txt["Hc"])
cat("E:", lab_E, "\nF:", lab_F, "\nG:", lab_G, "\nH:", gsub("\n", " | ", lab_H), "\n")

# STEP 3. Theme, every size absolute
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
        legend.text       = element_text(size = PT_LEGEND * FS, margin = margin(l = 1.5)),
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

# STEP 4. Panels
# E. G by CAF subtype, 6 proCAF vs 8 restCAF donors
n_caf <- table(dn$cafType)
pE <- ggplot(dn, aes(cafType, G)) +
  geom_boxplot(aes(fill = cafType), colour = "black", alpha = 0.30, width = 0.5,
               linewidth = PT_BOX * LW, fatten = PT_MEDIAN / PT_BOX,
               staplewidth = 0.4, outlier.shape = NA) +
  geom_point(aes(colour = cafType, shape = cafType), size = PT_POINT * PS, stroke = 0,
             position = position_jitter(width = 0.08, height = 0, seed = 1)) +
  scale_fill_manual(values = c(proCAF = PRO, restCAF = REST)) +
  scale_colour_manual(values = c(proCAF = PRO, restCAF = REST)) +
  scale_shape_manual(values = c(proCAF = 16, restCAF = 15)) +
  scale_x_discrete(labels = c(proCAF  = sprintf("proCAF\n(n = %d)",  n_caf["proCAF"]),
                              restCAF = sprintf("restCAF\n(n = %d)", n_caf["restCAF"]))) +
  scale_y_continuous(limits = G_LIMITS, breaks = G_BREAKS) +
  labs(x = NULL, y = "G (donor mean)", title = lab_E) +
  THEME + theme(legend.position = "none")

# F. G against local fiber order, one point per donor
pF <- ggplot(dn, aes(R_local, G, colour = cafType, shape = cafType)) +
  geom_point(size = PT_POINT * PS, stroke = 0) +
  scale_colour_manual(values = c(proCAF = PRO, restCAF = REST)) +
  scale_shape_manual(values = c(proCAF = 16, restCAF = 15)) +
  scale_x_continuous(breaks = c(0.91, 0.93, 0.95, 0.97)) +
  scale_y_continuous(limits = G_LIMITS, breaks = G_BREAKS) +
  labs(x = "Local fiber order R_local", y = "G (donor mean)", colour = NULL, shape = NULL,
       title = lab_F) +
  THEME

# G. the same donor matrix under basal-like and under classical lines; equal axis limits and a
# square panel (STEP 5), so the dashed line is the true line of equality
pG <- ggplot(dn2, aes(G_basal, G_classical, colour = cafType, shape = cafType)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = PT_DASH * LW,
              colour = "grey35") +
  geom_point(size = PT_POINT * PS, stroke = 0) +
  scale_colour_manual(values = c(proCAF = PRO, restCAF = REST)) +
  scale_shape_manual(values = c(proCAF = 16, restCAF = 15)) +
  scale_x_continuous(limits = G_LIMITS, breaks = G_BREAKS) +
  scale_y_continuous(limits = G_LIMITS, breaks = G_BREAKS) +
  labs(x = "G on basal-like lines (donor mean)", y = "G on classical lines (donor mean)",
       colour = NULL, shape = NULL, title = lab_G) +
  THEME

# H. G against donor R_local, separately for the two tumor subtypes
long_r <- rbind(data.frame(donor = dn2$donor, R_local = dn2$R_local, G = dn2$G_basal,     tumorType = "basal-like"),
                data.frame(donor = dn2$donor, R_local = dn2$R_local, G = dn2$G_classical, tumorType = "classical"))
pH <- ggplot(long_r, aes(R_local, G)) +
  geom_line(aes(group = donor), colour = "grey80", linewidth = PT_LINK * LW) +
  geom_smooth(aes(colour = tumorType), method = "lm", formula = y ~ x, se = FALSE,
              linewidth = PT_FIT * LW) +
  geom_point(aes(colour = tumorType), shape = 16, size = PT_POINT * PS, stroke = 0) +
  scale_colour_manual(values = c("basal-like" = BASAL, classical = CLASSICAL)) +
  scale_x_continuous(breaks = c(0.91, 0.93, 0.95, 0.97)) +
  scale_y_continuous(limits = G_LIMITS, breaks = G_BREAKS) +
  labs(x = "Local fiber order R_local", y = "G (donor mean)", colour = NULL, title = lab_H) +
  THEME

# STEP 5. Fix the plotting area to its slot size and write each panel
# the plotting area is fixed in pt and the page is made just large enough for it plus the axis text
panels <- list(E = pE, F = pF, G = pG, H = pH)
files  <- c(E = "fig4E_G_by_caf.pdf", F = "fig4F_G_vs_Rlocal.pdf",
            G = "fig4G_basal_vs_classical.pdf", H = "fig4H_G_vs_Rlocal_by_subtype.pdf")

for (k in names(panels)) {
  # ggplotGrob() measures text, so open a cairo device to nowhere first; the default pdf()
  # device does not know Arial (wrong text metrics) and leaves a stray Rplots.pdf
  cairo_pdf(nullfile(), family = FONT)
  g   <- ggplotGrob(panels[[k]])
  loc <- g$layout[g$layout$name == "panel", ]
  stopifnot(nrow(loc) == 1)
  g$widths[loc$l]  <- unit(PANEL_W[[k]], "bigpts")
  g$heights[loc$t] <- unit(PANEL_H,      "bigpts")
  page_w <- sum(convertWidth(g$widths,   "in", valueOnly = TRUE))
  page_h <- sum(convertHeight(g$heights, "in", valueOnly = TRUE))
  invisible(dev.off())

  cairo_pdf(file.path(OUT_DIR, files[[k]]), width = page_w, height = page_h, family = FONT)
  grid.draw(g)
  invisible(dev.off())
  cat(sprintf("%s  plotting area %5.1f x %4.1f pt | page %5.1f x %5.1f pt (%.3f x %.3f in)  %s\n",
              k, PANEL_W[[k]], PANEL_H, page_w * 72, page_h * 72, page_w, page_h, files[[k]]))
}
cat("written to", OUT_DIR, "\n")
