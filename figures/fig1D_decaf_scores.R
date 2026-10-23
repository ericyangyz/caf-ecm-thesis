# Thesis Fig 1D: DeCAF score of each thesis cfECM donor (lollipop colored by subtype call,
# dashed line at 0.5), drawn at final size to fit the 195 x 127 pt Fig 1D slot.
# In: data/data_decaf/decaf_scores_20260922.tsv
# Out: results/thesis/fig1/fig1D_decaf_scores_<n>donors.pdf (+ .csv of what was drawn)

suppressPackageStartupMessages({ library(ggplot2); library(grid) })

# USER INPUTS
ROOT      <- "."
DECAF_TSV <- file.path(ROOT, "data/data_decaf/decaf_scores_20260922.tsv")
OUT_DIR   <- file.path(ROOT, "results/thesis/fig1")

# the 17 thesis donors (Supp Table S1: 7 proCAF, 10 restCAF) and their S1 call
THESIS_DONORS <- c(P201015T1 = "proCAF", P220804N1 = "proCAF", P220804T1 = "proCAF",
                   P220815T1 = "proCAF", P220908N1 = "proCAF", P231207T1 = "proCAF",
                   P240229T1 = "proCAF",
                   P140807T1 = "restCAF", P210923T1 = "restCAF", P220201N1 = "restCAF",
                   P220212T1 = "restCAF", P220519T1 = "restCAF", P221108T1 = "restCAF",
                   P231207NA = "restCAF", P240509N1 = "restCAF", P240711N1 = "restCAF",
                   P220118N1 = "restCAF")
# P220815T1 and P220118N1 are missing from the DeCAF table; add as c(P220815T1 = 0.xx, ...)
EXTRA_SCORES <- c()
# donors ringed in Fig 1D (the two representative lines); c() for none
HIGHLIGHT <- c("P201015T1", "P220519T1")

# sizes in pt at final print size (same as Fig 4 / Fig 5)
ROW_PITCH   <- 6.6     # pt per donor row -> 17 rows = 112 pt; 15 rows = 99 pt
PANEL_W     <- 150     # x axis length, pt
PT_AXIS     <- 0.59;  PT_TICK <- 0.36;  PT_TICK_LEN <- 0.95
PT_TICKLAB  <- 5.4;   PT_TITLE <- 5.8;  PT_LEGEND <- 4.9;  PT_LEG_ROW <- 6.8
PT_POINT    <- 2.6     # lollipop head
PT_STEM     <- 0.37    # grey stem
PT_DASH     <- 0.37    # 0.5 threshold line
PT_RING     <- 0.5     # highlight ring line width
RING_DIAM   <- 5.0

PRO  <- "#E74095";  REST <- "#288D92"    # proCAF / restCAF colors
CAF_COLS <- c(proCAF = PRO, restCAF = REST)
FONT <- "Arial"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# pt -> ggplot units for line width (LW), font size (FS) and point size (PS)
LW <- 96 / (72 * .pt);  FS <- 72.27 / 72;  PS <- FS / (0.75 * .pt)
# shape-21 border: ggplot draws it with lwd = stroke * .stroke / 2, and 1 lwd = 1/96 in
RING_STROKE <- (PT_RING * 96 / 72) * 2 / .stroke

# STEP 1. Scores of the thesis donors
decaf <- read.delim(DECAF_TSV, stringsAsFactors = FALSE, na.strings = "")
stopifnot(!any(duplicated(decaf$Line)), all(decaf$DeCAF == ifelse(decaf$DeCAF_score >= 0.5, "proCAF", "restCAF")))
sc <- setNames(decaf$DeCAF_score, decaf$Line)
if (length(EXTRA_SCORES)) {
  stopifnot(!any(names(EXTRA_SCORES) %in% names(sc)), all(names(EXTRA_SCORES) %in% names(THESIS_DONORS)))
  sc <- c(sc, EXTRA_SCORES)
}
have    <- names(THESIS_DONORS)[names(THESIS_DONORS) %in% names(sc)]
missing <- setdiff(names(THESIS_DONORS), have)
d <- data.frame(donor = have, score = unname(sc[have]), call_S1 = unname(THESIS_DONORS[have]),
                stringsAsFactors = FALSE)
d$call <- ifelse(d$score >= 0.5, "proCAF", "restCAF")
# the call from the score must match the S1 call
stopifnot(all(d$call == d$call_S1))
d$alias <- paste0("P", substring(d$donor, 5))          # P201015T1 -> P015T1, as in Fig 4 / Fig 6
stopifnot(all(nchar(d$alias) == 6), !any(duplicated(d$alias)))
d <- d[order(d$score), ]                                # bottom -> top = low -> high
d$alias <- factor(d$alias, levels = d$alias)
d$call  <- factor(d$call, levels = c("proCAF", "restCAF"))
cat(sprintf("STEP 1  %d of %d thesis donors have a score (%d proCAF, %d restCAF); missing: %s\n",
            nrow(d), length(THESIS_DONORS), sum(d$call == "proCAF"), sum(d$call == "restCAF"),
            if (length(missing)) paste(missing, collapse = ", ") else "none"))
stopifnot(all(HIGHLIGHT %in% names(THESIS_DONORS)))
ring <- d[d$donor %in% HIGHLIGHT, ]

# STEP 2. The panel
n_pro <- sum(d$call == "proCAF"); n_rest <- sum(d$call == "restCAF")
p <- ggplot(d, aes(score, alias)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = PT_DASH * LW, colour = "grey35") +
  geom_segment(aes(x = 0, xend = score, yend = alias), colour = "grey70", linewidth = PT_STEM * LW) +
  geom_point(aes(colour = call), size = PT_POINT * PS, stroke = 0) +
  geom_point(data = ring, shape = 21, size = RING_DIAM * PS, stroke = RING_STROKE,
             colour = "black", fill = NA) +
  scale_colour_manual(values = CAF_COLS,
                      labels = c(proCAF = sprintf("proCAF (n = %d)", n_pro), restCAF = sprintf("restCAF (n = %d)", n_rest))) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"),
                     expand = expansion(mult = c(0, 0.03))) +
  scale_y_discrete(expand = expansion(add = 0.6)) +
  labs(x = "DeCAF score (proCAF probability)", y = NULL, colour = NULL) +
  theme_classic(base_size = PT_TITLE * FS, base_family = FONT) +
  theme(axis.line         = element_line(linewidth = PT_AXIS * LW, colour = "black", lineend = "square"),
        axis.ticks        = element_line(linewidth = PT_TICK * LW, colour = "grey20"),
        axis.ticks.length = unit(PT_TICK_LEN, "bigpts"),
        axis.text         = element_text(size = PT_TICKLAB * FS, colour = "black"),
        axis.text.x       = element_text(margin = margin(t = 1.5)),
        axis.text.y       = element_text(margin = margin(r = 1.5)),
        axis.title        = element_text(size = PT_TITLE * FS, colour = "black"),
        axis.title.x      = element_text(margin = margin(t = 2.5)),
        legend.text       = element_text(size = PT_LEGEND * FS, margin = margin(l = 1.5)),
        legend.key.height = unit(PT_LEG_ROW, "bigpts"),
        legend.key.width  = unit(7, "bigpts"),
        legend.key        = element_blank(),
        legend.background = element_blank(),
        legend.margin     = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.position   = "inside",
        legend.position.inside = c(0.98, 0.35),
        legend.justification   = c(1, 0.5),
        plot.margin       = margin(3, 3, 2, 2),
        plot.background   = element_blank(),
        panel.background  = element_blank())

# STEP 3. Fix the plotting area to (n rows x ROW_PITCH) by PANEL_W pt and write at 100%
panel_h <- nrow(d) * ROW_PITCH
cairo_pdf(nullfile(), family = FONT)
g   <- ggplotGrob(p)
loc <- g$layout[g$layout$name == "panel", ]
stopifnot(nrow(loc) == 1)
g$widths[loc$l]  <- unit(PANEL_W, "bigpts")
g$heights[loc$t] <- unit(panel_h, "bigpts")
page_w <- sum(convertWidth(g$widths,   "in", valueOnly = TRUE))
page_h <- sum(convertHeight(g$heights, "in", valueOnly = TRUE))
invisible(dev.off())
out <- file.path(OUT_DIR, sprintf("fig1D_decaf_scores_%ddonors.pdf", nrow(d)))
cairo_pdf(out, width = page_w, height = page_h, family = FONT)
grid.draw(g)
invisible(dev.off())
write.csv(d[, c("donor", "alias", "score", "call")], sub("\\.pdf$", ".csv", out), row.names = FALSE)
cat(sprintf("plotting area %5.1f x %5.1f pt | page %5.1f x %5.1f pt (slot in Fig 1 is 195 x 127 pt)  %s\n",
            PANEL_W, panel_h, page_w * 72, page_h * 72, out))
