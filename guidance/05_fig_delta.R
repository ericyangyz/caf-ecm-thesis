# Is delta (acute angle between a step and the fiber under it, 0-90 deg, 45 = unrelated) below chance?
# Tested against each field's own permutation null, not 45: field-wide drift of cells and fibers
# alone pushes delta below 45. Reads results/guidance/; writes fig_delta.pdf, fig_guidance_score.pdf.
library(ggplot2)

ROOT <- "./results/guidance"
PRO  <- "#E74095"
REST <- "#288D92"
REAL <- "#1a1a1a"
NUL  <- "#b8bcc0"
DOT  <- "#E70495"          # proCAF dots in b
REAL_FILL <- "#D9D2E9"     # real guidance score distribution

pool <- read.csv(file.path(ROOT, "delta_pooled.csv"))
fld  <- read.csv(file.path(ROOT, "delta_per_field.csv"))
st   <- read.csv(file.path(ROOT, "delta_stats.csv"))
gsc  <- read.csv(file.path(ROOT, "guidance_score_pooled.csv"))
fld$cafType <- factor(ifelse(fld$cafType == "pro", "proCAF", "restCAF"))

# a: pooled delta distribution, real vs shuffled
long <- rbind(
  data.frame(mid = pool$mid, v = pool$real,     grp = "real: the fiber under the cell"),
  data.frame(mid = pool$mid, v = pool$shuffled, grp = "shuffled: another cell's fiber"))
long$grp <- factor(long$grp, levels = c("real: the fiber under the cell",
                                        "shuffled: another cell's fiber"))

pa <- ggplot(long, aes(mid, v, color = grp)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.4, color = "grey35") +
  geom_vline(xintercept = 45, linetype = "dotted", linewidth = 0.4, color = "grey55") +
  geom_step(linewidth = 0.8) +
  annotate("text", x = 88, y = 1.03, label = "1.0 = no relationship",
           hjust = 1, vjust = 0, size = 3, color = "grey35") +
  annotate("text", x = 44.2, y = 0.79, label = "45° = chance",
           hjust = 0, angle = 90, size = 3, color = "grey45") +
  scale_color_manual(values = c(REAL, NUL)) +
  scale_x_continuous(breaks = seq(0, 90, 15), limits = c(0, 90), expand = c(0.01, 0)) +
  labs(x = "δ, angle between the step and the fiber beneath it (°)",
       y = "Relative frequency (1.0 = flat)", color = NULL,
       title = "Steps lean toward the fiber",
       subtitle = sprintf("n = %s steps · mean δ = %.1f° vs %.1f° shuffled",
                          format(st$n_steps, big.mark = ","), st$mean_real, st$mean_null)) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey30"),
        legend.position = "inside", legend.position.inside = c(0.98, 0.97),
        legend.justification = c(1, 1), legend.background = element_blank(),
        legend.key.height = unit(11, "pt"), legend.text = element_text(size = 8.5))

# b: each field's mean delta vs the 95% range of its 200 shuffles, ranked
fld <- fld[order(fld$delta_real), ]
fld$rank <- seq_len(nrow(fld))

pb <- ggplot(fld) +
  geom_vline(xintercept = 45, linetype = "dotted", linewidth = 0.4, color = "black") +
  geom_segment(aes(x = null_lo, xend = null_hi, y = rank, yend = rank,
                   linetype = "shuffled, 95% of 200 permutations"),
               color = NUL, linewidth = 2.1) +
  geom_point(aes(delta_real, rank, color = cafType, shape = cafType), size = 2.1) +
  annotate("text", x = 44.4, y = nrow(fld) * 0.30, label = "45° = chance",
           hjust = 0, angle = 90, size = 3, color = "black") +
  scale_color_manual(values = c(proCAF = DOT, restCAF = REST)) +
  scale_shape_manual(values = c(proCAF = 16, restCAF = 15)) +
  scale_linetype_manual(values = c("solid"), name = NULL) +
  guides(linetype = guide_legend(override.aes = list(color = NUL, linewidth = 2.1))) +
  labs(x = "Mean δ per field (°)", y = sprintf("the %d fields, ranked", nrow(fld)),
       color = NULL, shape = NULL,
       title = sprintf("%d of %d fields beat their own shuffle", st$n_sig, st$n_fields),
       subtitle = sprintf("%d / %d donors below their null · Wilcoxon p = %.4f",
                          st$n_donor_below, st$n_donors, st$p_donor)) +
  theme_classic(base_size = 11, base_family = "Arial") +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        plot.subtitle = element_text(hjust = 0.5, size = 8.5, color = "grey30"),
        axis.title = element_text(size = 10), axis.text = element_text(size = 9),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "inside", legend.position.inside = c(0.02, 0.98),
        legend.justification = c(0, 1), legend.background = element_blank(),
        legend.spacing.y = unit(0, "pt"), legend.key.height = unit(12, "pt"),
        legend.text = element_text(size = 8.5))

# c: guidance score A = mean cos 2*delta per field, real vs shuffled; G = A_real - A_null.
# per field, not per step: a per-step cos 2*delta histogram is U-shaped whatever the truth
# (short noisy steps pile up at both ends)
gpf <- read.csv(file.path(ROOT, "guidance_per_field.csv"))
gl  <- rbind(
  data.frame(v = gpf$A_real, grp = "real: the fiber under the cell"),
  data.frame(v = gpf$A_null, grp = "shuffled: another cell's fiber"))
gl$grp <- factor(gl$grp, levels = levels(long$grp))
m_real <- mean(gpf$A_real); m_null <- mean(gpf$A_null)

pc <- ggplot(gl, aes(v, fill = grp, color = grp)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey35") +
  geom_density(alpha = 0.6, linewidth = 0.7, bw = 0.025) +
  geom_rug(aes(color = grp), linewidth = 0.4, alpha = 0.8, show.legend = FALSE) +
  annotate("segment", x = m_null, xend = m_real, y = 11.5, yend = 11.5,
           arrow = arrow(length = unit(4, "pt"), ends = "both"), linewidth = 0.4) +
  annotate("text", x = (m_real + m_null) / 2, y = 12.6,
           label = sprintf("mean G = %.3f", m_real - m_null), size = 3.2,
           fontface = "bold") +
  annotate("text", x = 0.012, y = 6, label = "0 = no relationship",
           hjust = 0, size = 3, color = "grey35") +
  scale_fill_manual(values = c(REAL_FILL, NUL)) +
  scale_color_manual(values = c(REAL, NUL)) +
  scale_x_continuous(limits = c(-0.2, 0.55), breaks = seq(-0.2, 0.5, 0.1)) +
  labs(x = "Axial order parameter  A = mean cos 2δ  per field",
       y = "Density", fill = NULL, color = NULL,
       title = "Guidance score, one value per field",
       subtitle = sprintf("%d fields · mean %.3f real vs %.3f shuffled",
                          nrow(gpf), m_real, m_null)) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey30"),
        legend.position = "inside", legend.position.inside = c(0.98, 0.97),
        legend.justification = c(1, 1), legend.background = element_blank(),
        legend.text = element_text(size = 8.5))

ggsave(file.path(ROOT, "fig_guidance_score.pdf"), pc, width = 6, height = 4.4,
       device = cairo_pdf)

# panel a (pa) is not saved; it is not part of the thesis figures
ggsave(file.path(ROOT, "fig_delta.pdf"), pb, width = 4.8, height = 4.8, device = cairo_pdf)
cat("written\n")
