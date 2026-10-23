# UMA vs Mathematica alignment, same 15 old Keyence images -> results/fig13_uma_vs_mathematica.pdf/.png
# UMA         = sd_around_0_deg, |grad I|^2-weighted, deg
# Mathematica = Variance[oCentered], unweighted, rad^2 -> deg here
library(ggplot2)
library(ggrepel)

ROOT <- "."

# ---- Mathematica values, in the notebook's proSamples / restSamples order ----
m <- data.frame(
  sample = c("015", "229T", "804N", "804T", "815", "908", "1207T",
             "0118N", "201N", "212", "509", "519", "711", "807", "1108"),
  var_rad2 = c(0.657998, 0.618538, 0.51705, 0.772362, 0.407781, 0.583037, 0.589696,
               0.66782, 0.465514, 0.457228, 0.42473, 0.391601, 0.466767, 0.452664, 0.430305),
  # de-identified alias, for the figure
  alias = c("P015T1", "P229T1", "P804N1", "P804T1", "P815T1", "P908N1", "P207T1",
            "P118N1", "P201N1", "P212T1", "P509N1", "P519T1", "P711N1", "P807T1", "P108T1"))

m$math_sd <- sqrt(m$var_rad2) * 180 / pi        # rad^2 -> deg

# ---- join to UMA ----
uma <- read.csv(file.path(ROOT, "results/old_keyence_alignment_per_sample.csv"))
d <- merge(uma[, c("group", "sample", "sd_around_0_deg")], m, by = "sample")
stopifnot(nrow(d) == 15)

# ---- correlation ----
ct <- cor.test(d$sd_around_0_deg, d$math_sd)
r_lab <- sprintf("r = %.2f, p = %.4f", ct$estimate, ct$p.value)
cat(r_lab, "\n")

# ---- plot ----
p <- ggplot(d, aes(sd_around_0_deg, math_sd, color = group, shape = group)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = alias), size = 3, show.legend = FALSE, seed = 1) +
  annotate("text", x = -Inf, y = Inf, label = r_lab,
           hjust = -0.15, vjust = 1.8, size = 4) +
  scale_color_manual(values = c(proCAF = "#E74095", restCAF = "#288D92")) +
  scale_shape_manual(values = c(proCAF = 16, restCAF = 15)) +
  labs(x = "UMA  sd_around_0_deg  (deg)",
       y = "Mathematica  sqrt(Variance)  (deg)",
       color = NULL, shape = NULL) +
  theme_classic(base_size = 12) +
  theme(legend.position = c(0.88, 0.12))

ggsave(file.path(ROOT, "results/fig13_uma_vs_mathematica.pdf"), p, width = 6, height = 5)
ggsave(file.path(ROOT, "results/fig13_uma_vs_mathematica.png"), p, width = 6, height = 5, dpi = 200)
