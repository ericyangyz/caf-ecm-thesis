# ECM thickness and nuclei per sample, proCAF vs restCAF (data/uma_outputs_thickness, data/uma_outputs_nuclei).
# out: results/thickness_nuclei_by_sample.csv, thickness_nuclei_tests.csv, fig14_thickness, fig15_nuclei
# group comes from the key below, not the export folder: 1207NA is restCAF but was exported as proCAF
library(tidyverse)
library(patchwork)

ROOT <- "."
PRO  <- "#E74095"
REST <- "#288D92"

# ---- sample -> group + de-identified alias ----
key <- tribble(
  ~sample,   ~group,    ~alias,
  "015",     "proCAF",  "P015T1",
  "229T1",   "proCAF",  "P229T1",
  "804N1",   "proCAF",  "P804N1",
  "804T1",   "proCAF",  "P804T1",
  "815",     "proCAF",  "P815T1",
  "908N1",   "proCAF",  "P908N1",
  "1207T1",  "proCAF",  "P207T1",
  "1108T1",  "restCAF", "P108T1",
  "1207NA",  "restCAF", "P207NA",     # was mislabeled proCAF
  "201N1",   "restCAF", "P201N1",
  "212T1",   "restCAF", "P212T1",
  "519T1",   "restCAF", "P519T1",
  "711N1",   "restCAF", "P711N1",
  "807",     "restCAF", "P807T1",
  "923",     "restCAF", "P923T1")

# ---- thickness: both files, group assigned afterwards ----
thick <- bind_rows(
    read_csv(file.path(ROOT, "data/uma_outputs_thickness/pro_Thickness_Summary.csv"),  show_col_types = FALSE),
    read_csv(file.path(ROOT, "data/uma_outputs_thickness/rest_Thickness_Summary.csv"), show_col_types = FALSE)) %>%
  # 1108T1.tiff has no "_", so strip the extension before the ROI number
  mutate(sample = str_remove(str_remove(File_Name, "\\.tiff$"), "_[0-9]+$")) %>%
  left_join(key, by = "sample")
stopifnot(nrow(thick) == 35, !anyNA(thick$group))

# ---- nuclei: pro_clustering_summary.csv is the combined export (all 35 ROIs) ----
nuc <- read_csv(file.path(ROOT, "data/uma_outputs_nuclei/pro_clustering_summary.csv"),
                show_col_types = FALSE) %>%
  mutate(sample = str_remove(str_remove(Image, "_nuclei_analysis$"), "_[0-9]+$"),
         pct_unclustered = XZ_Unclustered / Total_Nuclei * 100) %>%
  left_join(key, by = "sample")
stopifnot(nrow(nuc) == 35, !anyNA(nuc$group))

# ---- average ROIs within a sample (ROIs are not independent) ----
thick_s <- thick %>% group_by(sample, alias, group) %>%
  summarise(n_ROIs = n(), median_thickness = mean(Median), max_thickness = mean(Max), .groups = "drop")
nuc_s <- nuc %>% group_by(sample, alias, group) %>%
  summarise(n_ROIs = n(), total_nuclei = mean(Total_Nuclei), layers = mean(XZ_Clusters),
            pct_unclustered = mean(pct_unclustered), .groups = "drop")

by_sample <- full_join(thick_s, nuc_s, by = c("sample", "alias", "group", "n_ROIs"))
stopifnot(nrow(by_sample) == 15)
write_csv(by_sample, file.path(ROOT, "results/thickness_nuclei_by_sample.csv"))
print(count(by_sample, group))

# ---- Welch t-test + Wilcoxon on every metric ----
metrics <- c("median_thickness", "max_thickness", "total_nuclei", "layers", "pct_unclustered")
tests <- map_dfr(metrics, function(m) {
  a <- by_sample[[m]][by_sample$group == "proCAF"]
  b <- by_sample[[m]][by_sample$group == "restCAF"]
  tt <- t.test(a, b); wt <- suppressWarnings(wilcox.test(a, b))
  tibble(metric = m, n_pro = length(a), n_rest = length(b),
         mean_pro = mean(a), mean_rest = mean(b),
         p_welch = tt$p.value, p_wilcoxon = wt$p.value)
})
print(tests)
write_csv(tests, file.path(ROOT, "results/thickness_nuclei_tests.csv"))

# ---- boxplots ----
# one panel per metric so each keeps its own vertical y label; Welch p on top
by_sample <- by_sample %>%
  mutate(grp = factor(paste0(group, "\n(n = ", ifelse(group == "proCAF", 7, 8), ")"),
                      levels = c("proCAF\n(n = 7)", "restCAF\n(n = 8)")))

panel <- function(col, ylab) {
  pv <- tests$p_welch[tests$metric == col]
  top <- max(by_sample[[col]])
  rng <- diff(range(by_sample[[col]]))
  ggplot(by_sample, aes(grp, .data[[col]], fill = group)) +
    geom_boxplot(alpha = 0.35, outlier.shape = NA, width = 0.5) +
    geom_jitter(aes(color = group), width = 0.12, size = 2) +
    annotate("text", x = 1.5, y = top + 0.16 * rng, size = 3.1,
             label = sprintf("Welch's t-test, p = %.3f", pv)) +
    scale_fill_manual(values = c(proCAF = PRO, restCAF = REST)) +
    scale_color_manual(values = c(proCAF = PRO, restCAF = REST)) +
    expand_limits(y = top + 0.26 * rng) +
    labs(x = NULL, y = ylab) +
    theme_classic(base_size = 12) + theme(legend.position = "none")
}

p_thick <- panel("median_thickness", "Median thickness (px)") |
           panel("max_thickness",    "Max thickness (px)")
ggsave(file.path(ROOT, "results/fig14_thickness.pdf"), p_thick, width = 6.5, height = 4)
ggsave(file.path(ROOT, "results/fig14_thickness.png"), p_thick, width = 6.5, height = 4, dpi = 200)

p_nuc <- panel("total_nuclei",    "Total nuclei per ROI") |
         panel("layers",          "Layers (XZ clusters)") |
         panel("pct_unclustered", "Unclustered nuclei (%)")
ggsave(file.path(ROOT, "results/fig15_nuclei.pdf"), p_nuc, width = 9, height = 4)
ggsave(file.path(ROOT, "results/fig15_nuclei.png"), p_nuc, width = 9, height = 4, dpi = 200)
