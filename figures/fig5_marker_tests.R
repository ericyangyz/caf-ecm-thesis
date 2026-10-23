# Fig 5F: do the 10 CAF marker proteins differ by subtype or track fiber spread? BH over 10 markers per test.
# (a) proCAF vs restCAF, Welch t and Wilcoxon on log2(intensity + 1): 13 lines (6 pro, 7 rest; P220804T1 has no markers)
# (b) Spearman vs fiber spread (SD around dominant direction, deg): the 12 lines with fiber data (heatmap columns)
# Out: results/thesis/fig5/fig5F_marker_tests.csv
suppressPackageStartupMessages(library(readxl))
ROOT <- "."
stiff <- as.data.frame(read_excel(file.path(ROOT, "data/data_stiffness/items.xlsx"), sheet = "Sheet1"))
decaf <- read.delim(file.path(ROOT, "data/data_decaf/decaf_scores_20260922.tsv"), stringsAsFactors = FALSE, na.strings = "")
f5    <- read.csv(file.path(ROOT, "results/thesis/fig5/fig5F_marker_input.csv"), stringsAsFactors = FALSE)
MARKERS <- names(stiff)[5:14]
mk <- stiff[complete.cases(stiff[, MARKERS]), c("Line", MARKERS)]
mk$cafType <- factor(decaf$DeCAF[match(mk$Line, decaf$Line)], levels = c("proCAF", "restCAF"))
stopifnot(nrow(mk) == 13, !anyNA(mk$cafType), all(table(mk$cafType) == c(6, 7)), nrow(f5) == 12)

out <- data.frame()
for (m in MARKERS) {
  x <- log2(mk[[m]] + 1)
  tt <- t.test(x ~ mk$cafType); wt <- wilcox.test(x ~ mk$cafType, exact = FALSE)
  sp <- cor.test(f5$sd_around_0_deg, log2(f5[[m]] + 1), method = "spearman", exact = FALSE)
  out <- rbind(out, data.frame(
    marker = m,
    n_zero = sum(mk[[m]] == 0),
    median_proCAF_au  = median(mk[[m]][mk$cafType == "proCAF"]),
    median_restCAF_au = median(mk[[m]][mk$cafType == "restCAF"]),
    log2FC_pro_vs_rest = mean(x[mk$cafType == "proCAF"]) - mean(x[mk$cafType == "restCAF"]),
    p_welch = tt$p.value, p_wilcoxon = wt$p.value,
    rho_vs_spread = unname(sp$estimate), p_spearman = sp$p.value))
}
out$fdr_welch    <- p.adjust(out$p_welch, "BH")
out$fdr_wilcoxon <- p.adjust(out$p_wilcoxon, "BH")
out$fdr_spearman <- p.adjust(out$p_spearman, "BH")
out <- out[order(out$p_welch), ]
write.csv(out, file.path(ROOT, "results/thesis/fig5/fig5F_marker_tests.csv"), row.names = FALSE)
print(format(out[, c("marker", "n_zero", "median_proCAF_au", "median_restCAF_au", "log2FC_pro_vs_rest",
                     "p_welch", "fdr_welch", "p_wilcoxon", "rho_vs_spread", "p_spearman", "fdr_spearman")],
             digits = 3, scientific = FALSE), row.names = FALSE)
