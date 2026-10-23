# Which FDR-significant proteins of the full 7578-protein run are matrisome-annotated.
# FDR is controlled over the whole proteome; matrisome membership only describes the hit list.
# In: RUN_DIR/de_continuous.csv. Out: RUN_DIR/matrisome_intersection.csv and .png

library(ggplot2)

# USER INPUTS
RUN_DIR        <- "results/proteomics/run_old_keyence"
matrisome_file <- "data/data_proteomics/matrisome.list.rda"
P_THRESHOLD    <- 0.05

GROUP_COLORS <- c(proCAF = "#E74095", restCAF = "#288D92")   # proCAF / restCAF colors

# STEP 1. The full-proteome DE result
de <- read.csv(file.path(RUN_DIR, "de_continuous.csv"))
stopifnot(nrow(de) == 7578, all(c("Genes", "adj.P.Val", "PropObs") %in% colnames(de)))

hits <- de[de$adj.P.Val < P_THRESHOLD, ]
cat("STEP 1:", nrow(de), "proteins tested,", nrow(hits), "FDR-significant\n")

# STEP 2. Matrisome annotation
# cleaning as in unsupervised_discovery_analysis.qmd (CAF_multilayer_proteomics, L1529-1549):
# drop non-symbol ids; Division = category, Category = family
load(matrisome_file)
hs <- matrisome.list$human
hs <- hs[!is.na(hs$gene) & nzchar(hs$gene) &
         !grepl("^[0-9]+$", hs$gene) & !grepl("^ENSG[0-9]", hs$gene) &
         !grepl("^[NXY][PM]_", hs$gene) & !grepl("\\.", hs$gene) &
         !grepl("^[A-Z][0-9][A-Z0-9]{3}[0-9]$", hs$gene), ]
hs <- hs[!duplicated(hs$gene), ]

de$Division <- hs$category[match(de$Genes, hs$gene)]
de$Category <- hs$family[match(de$Genes, hs$gene)]

cat("STEP 2:", sum(!is.na(de$Category)), "of", nrow(de), "tested proteins are matrisome\n")

# STEP 3. The intersection
ov <- de[de$adj.P.Val < P_THRESHOLD & !is.na(de$Category), ]
ov <- ov[order(ov$adj.P.Val), ]

cat("STEP 3:", nrow(ov), "of the", nrow(hits), "hits are matrisome\n\n")
print(ov[, c("Genes", "Division", "Category", "logFC_range", "adj.P.Val",
             "PropObs", "direction")], row.names = FALSE, digits = 3)

write.csv(ov[, c("UNIPROT", "Genes", "Division", "Category", "logFC", "logFC_range",
                 "P.Value", "adj.P.Val", "PropObs", "direction")],
          file.path(RUN_DIR, "matrisome_intersection.csv"), row.names = FALSE)

# is the matrisome over-represented among the hits?
tab <- table(matrisome = !is.na(de$Category), hit = de$adj.P.Val < P_THRESHOLD)
cat("\n"); print(tab)
ft <- fisher.test(tab)
cat(sprintf("\nFisher odds ratio %.2f (95%% CI %.2f - %.2f), p = %.3f\n",
            ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value))
cat(sprintf("hit rate: matrisome %.1f%% vs rest %.1f%%\n",
            100 * tab["TRUE","TRUE"] / sum(tab["TRUE",]),
            100 * tab["FALSE","TRUE"] / sum(tab["FALSE",])))

# STEP 4. Figure
ov$Genes <- factor(ov$Genes, levels = rev(ov$Genes))

p <- ggplot(ov, aes(x = logFC_range, y = Genes, fill = Division)) +
  geom_col() +
  geom_vline(xintercept = 0, color = "grey40") +
  scale_fill_manual(values = c("Core matrisome" = "#E74095",
                               "Matrisome-associated" = "#288D92")) +
  labs(title = paste(nrow(ov), "matrisome proteins among the",
                     nrow(hits), "FDR-significant hits"),
       subtitle = "Positive = higher in a DISORGANIZED matrix. Ordered by adjusted p.",
       x = "log2 FC over the observed linearity range", y = NULL, fill = NULL) +
  theme_bw()

ggsave(file.path(RUN_DIR, "matrisome_intersection.png"), p,
       width = 7.5, height = 6, dpi = 300)
cat("\nWrote matrisome_intersection.csv and .png to", RUN_DIR, "\n")
