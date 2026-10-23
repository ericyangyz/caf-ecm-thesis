# PDF figures of the matrisome hits (from 03_matrisome_intersection.R) of the old-Keyence run.
# The contrast is matrix organization (organized - disorganized), not proCAF vs restCAF;
# the two are confounded (proCAF 14 disorganized / 3 organized; restCAF 2 / 16).
# Out: results/proteomics/run_old_keyence/fig_*.pdf

library(tidyverse)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)

# USER INPUTS
RUN_DIR      <- "results/proteomics/run_old_keyence"
P_THRESHOLD  <- 0.05
LOGFC_THRESH <- 0.5
N_SCATTER    <- 6

SIG_COLORS <- c("Up in organized" = "#377EB8", "Up in disorganized" = "#E41A1C",
                NS = "grey78")
# Division: violet / ochre, unrelated to the CAF colors, differ in lightness for color-blind readers
DIV_COLORS <- c("Core matrisome" = "#5E4B8B", "Matrisome-associated" = "#D9A441")
BATCH_COLORS <- c(PC1514 = "#BDBDBD", PC1556 = "#6E6E6E")
GROUP_COLORS <- c(proCAF = "#E74095", restCAF = "#288D92")     # proCAF / restCAF colors

# de-identified donor aliases for figures
ALIAS <- c(P201015T1 = "P015T1", P220804N1 = "P804N1", P220804T1 = "P804T1",
           P220815T1 = "P815T1", P220908N1 = "P908N1", P231207T1 = "P207T1",
           P240229T1 = "P229T1", P140807T1 = "P807T1", P210923T1 = "P923T1",
           P220201N1 = "P201N1", P220212T1 = "P212T1", P220519T1 = "P519T1",
           P221108T1 = "P108T1", P231207NA = "P207NA", P240711N1 = "P711N1",
           P220118N1 = "P118N1", P240509N1 = "P509N1")

base_thm <- theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey30"))

# STEP 1. Load the run and the matrisome intersection
d       <- readRDS(file.path(RUN_DIR, "design.rds"))
fits    <- readRDS(file.path(RUN_DIR, "de_fits.rds"))
samples <- d$samples
quant   <- d$mat_quant_use
cont_de <- fits$cont_de
cat_de  <- fits$cat_de
ov      <- read.csv(file.path(RUN_DIR, "matrisome_intersection.csv"))

stopifnot(nrow(cont_de) == 7578, nrow(ov) == 23,
          identical(colnames(quant$E), samples$sample),
          all(ov$UNIPROT %in% rownames(quant$E)))
cat(sprintf("STEP 1: %d proteins tested, %d FDR-significant, %d of those matrisome\n",
            nrow(cont_de), sum(cont_de$sig_adj), nrow(ov)))

# annotate every tested protein (denominator for the enrichment test), same cleaning as
# 03_matrisome_intersection.R
load("data/data_proteomics/matrisome.list.rda")
hs <- matrisome.list$human
hs <- hs[!is.na(hs$gene) & nzchar(hs$gene) &
         !grepl("^[0-9]+$", hs$gene) & !grepl("^ENSG[0-9]", hs$gene) &
         !grepl("^[NXY][PM]_", hs$gene) & !grepl("\\.", hs$gene) &
         !grepl("^[A-Z][0-9][A-Z0-9]{3}[0-9]$", hs$gene), ]
hs <- hs[!duplicated(hs$gene), ]
cont_de$MatDivision <- hs$category[match(cont_de$Genes, hs$gene)]
de_all_matrisome <- cont_de[!is.na(cont_de$MatDivision), ]
cat(sprintf("STEP 2: %d of %d tested proteins are matrisome-annotated\n",
            nrow(de_all_matrisome), nrow(cont_de)))

cont_de$Division <- ov$Division[match(cont_de$UNIPROT, ov$UNIPROT)]
cat_de$Division  <- ov$Division[match(cat_de$UNIPROT,  ov$UNIPROT)]

# STEP 2. Volcanoes, matrisome hits drawn on top and labeled
# x = logFC across the observed organization span (not per degree), negated so
# positive = higher in organized
cont_de$x <- -cont_de$logFC_range
lab_cont  <- cont_de[!is.na(cont_de$Division), ]

p1 <- ggplot(cont_de, aes(x, -log10(adj.P.Val))) +
  geom_point(aes(color = direction), alpha = 0.45, size = 1) +
  geom_point(data = lab_cont, aes(fill = Division), shape = 21, size = 2.6,
             color = "white", stroke = 0.5) +
  geom_hline(yintercept = -log10(P_THRESHOLD), linetype = "dashed", color = "grey55") +
  geom_vline(xintercept = c(-LOGFC_THRESH, LOGFC_THRESH), linetype = "dashed",
             color = "grey55") +
  geom_text_repel(data = lab_cont, aes(label = Genes), size = 3, seed = 1,
                  max.overlaps = Inf, min.segment.length = 0.2, color = "grey15") +
  scale_color_manual(values = SIG_COLORS, guide = "none") +
  scale_fill_manual(values = DIV_COLORS) +
  labs(x = "log2 FC across the observed organization range (organized - disorganized)",
       y = "-log10 adjusted P", fill = NULL,
       title = "Matrix organization, continuous model",
       subtitle = sprintf("%d of %d FDR-significant proteins are matrisome",
                          nrow(ov), sum(cont_de$sig_adj))) +
  base_thm + theme(legend.position = "top")
ggsave(file.path(RUN_DIR, "fig_volcano_continuous.pdf"), p1, width = 7.5, height = 6.4,
       device = cairo_pdf)

lab_cat <- cat_de[!is.na(cat_de$Division) & cat_de$sig_adj, ]
p2 <- ggplot(cat_de, aes(logFC, -log10(adj.P.Val))) +
  geom_point(aes(color = direction), alpha = 0.45, size = 1) +
  geom_point(data = cat_de[!is.na(cat_de$Division), ], aes(fill = Division),
             shape = 21, size = 2.4, color = "white", stroke = 0.5) +
  geom_hline(yintercept = -log10(P_THRESHOLD), linetype = "dashed", color = "grey55") +
  geom_vline(xintercept = c(-LOGFC_THRESH, LOGFC_THRESH), linetype = "dashed",
             color = "grey55") +
  geom_text_repel(data = lab_cat, aes(label = Genes), size = 3, seed = 1,
                  max.overlaps = Inf, color = "grey15") +
  scale_color_manual(values = SIG_COLORS, guide = "none") +
  scale_fill_manual(values = DIV_COLORS) +
  labs(x = "log2 FC (organized - disorganized)", y = "-log10 adjusted P", fill = NULL,
       title = "Matrix organization, median-split model",
       subtitle = sprintf("%d FDR-significant proteins", sum(cat_de$sig_adj))) +
  base_thm + theme(legend.position = "top")
ggsave(file.path(RUN_DIR, "fig_volcano_categorical.pdf"), p2, width = 7.5, height = 6.4,
       device = cairo_pdf)

# STEP 3. Heatmap of the 23 matrisome hits
# columns = the 12 donors (mean of their MS runs), in organization order, not clustered
donor_ids <- levels(samples$donor)
E_sig   <- quant$E[ov$UNIPROT, , drop = FALSE]
E_donor <- sapply(donor_ids, function(dd) rowMeans(E_sig[, samples$donor == dd, drop = FALSE]))

info <- distinct(samples, donor, org_sd_deg, batch, group)
info <- info[match(donor_ids, info$donor), ]
ord  <- order(info$org_sd_deg)

hm <- t(scale(t(E_donor[, ord, drop = FALSE])))
rownames(hm) <- ifelse(is.na(ov$Genes) | ov$Genes == "", ov$UNIPROT, ov$Genes)
stopifnot(all(colnames(hm) %in% names(ALIAS)))      # every donor needs an alias
colnames(hm) <- ALIAS[colnames(hm)]

col_annot <- data.frame(`CAF type` = info$group[ord], Batch = info$batch[ord],
                        `SD around 0 (deg)` = info$org_sd_deg[ord],
                        row.names = colnames(hm), check.names = FALSE)
row_annot <- data.frame(Division = ov$Division, row.names = rownames(hm))

pheatmap(hm,
         annotation_col = col_annot, annotation_row = row_annot,
         annotation_colors = list(Division = DIV_COLORS, `CAF type` = GROUP_COLORS,
                                  Batch = BATCH_COLORS,
                                  `SD around 0 (deg)` = brewer.pal(9, "Greys")[2:7]),
         cluster_cols = FALSE, cluster_rows = TRUE,
         clustering_distance_rows = "euclidean", clustering_method = "ward.D2",
         show_rownames = TRUE, show_colnames = TRUE, fontsize_row = 8,
         main = "Matrisome proteins among the organization-associated hits",
         color = colorRampPalette(rev(brewer.pal(7, "RdBu")))(100),
         breaks = seq(-2, 2, length.out = 101), border_color = NA,
         filename = file.path(RUN_DIR, "fig_matrisome_heatmap.pdf"),
         width = 8.5, height = 6.5)

# STEP 4. Division composition of each direction
comp <- ov %>%
  count(direction, Division) %>%
  group_by(direction) %>%
  mutate(prop = n / sum(n), tot = sum(n)) %>%
  ungroup() %>%
  mutate(direction = factor(direction,
           levels = c("Up in disorganized", "Up in organized"),
           labels = c(sprintf("Up in disorganized\n(n = %d)",
                              sum(ov$direction == "Up in disorganized")),
                      sprintf("Up in organized\n(n = %d)",
                              sum(ov$direction == "Up in organized")))))

p4 <- ggplot(comp, aes(direction, prop, fill = Division)) +
  geom_col(width = 0.55) +
  geom_text(aes(label = ifelse(prop > 0.06, sprintf("%d\n%.0f%%", n, 100 * prop), "")),
            position = position_stack(vjust = 0.5), size = 3.4, color = "white",
            fontface = "bold", lineheight = 0.95) +
  scale_fill_manual(values = DIV_COLORS) +
  scale_y_continuous(labels = scales::percent, expand = expansion(c(0, 0.02))) +
  labs(x = NULL, y = "Proportion of matrisome hits", fill = NULL,
       title = "Matrisome composition of each side",
       subtitle = "the organized side carries only two proteins, so its proportions are fragile") +
  base_thm + theme(legend.position = "top")
ggsave(file.path(RUN_DIR, "fig_matrisome_composition.pdf"), p4, width = 5.6, height = 5.2,
       device = cairo_pdf)

# STEP 5. Per-donor values behind the top matrisome hits
top <- head(ov[order(ov$adj.P.Val), ], N_SCATTER)
raw <- lapply(seq_len(nrow(top)), function(i) {
  data.frame(gene = top$Genes[i], Division = top$Division[i],
             donor = as.character(samples$donor),
             org = samples$org_sd_deg, y = quant$E[top$UNIPROT[i], ])
}) %>% bind_rows() %>%
  group_by(gene, Division, donor, org) %>% summarise(y = mean(y), .groups = "drop") %>%
  mutate(gene = factor(gene, levels = top$Genes))

p5 <- ggplot(raw, aes(org, y)) +
  geom_smooth(method = "lm", se = FALSE, color = "grey65", linewidth = 0.6) +
  geom_point(aes(color = Division), size = 2.2) +
  facet_wrap(~gene, scales = "free_y", nrow = 2) +
  scale_color_manual(values = DIV_COLORS) +
  labs(x = "Fiber spread, SD around dominant direction (°)   — lower = more organized",
       y = "log2 intensity, donor mean", color = NULL,
       title = "The strongest matrisome hits, one point per donor",
       subtitle = "each panel is 12 donors; a single extreme donor at 51° anchors the slope") +
  base_thm + theme(legend.position = "top", strip.background = element_blank(),
                   strip.text = element_text(face = "bold"))
ggsave(file.path(RUN_DIR, "fig_top_hits_scatter.pdf"), p5, width = 9, height = 5.8,
       device = cairo_pdf)

cat("STEP 5: wrote 5 PDFs to", RUN_DIR, "\n")
print(comp[, c("direction", "Division", "n", "prop")])

# STEP 6. Matrisome enrichment among the hits (Fisher test)
cont_de$is_matrisome <- cont_de$UNIPROT %in% ov$UNIPROT |
  cont_de$UNIPROT %in% de_all_matrisome$UNIPROT
tab <- table(matrisome = cont_de$is_matrisome, hit = cont_de$sig_adj)
ft  <- fisher.test(tab)

enr <- data.frame(
  grp  = c("Matrisome", "Not matrisome"),
  hits = c(tab["TRUE", "TRUE"], tab["FALSE", "TRUE"]),
  tot  = c(sum(tab["TRUE", ]), sum(tab["FALSE", ])))
enr$rate <- enr$hits / enr$tot
# Wilson interval; the normal approximation fails at a 2% rate
z <- qnorm(0.975)
enr$lo <- (enr$hits + z^2/2)/(enr$tot + z^2) -
  z/(enr$tot + z^2) * sqrt(enr$hits*(enr$tot - enr$hits)/enr$tot + z^2/4)
enr$hi <- (enr$hits + z^2/2)/(enr$tot + z^2) +
  z/(enr$tot + z^2) * sqrt(enr$hits*(enr$tot - enr$hits)/enr$tot + z^2/4)
enr$grp <- factor(enr$grp, levels = c("Matrisome", "Not matrisome"))

# print p as "1.9 x 10^-5" rather than "1.89e-05"
sup <- c("\u2070","\u00b9","\u00b2","\u00b3","\u2074","\u2075","\u2076","\u2077","\u2078","\u2079")
ex  <- floor(log10(ft$p.value))
man <- ft$p.value / 10^ex
p_lab <- paste0(sprintf("%.1f \u00d7 10\u207b", man),
                paste(sup[as.integer(strsplit(as.character(abs(ex)), "")[[1]]) + 1], collapse = ""))

# light neutral greys: the error bar stays visible and the bars do not read as Division colors
p6 <- ggplot(enr, aes(grp, rate, fill = grp)) +
  geom_col(width = 0.5, color = "grey35", linewidth = 0.3) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.12, linewidth = 0.55, color = "#1A1A1A") +
  geom_text(aes(y = hi, label = sprintf("%.1f%%", 100 * rate)),
            vjust = -1.5, size = 3.6, fontface = "bold") +
  geom_text(aes(y = hi, label = sprintf("%s / %s", trimws(format(hits, big.mark = ",")),
                                        trimws(format(tot, big.mark = ",")))),
            vjust = -0.3, size = 3.1, color = "grey35") +
  scale_fill_manual(values = c(Matrisome = "#8A8F98", `Not matrisome` = "#D6D9DD"),
                    guide = "none") +
  scale_y_continuous(labels = scales::percent, expand = expansion(c(0, 0.22))) +
  labs(x = NULL, y = "FDR-significant proteins",
       title = "The matrisome is enriched among\norganization-associated proteins",
       subtitle = sprintf("Fisher odds ratio %.2f (95%% CI %.2f\u2013%.2f), p = %s",
                          ft$estimate, ft$conf.int[1], ft$conf.int[2], p_lab)) +
  base_thm
ggsave(file.path(RUN_DIR, "fig_matrisome_enrichment.pdf"), p6,
       width = 4.4, height = 5, device = cairo_pdf)

cat(sprintf("STEP 6: matrisome %d/%d = %.1f%%, other %d/%d = %.1f%%, OR %.2f (%.2f-%.2f), p = %.2e\n",
            enr$hits[1], enr$tot[1], 100*enr$rate[1], enr$hits[2], enr$tot[2], 100*enr$rate[2],
            ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value))
