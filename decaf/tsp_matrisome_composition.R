# How many of the 18 DeCAF TSP genes are matrisome / ECM genes (why we study the ECM).
# Out: results/decaf/tsp_ecm_genes.pdf/.png, tsp_ecm_composition.pdf/.png

library(ggplot2)
library(patchwork)

out <- "./results/decaf"
dir.create(out, showWarnings = FALSE)

pro  <- "#E74095"    # proCAF / restCAF colors
rest <- "#288D92"
pal  <- c("Core matrisome"       = "#B99338",
          "Matrisome-associated" = "#E0C97F",
          "ECM effector"         = "#F2E6C4",
          "Non-ECM"              = "#EFEFEC")

theme_set(theme_classic(base_size = 12, base_family = "Arial") +
            theme(plot.title    = element_text(hjust = 0.5),
                  plot.subtitle = element_text(hjust = 0.5, size = 9, colour = "grey30"),
                  plot.tag      = element_text(face = "bold")))

# the 9 DeCAF TSP pairs; TSP8 is CNIH3 (written GNIH3 in the DeCAF table)
decaf <- data.frame(
  tsp  = rep(paste0("TSP", 1:9), each = 2),
  side = rep(c("proCAF", "restCAF"), 9),
  gene = c("IGFL2", "CHRDL1", "NOX4", "OGN", "VSNL1", "PI16", "BICD1", "ANK2", "NPR3", "ABCA8",
           "ETV1", "TGFBR3", "ITGA11", "FBLN5", "CNIH3", "SCARA5", "COL11A1", "KIAA1217"))

# core / associated from the matrisome list; ITGA11 and TGFBR3 are not in it but are
# GO-annotated to cell-matrix adhesion and ECM organization
load("./data/data_proteomics/matrisome.list.rda")
hs <- matrisome.list$human
decaf$class <- hs$category[match(decaf$gene, hs$gene)]
decaf$class[decaf$gene %in% c("ITGA11", "TGFBR3")] <- "ECM effector"
decaf$class[is.na(decaf$class)] <- "Non-ECM"
decaf$class <- factor(decaf$class, levels = names(pal))
table(decaf$class, decaf$side)

ecm <- decaf[decaf$class != "Non-ECM", ]

# ---- table of the 18 genes ----
p_tab <- ggplot(decaf, aes(side, factor(tsp, levels = paste0("TSP", 9:1)), fill = class)) +
  geom_tile(colour = "white", linewidth = 1.4) +
  geom_text(aes(label = gene, fontface = ifelse(class == "Non-ECM", "plain", "bold")),
            size = 3.1, colour = "#16211F", family = "Arial") +
  scale_fill_manual(values = pal, name = NULL) +
  scale_x_discrete(position = "top", labels = c("proCAF gene", "restCAF gene")) +
  labs(x = NULL, y = NULL, title = "The 18 DeCAF TSP genes",
       subtitle = paste0(nrow(ecm), " of 18 are matrisome or ECM effector genes (bold)")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(face = "bold", colour = "black"),
        legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2))

ggsave(file.path(out, "tsp_ecm_genes.pdf"), p_tab, width = 4.2, height = 5.2, device = cairo_pdf)
ggsave(file.path(out, "tsp_ecm_genes.png"), p_tab, width = 4.2, height = 5.2, dpi = 300)

# ---- composition ----
n <- table(ecm$side)
comp <- rbind(data.frame(group = "All genes\n(n = 18)", class = decaf$class),
              data.frame(group = paste0(ecm$side, " side\n(n = ", n[ecm$side], ")"), class = ecm$class))
comp$group <- factor(comp$group, levels = c("All genes\n(n = 18)", paste0(names(n), " side\n(n = ", n, ")")))

p_a <- ggplot(comp, aes(group, fill = class)) +
  geom_bar(position = "fill", width = 0.55, colour = "grey35", linewidth = 0.3) +
  geom_text(stat = "count", position = position_fill(vjust = 0.5), size = 3.2, colour = "#16211F",
            fontface = "bold", family = "Arial",
            aes(label = after_stat(paste0(count, " (", round(100 * count / ave(count, x, FUN = sum)), "%)")))) +
  scale_fill_manual(values = pal, name = NULL) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Proportion of genes", title = "Gene type",
       subtitle = "all 18 genes, then the ECM-related genes on each side")

p_b <- ggplot(ecm, aes(paste0("ECM-related\n(n = ", nrow(ecm), ")"), fill = side)) +
  geom_bar(position = "fill", width = 0.55, colour = "grey35", linewidth = 0.3) +
  geom_text(stat = "count", position = position_fill(vjust = 0.5), size = 3.2, colour = "white",
            fontface = "bold", family = "Arial",
            aes(label = after_stat(paste0(count, " (", round(100 * count / sum(count)), "%)")))) +
  scale_fill_manual(values = c(proCAF = pro, restCAF = rest), name = NULL) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = NULL, title = "Side of the pair", subtitle = "ECM-related genes only")

p <- p_a + p_b + plot_layout(widths = c(3, 1)) + plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom")

ggsave(file.path(out, "tsp_ecm_composition.pdf"), p, width = 8.5, height = 5.2, device = cairo_pdf)
ggsave(file.path(out, "tsp_ecm_composition.png"), p, width = 8.5, height = 5.2, dpi = 300)
