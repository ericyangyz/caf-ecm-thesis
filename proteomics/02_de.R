# Differential expression of the matrisome proteome against matrix organization, with a
# continuous model (sd_around_0_deg) and a median-split model, set up as in
# CAF_multilayer_proteomics/regenerate_pool_quant.R
# In: RUN_DIR/design.rds. Out: RUN_DIR/de_continuous.csv, de_categorical.csv, de_fits.rds

library(limpa)
library(limma)
library(tidyverse)
library(org.Hs.eg.db)

here::i_am("proteomics/02_de.R")

# USER INPUTS (thresholds as in individual_layer_analysis.qmd)
# RUN_DIR is the only setting that differs between the two runs
RUN_DIR      <- "results/proteomics/run_old_keyence"
P_THRESHOLD  <- 0.05
LOGFC_THRESH <- 0.5

# STEP 1. Load the design, map UniProt to gene symbols

cat("Loading design object...\n")
d             <- readRDS(here::here(RUN_DIR, "design.rds"))
mat_quant_use <- d$mat_quant_use
samples       <- d$samples
org_range     <- d$org_range

stopifnot(
  identical(colnames(mat_quant_use$E), samples$sample),
  !is.null(mat_quant_use$other$standard.error),
  !any(grepl("^Cont_", rownames(mat_quant_use$E)))
)

all_uniprot <- rownames(mat_quant_use$E)
all_uniprot <- all_uniprot[!grepl(";", all_uniprot)]

uniprot_to_gene <- tryCatch(
  AnnotationDbi::select(org.Hs.eg.db,
                        keys   = unique(all_uniprot),
                        columns = "SYMBOL",
                        keytype = "UNIPROT") %>%
    dplyr::distinct(UNIPROT, .keep_all = TRUE) %>%
    dplyr::rename(Genes = SYMBOL),
  error = function(e) {
    cat("UniProt mapping failed:", e$message, "\n")
    data.frame(UNIPROT = character(), Genes = character())
  }
)
cat("UniProt mapping:", nrow(uniprot_to_gene), "/",
    length(unique(all_uniprot)), "mapped\n")

donor_block <- samples$donor

# STEP 2. Continuous model

cat("\n--- Continuous ---\n")

# dpcDE instead of lmFit: uses $dpc and the standard errors as weights, so imputed values
# count less; block = donor fits the within-donor correlation (duplicateCorrelation)
cont_fit  <- dpcDE(mat_quant_use, design = d$design_cont,
                   block = donor_block, plot = FALSE)
cont_fit2 <- eBayes(cont_fit)

cont_de <- topTable(cont_fit2, coef = "org_sd_deg", number = Inf,
                    sort.by = "none") %>%
  rownames_to_column(var = "UNIPROT") %>%
  left_join(uniprot_to_gene, by = "UNIPROT") %>%
  mutate(
    layer       = "Matrisome",
    model       = "continuous",
    sig_adj     = adj.P.Val < P_THRESHOLD,
    sig_nom     = P.Value < P_THRESHOLD,
    abs_logFC   = abs(logFC),
    logFC_range = logFC * org_range,
    # lower sd = more organized, so a negative slope means up in organized
    direction   = case_when(
      sig_adj & logFC_range < -LOGFC_THRESH ~ "Up in organized",
      sig_adj & logFC_range >  LOGFC_THRESH ~ "Up in disorganized",
      TRUE ~ "NS"
    )
  )

write_csv(cont_de, here::here(RUN_DIR, "de_continuous.csv"))
cat("DE continuous: adj-sig", sum(cont_de$sig_adj),
    "| nom-sig", sum(cont_de$sig_nom), "\n")

# STEP 3. Categorical model (median split)

cat("\n--- Categorical (median split) ---\n")

cat_fit  <- dpcDE(mat_quant_use, design = d$design_cat,
                  block = donor_block, plot = FALSE)

contrast_cat <- makeContrasts(organized - disorganized, levels = d$design_cat)
cat_fit2 <- contrasts.fit(cat_fit, contrast_cat)
cat_fit2 <- eBayes(cat_fit2)

cat_de <- topTable(cat_fit2, number = Inf, sort.by = "none") %>%
  rownames_to_column(var = "UNIPROT") %>%
  left_join(uniprot_to_gene, by = "UNIPROT") %>%
  mutate(
    layer     = "Matrisome",
    model     = "categorical",
    sig_adj   = adj.P.Val < P_THRESHOLD,
    sig_nom   = P.Value < P_THRESHOLD,
    abs_logFC = abs(logFC),
    direction = case_when(
      sig_adj & logFC >  LOGFC_THRESH ~ "Up in organized",
      sig_adj & logFC < -LOGFC_THRESH ~ "Up in disorganized",
      TRUE ~ "NS"
    )
  )

write_csv(cat_de, here::here(RUN_DIR, "de_categorical.csv"))
cat("DE categorical: adj-sig", sum(cat_de$sig_adj),
    "| nom-sig", sum(cat_de$sig_nom), "\n")

# STEP 4. Control: the continuous fit without the donor block

cat("\n--- Control: donor block dropped ---\n")

# not a result: without the block the 2-4 MS runs of a donor count as independent
noblock_fit <- eBayes(dpcDE(mat_quant_use, design = d$design_cont,
                            plot = FALSE))
noblock_de  <- topTable(noblock_fit, coef = "org_sd_deg", number = Inf,
                        sort.by = "none")

cat("Ignoring donor:", sum(noblock_de$adj.P.Val < P_THRESHOLD),
    "adj-sig vs", sum(cont_de$sig_adj), "with the block\n")

# STEP 5. Save fits

write_rds(list(cont_fit = cont_fit2, cat_fit = cat_fit2,
               cont_de = cont_de, cat_de = cat_de,
               noblock_de = noblock_de),
          here::here(RUN_DIR, "de_fits.rds"))
cat("\nSaved:", here::here(RUN_DIR, "de_fits.rds"), "\n")

cont_de %>%
  filter(sig_adj) %>%
  arrange(adj.P.Val) %>%
  dplyr::select(Genes, UNIPROT, logFC, logFC_range, adj.P.Val,
                PropObs, direction) %>%
  print(row.names = FALSE, digits = 3)
