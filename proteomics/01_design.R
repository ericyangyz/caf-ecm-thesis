# Attach a per-donor fiber organization value to each proteomics sample and build the DE designs.
# Organization = sd_around_0_deg (deg); lower sd = more aligned = more organized.
# In:  results/old_keyence_alignment_per_sample.csv (or results/fib), data/data_proteomics/
# Out: RUN_DIR/design.rds, donor_organization.csv, design_samples.csv

setwd(".")

# USER INPUTS

# ORG_SOURCE is the only difference between the two runs (STEP 4 onward is shared):
# "fib" = large-field movies, 3-4 fields x 25 ROIs per donor; "old_keyence" = old Keyence
# 20x stills, 1 image per donor (the run used for thesis Fig 6)
ORG_SOURCE <- "old_keyence"
RUN_DIR    <- "results/proteomics/run_old_keyence"

# TRUE keeps only matrisome proteins before fitting (smaller BH denominator)
MATRISOME_ONLY <- FALSE

fib_dir    <- "results/fib"
align_file <- "results/guidance/alignment_scales.csv"
old_file   <- "./results/old_keyence_alignment_per_sample.csv"
matrisome_file <- "data/data_proteomics/matrisome.list.rda"
quant_file <- "data/data_proteomics/matrisome_quant.rds"
meta_file  <- "data/data_proteomics/metadata_maps.rds"

ORG_METRIC <- "sd_around_0_deg"
DROP_FIELD <- "P170713_field_13"

# expected counts for the stopifnot checks below
if (ORG_SOURCE == "fib") {
  EXP_UNITS   <- 54          # fields read
  EXP_DONORS  <- 14          # donors with imaging
  EXP_NFIELDS <- c(3, 4)     # fields averaged per donor
} else {
  EXP_UNITS   <- 13          # images with a cfECM donor
  EXP_DONORS  <- 13
  EXP_NFIELDS <- 1
}
EXP_MATRISOME    <- 396      # proteins kept when MATRISOME_ONLY is TRUE
EXP_SAMPLES      <- 35       # MS runs left after joining to proteomics
EXP_PROT_DONORS  <- 12       # donors with both imaging and proteomics
EXP_SPLIT        <- 6        # donors on the organized side of the median

if (ORG_SOURCE == "fib") {

  # STEP 1. Read the ROI summaries of the 55 fields
  # glob the fixed results/fib/<category>/<field>/ layout; list.files(recursive = TRUE)
  # would walk all 341,729 per-frame files
  sum_files <- Sys.glob(file.path(fib_dir, "*", "*", "FINAL_summary_*.csv"))
  cat("found", length(sum_files), " FINAL_summary files\n")

  fields <- do.call(rbind, lapply(sum_files, read.csv))
  dim(fields)

  # STEP 2. Attach the cfECM donor, drop the duplicate field

  fields$key <- sub("^PanT6", "PancT6", fields$sample_name) # spelling used in the alignment file

  align        <- read.csv(align_file)
  fields$donor <- align$donor[match(fields$key, align$field)]

  cat("Fields with no donor:", fields$sample_name[is.na(fields$donor)], "\n")

  fields <- fields[fields$sample_name != DROP_FIELD, ]
  fields <- fields[!is.na(fields$donor), ]

  stopifnot(
    nrow(fields) == EXP_UNITS,                     # 55 files minus the duplicate
    length(unique(fields$donor)) == EXP_DONORS,
    all(fields$n_roi == 25),
    is.numeric(fields[[ORG_METRIC]]),
    !anyNA(fields[[ORG_METRIC]])
  )
  cat("STEP 2 OK:", nrow(fields), "fields,",
      length(unique(fields$donor)), "donors\n")

  # STEP 3. Average the 3-4 fields of each donor (proteomics has one matrix per donor)
  org <- fields[[ORG_METRIC]]

  donors <- aggregate(list(org_sd_deg = org),
                      by  = list(donor = fields$donor),
                      FUN = mean)

  donors$n_fields <- as.vector(table(fields$donor)[donors$donor])

  donors <- donors[order(donors$org_sd_deg), ]

  # share of variance between vs within donors
  ss          <- summary(aov(org ~ fields$donor))[[1]][["Sum Sq"]]
  pct_between <- 100 * ss[1] / sum(ss)
  cat(sprintf("Between-donor variance: %.1f%% of total\n", pct_between))

  stopifnot(
    nrow(donors) == EXP_DONORS,
    all(donors$n_fields %in% EXP_NFIELDS),
    !anyNA(donors$org_sd_deg),
    pct_between > 50            # fields of one donor must agree with each other
  )
  cat("STEP 3 OK:", nrow(donors), "donors,",
      sprintf("%.1f - %.1f deg\n", min(donors$org_sd_deg), max(donors$org_sd_deg)))

} else {

  # STEP 1-3 (old Keyence). One 20x still per donor, read from the per-sample table
  old <- read.csv(old_file)
  cat("found", nrow(old), "old-Keyence samples\n")

  # images with no cfECM donor in tracking_metadata.xlsx cannot join the proteomics
  cat("No donor, dropped:", old$sample[is.na(old$donor)], "\n")
  old <- old[!is.na(old$donor), ]

  stopifnot(
    nrow(old) == EXP_UNITS,
    !anyDuplicated(old$donor),              # one image per donor
    all(abs(old$total_perc - 100) < 1e-6),  # histogram weights are percentages
    is.numeric(old[[ORG_METRIC]]),
    !anyNA(old[[ORG_METRIC]])
  )
  cat("STEP 2 OK:", nrow(old), "images,", length(unique(old$donor)), "donors\n")

  donors <- data.frame(donor      = old$donor,
                       org_sd_deg = old[[ORG_METRIC]],
                       n_fields   = 1L)
  donors <- donors[order(donors$org_sd_deg), ]

  # no between/within-donor variance check: one image per donor
  stopifnot(
    nrow(donors) == EXP_DONORS,
    all(donors$n_fields %in% EXP_NFIELDS),
    !anyNA(donors$org_sd_deg)
  )
  cat("STEP 3 OK:", nrow(donors), "donors,",
      sprintf("%.1f - %.1f deg\n", min(donors$org_sd_deg), max(donors$org_sd_deg)))
}

# STEP 4. Attach the organization value to each proteomics sample
q  <- readRDS(quant_file)
m  <- readRDS(meta_file)
mm <- m$matrisome

# matrisome_quant.rds still contains the cRAP contaminants (Cont_ rows); drop them so
# the BH denominator is the real proteome
keep_prot <- !grepl("^Cont_", rownames(q$E))
cat("Dropped", sum(!keep_prot), "Cont_ rows,", sum(keep_prot), "proteins left\n")
stopifnot(sum(keep_prot) == 7578)
q <- q[keep_prot, ]

# optional matrisome filter: strip non-symbol ids from matrisome.list (same five regexes as
# unsupervised_discovery_analysis.qmd in CAF_multilayer_proteomics), map UniProt -> symbol,
# then join symbol to symbol
if (MATRISOME_ONLY) {
  suppressPackageStartupMessages({
    library(org.Hs.eg.db); library(AnnotationDbi)
  })
  load(matrisome_file)                       # supplies matrisome.list
  hs <- matrisome.list$human
  stopifnot(all(c("gene", "category", "family") %in% colnames(hs)))

  hs <- hs[!is.na(hs$gene) & nzchar(hs$gene) &
           !grepl("^[0-9]+$",                     hs$gene) &   # Entrez numbers
           !grepl("^ENSG[0-9]",                   hs$gene) &   # Ensembl gene ids
           !grepl("^[NXY][PM]_",                  hs$gene) &   # RefSeq accessions
           !grepl("\\.",                          hs$gene) &   # versioned ids
           !grepl("^[A-Z][0-9][A-Z0-9]{3}[0-9]$", hs$gene), ]  # UniProt accessions
  cat("Matrisome DB:", nrow(hs), "symbols,",
      length(unique(hs$category)), "divisions,",
      length(unique(hs$family)), "categories\n")

  acc <- rownames(q$E)
  sym_map <- suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db, keys = acc, columns = "SYMBOL", keytype = "UNIPROT"))
  sym_map <- sym_map[!duplicated(sym_map$UNIPROT), ]
  sym     <- sym_map$SYMBOL[match(acc, sym_map$UNIPROT)]

  keep_mat <- !is.na(sym) & sym %in% hs$gene
  cat("Matrisome filter:", sum(keep_mat), "of", length(acc), "proteins kept\n")
  stopifnot(sum(keep_mat) == EXP_MATRISOME)
  q <- q[keep_mat, ]
}

# attach donor; metadata_maps covers the 40 biological samples only, so the 6 REF_POOL
# channels come back NA and only non-pool samples must match
idx <- match(q$targets$sample, mm$expr_col)
stopifnot(!anyNA(idx[!q$targets$is_pool]))
q$targets$donor <- mm$line_clean[idx]

# drop pooled samples and donors without imaging
keep <- !q$targets$is_pool & q$targets$donor %in% donors$donor

cat("Proteomics but no imaging:",
    setdiff(unique(q$targets$donor[!q$targets$is_pool]), donors$donor), "\n")
cat("Imaging but no proteomics:",
    setdiff(donors$donor, q$targets$donor[!q$targets$is_pool]), "\n")

q_use   <- q[, keep]                       # subsets E, targets, weights, dpc
samples <- q_use$targets

samples$org_sd_deg <- donors$org_sd_deg[match(samples$donor, donors$donor)]

stopifnot(
  ncol(q_use$E) == nrow(samples),
  identical(colnames(q_use$E), samples$sample),      # same order
  nrow(samples) == EXP_SAMPLES,
  length(unique(samples$donor)) == EXP_PROT_DONORS,
  !anyNA(samples$org_sd_deg)
)
cat("STEP 4 OK:", nrow(samples), "samples,",
    length(unique(samples$donor)), "donors,",
    nrow(q_use$E), "proteins\n")

# STEP 5. Design matrices
# median split over donors, not over MS runs
donor_level <- unique(samples[, c("donor", "org_sd_deg")])
cut_deg     <- median(donor_level$org_sd_deg)
cat(sprintf("Median split at %.3f deg over %d donors\n",
            cut_deg, nrow(donor_level)))

samples$org_group <- factor(
  ifelse(samples$org_sd_deg <= cut_deg, "organized", "disorganized"),
  levels = c("disorganized", "organized"))

samples$batch <- factor(samples$batch)
samples$donor <- factor(samples$donor)

# both designs adjust for batch
design_cont <- model.matrix(~ org_sd_deg + batch, data = samples)  # primary
# no intercept, so 02_de.R can name both levels in makeContrasts()
design_cat  <- model.matrix(~ 0 + org_group + batch, data = samples)  # secondary
colnames(design_cat) <- sub("^org_group", "", colnames(design_cat))

stopifnot(
  nrow(design_cont) == ncol(q_use$E),           # one row per sample
  qr(design_cont)$rank == ncol(design_cont),    # no redundant column
  qr(design_cat)$rank  == ncol(design_cat),
  sum(donor_level$org_sd_deg <= cut_deg) == EXP_SPLIT   # even split
)
cat("STEP 5 OK: design_cont", dim(design_cont),
    "| design_cat", dim(design_cat), "\n")

# STEP 6. Save everything 02_de.R needs
dir.create(RUN_DIR, showWarnings = FALSE, recursive = TRUE)

write.csv(donors,  file.path(RUN_DIR, "donor_organization.csv"), row.names = FALSE)
write.csv(samples, file.path(RUN_DIR, "design_samples.csv"),     row.names = FALSE)

# observed organization span (deg); 02_de.R multiplies the per-degree logFC by it
org_range <- diff(range(donor_level$org_sd_deg))

saveRDS(list(mat_quant_use = q_use,
             samples       = samples,
             design_cont   = design_cont,
             design_cat    = design_cat,
             cut_deg       = cut_deg,
             org_range     = org_range),
        file.path(RUN_DIR, "design.rds"))

cat("Saved to", RUN_DIR, "\n")
