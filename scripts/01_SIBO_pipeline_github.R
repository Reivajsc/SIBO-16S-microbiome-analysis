# ==========================================================
# TFM SIBO - Pipeline metataxonómico 16S rRNA
# BioProject: PRJNA1273190
# Cohorte analizada: 123 muestras
#
# Ejecución:
#   1. Situar el directorio de trabajo en la raíz del proyecto.
#   2. Comprobar que existen las carpetas de entrada descritas en la sección 1.
#   3. Ejecutar el script completo desde una sesión limpia de R.
#
# Este script integra:
#   - DADA2
#   - metadatos clínicos (SIBO + hipotiroidismo)
#   - diversidad alfa
#   - diversidad beta (Bray-Curtis + PCoA)
#   - PERMDISP
#   - PERMANOVA sin ajustar y ajustada por hipotiroidismo
#   - composición taxonómica a nivel de género
#   - ANCOM-BC2 sin ajustar y ajustado por hipotiroidismo
#   - tablas y figuras definitivas
#   - sessionInfo()
#
# IMPORTANTE:
# El pipeline reproduce el análisis presentado en el TFM y mantiene los
# parámetros de procesamiento y análisis empleados para los resultados finales.
# ==========================================================


# ==========================================================
# 0. PAQUETES
# ==========================================================

required_packages <- c(
  "dada2",
  "vegan",
  "ggplot2",
  "phyloseq",
  "ANCOMBC",
  "microbiome"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Faltan los siguientes paquetes: ",
      paste(missing_packages, collapse = ", "),
      ". Instálelos antes de ejecutar el pipeline."
    )
  )
}

library(dada2)
library(vegan)
library(ggplot2)
library(phyloseq)
library(ANCOMBC)
library(microbiome)


# ==========================================================
# 1. DIRECTORIOS DEL PROYECTO
# ==========================================================

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)

# FASTQ
path_pos <- file.path(project_dir, "data", "fastq", "SIBO_POS")
path_neg <- file.path(project_dir, "data", "fastq", "SIBO_NEG")
path_control <- file.path(project_dir, "data", "fastq", "CONTROL")

# Metadatos
metadata_file <- file.path(
  project_dir,
  "data",
  "metadata",
  "metadata_123_samples.csv"
)

# Bases de datos SILVA
silva_train <- file.path(
  project_dir,
  "databases",
  "SILVA",
  "silva_nr99_v138.1_train_set.fa.gz"
)

silva_species <- file.path(
  project_dir,
  "databases",
  "SILVA",
  "silva_species_assignment_v138.fa.gz"
)

# Datos procesados
filt_path <- file.path(
  project_dir,
  "data",
  "processed",
  "filtered"
)

# Resultados definitivos
results_dir <- file.path(project_dir, "results")
results_tables <- file.path(results_dir, "tables")
results_figures <- file.path(results_dir, "figures")
results_rds <- file.path(results_dir, "rds")

dirs_to_create <- c(
  filt_path,
  results_dir,
  results_tables,
  results_figures,
  results_rds
)

invisible(
  lapply(
    dirs_to_create,
    dir.create,
    showWarnings = FALSE,
    recursive = TRUE
  )
)

# Comprobar que las entradas necesarias están disponibles
input_dirs <- c(
  path_pos,
  path_neg,
  path_control,
  dirname(metadata_file),
  dirname(silva_train)
)

missing_dirs <- input_dirs[!dir.exists(input_dirs)]

if (length(missing_dirs) > 0) {
  stop(
    paste0(
      "No se encontraron los siguientes directorios de entrada:\n",
      paste(missing_dirs, collapse = "\n")
    )
  )
}

input_files <- c(
  metadata_file,
  silva_train,
  silva_species
)

missing_files <- input_files[!file.exists(input_files)]

if (length(missing_files) > 0) {
  stop(
    paste0(
      "No se encontraron los siguientes archivos de entrada:\n",
      paste(missing_files, collapse = "\n")
    )
  )
}


# ==========================================================
# 2. IMPORTACIÓN Y COMPROBACIÓN DE FASTQ
# ==========================================================

get_fastq_pairs <- function(path) {

  fwd <- sort(
    list.files(
      path,
      pattern = "_1\\.fastq$",
      full.names = TRUE
    )
  )

  rev <- sort(
    list.files(
      path,
      pattern = "_2\\.fastq$",
      full.names = TRUE
    )
  )

  if (length(fwd) != length(rev)) {
    stop(
      paste0(
        "Número diferente de forward y reverse en: ",
        path
      )
    )
  }

  list(F = fwd, R = rev)
}

fastq_pos <- get_fastq_pairs(path_pos)
fastq_neg <- get_fastq_pairs(path_neg)
fastq_control <- get_fastq_pairs(path_control)

fnFs_pos <- fastq_pos$F
fnRs_pos <- fastq_pos$R

fnFs_neg <- fastq_neg$F
fnRs_neg <- fastq_neg$R

fnFs_control <- fastq_control$F
fnRs_control <- fastq_control$R

# Comprobación esperada de la cohorte analizada
if (length(fnFs_pos) != 54) {
  warning("No se detectaron exactamente 54 muestras en SIBO_POS.")
}

if (length(fnFs_neg) != 29) {
  warning("No se detectaron exactamente 29 muestras en SIBO_NEG.")
}

if (length(fnFs_control) != 40) {
  warning("No se detectaron exactamente 40 muestras en CONTROL.")
}

fnFs <- c(fnFs_pos, fnFs_neg, fnFs_control)
fnRs <- c(fnRs_pos, fnRs_neg, fnRs_control)

sample.names <- sub(
  "_1\\.fastq$",
  "",
  basename(fnFs)
)

if (anyDuplicated(sample.names)) {
  stop("Existen identificadores de muestra duplicados.")
}

if (length(fnFs) != 123) {
  warning(
    paste0(
      "Se detectaron ",
      length(fnFs),
      " muestras, no 123."
    )
  )
}

cat("\n============================================\n")
cat("FASTQ DETECTADOS\n")
cat("============================================\n")
cat("SIBO_POS: ", length(fnFs_pos), "\n", sep = "")
cat("SIBO_NEG: ", length(fnFs_neg), "\n", sep = "")
cat("CONTROL: ", length(fnFs_control), "\n", sep = "")
cat("TOTAL: ", length(fnFs), "\n", sep = "")


# ==========================================================
# 3. PERFILES DE CALIDAD PREVIOS AL FILTRADO
# ==========================================================

# Se muestran muestras representativas de los tres conjuntos.
quality_F <- c(
  fnFs_pos[1],
  fnFs_neg[1],
  fnFs_control[1]
)

quality_R <- c(
  fnRs_pos[1],
  fnRs_neg[1],
  fnRs_control[1]
)

png(
  file.path(results_figures, "quality_profile_forward.png"),
  width = 2400,
  height = 1600,
  res = 200
)
plotQualityProfile(quality_F)
dev.off()

png(
  file.path(results_figures, "quality_profile_reverse.png"),
  width = 2400,
  height = 1600,
  res = 200
)
plotQualityProfile(quality_R)
dev.off()


# ==========================================================
# 4. FILTRADO DE CALIDAD MEDIANTE DADA2
# ==========================================================

filtFs <- file.path(
  filt_path,
  paste0(sample.names, "_F_filt.fastq.gz")
)

filtRs <- file.path(
  filt_path,
  paste0(sample.names, "_R_filt.fastq.gz")
)

out <- filterAndTrim(
  fnFs,
  filtFs,
  fnRs,
  filtRs,
  truncLen = c(280, 220),
  maxN = 0,
  maxEE = c(2, 2),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = FALSE
)

write.csv(
  out,
  file.path(results_tables, "filtering_summary.csv"),
  row.names = TRUE
)


# ==========================================================
# 5. ESTIMACIÓN DE LOS MODELOS DE ERROR
# ==========================================================

errF <- learnErrors(
  filtFs,
  multithread = FALSE
)

errR <- learnErrors(
  filtRs,
  multithread = FALSE
)

png(
  file.path(results_figures, "dada2_error_model_forward.png"),
  width = 3600,
  height = 2400,
  res = 300
)
plotErrors(errF, nominalQ = TRUE)
dev.off()

png(
  file.path(results_figures, "dada2_error_model_reverse.png"),
  width = 3600,
  height = 2400,
  res = 300
)
plotErrors(errR, nominalQ = TRUE)
dev.off()


# ==========================================================
# 6. INFERENCIA DE ASVs
# ==========================================================

dadaFs <- dada(
  filtFs,
  err = errF,
  multithread = FALSE
)

dadaRs <- dada(
  filtRs,
  err = errR,
  multithread = FALSE
)


# ==========================================================
# 7. ENSAMBLAJE FORWARD + REVERSE
# ==========================================================

mergers <- mergePairs(
  dadaFs,
  filtFs,
  dadaRs,
  filtRs,
  verbose = TRUE
)

seqtab <- makeSequenceTable(mergers)


# ==========================================================
# 8. ELIMINACIÓN DE QUIMERAS
# ==========================================================

seqtab.nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = FALSE,
  verbose = TRUE
)

chimera_summary <- data.frame(
  ASVs_initial = ncol(seqtab),
  ASVs_chimeric = ncol(seqtab) - ncol(seqtab.nochim),
  ASVs_final = ncol(seqtab.nochim),
  reads_retained_percent =
    100 * sum(seqtab.nochim) / sum(seqtab)
)

write.csv(
  chimera_summary,
  file.path(results_tables, "chimera_removal_summary.csv"),
  row.names = FALSE
)


# ==========================================================
# 8.1. TRACKING DE LECTURAS
# ==========================================================

getN <- function(x) {
  sum(getUniques(x))
}

track <- cbind(
  input = out[, 1],
  filtered = out[, 2],
  denoisedF = sapply(dadaFs, getN),
  denoisedR = sapply(dadaRs, getN),
  merged = sapply(mergers, getN),
  nonchim = rowSums(seqtab.nochim)
)

rownames(track) <- sample.names

write.csv(
  track,
  file.path(results_tables, "read_tracking_by_sample.csv"),
  row.names = TRUE
)

tracking_total <- data.frame(
  Stage = c(
    "Input",
    "Filtered",
    "Denoised_forward",
    "Denoised_reverse",
    "Merged",
    "Non_chimeric"
  ),
  Reads = c(
    sum(track[, "input"]),
    sum(track[, "filtered"]),
    sum(track[, "denoisedF"]),
    sum(track[, "denoisedR"]),
    sum(track[, "merged"]),
    sum(track[, "nonchim"])
  )
)

tracking_total$Percent_input <-
  100 * tracking_total$Reads / tracking_total$Reads[1]

write.csv(
  tracking_total,
  file.path(results_tables, "read_tracking_total.csv"),
  row.names = FALSE
)


# ==========================================================
# 9. ASIGNACIÓN TAXONÓMICA CON SILVA v138.1
# ==========================================================

taxa <- assignTaxonomy(
  seqtab.nochim,
  silva_train,
  multithread = FALSE
)

taxa <- addSpecies(
  taxa,
  silva_species
)


# ==========================================================
# 10. NORMALIZACIÓN DE LOS IDENTIFICADORES DE MUESTRA
# ==========================================================

sample_ids <- rownames(seqtab.nochim)

sample_ids <- sub(
  "_F_filt\\.fastq\\.gz$",
  "",
  sample_ids
)

rownames(seqtab.nochim) <- sample_ids


# ==========================================================
# 11. METADATOS CLÍNICOS
#     SIBO + HIPOTIROIDISMO DESDE EL INICIO DEL ANÁLISIS
# ==========================================================

metadata_sra <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_metadata_cols <- c(
  "Run",
  "Group_hypothyroidism"
)

if (!all(required_metadata_cols %in% colnames(metadata_sra))) {
  stop(
    "SraRunTable.csv no contiene las columnas Run y/o Group_hypothyroidism."
  )
}

idx <- match(
  rownames(seqtab.nochim),
  metadata_sra$Run
)

if (any(is.na(idx))) {
  stop(
    paste(
      "Muestras sin metadatos:",
      paste(
        rownames(seqtab.nochim)[is.na(idx)],
        collapse = ", "
      )
    )
  )
}

metadata <- data.frame(
  sample_id = rownames(seqtab.nochim),
  original_group =
    metadata_sra$Group_hypothyroidism[idx],
  stringsAsFactors = FALSE
)

metadata$SIBO <- ifelse(
  metadata$original_group %in%
    c(
      "SIBO_noHypothyroidism",
      "SIBO_Hypothyroidism"
    ),
  "SIBO_POS",
  "SIBO_NEG"
)

metadata$hypothyroidism <- ifelse(
  metadata$original_group %in%
    c(
      "noSIBO_Hypothyroidism",
      "SIBO_Hypothyroidism"
    ),
  "Hypo_POS",
  "Hypo_NEG"
)

valid_groups <- c(
  "Control",
  "SIBO_noHypothyroidism",
  "noSIBO_Hypothyroidism",
  "SIBO_Hypothyroidism"
)

if (!all(metadata$original_group %in% valid_groups)) {
  stop("Se detectaron categorías clínicas no reconocidas.")
}

metadata$SIBO <- factor(
  metadata$SIBO,
  levels = c("SIBO_NEG", "SIBO_POS")
)

metadata$hypothyroidism <- factor(
  metadata$hypothyroidism,
  levels = c("Hypo_NEG", "Hypo_POS")
)

rownames(metadata) <- metadata$sample_id

stopifnot(
  identical(
    rownames(seqtab.nochim),
    rownames(metadata)
  )
)

subgroup_counts <- table(
  SIBO = metadata$SIBO,
  Hypothyroidism = metadata$hypothyroidism
)

cat("\n============================================\n")
cat("DISTRIBUCIÓN FINAL DE LA COHORTE\n")
cat("============================================\n")
print(subgroup_counts)

write.csv(
  metadata,
  file.path(results_tables, "metadata_complete.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(subgroup_counts),
  file.path(results_tables, "sample_distribution.csv"),
  row.names = FALSE
)


# ==========================================================
# 12. GUARDADO DE MATRICES PRINCIPALES
# ==========================================================

saveRDS(
  seqtab.nochim,
  file.path(results_rds, "seqtab_nochim.rds")
)

saveRDS(
  taxa,
  file.path(results_rds, "taxa_silva.rds")
)

write.csv(
  seqtab.nochim,
  file.path(results_tables, "ASV_table.csv"),
  row.names = TRUE
)

write.csv(
  taxa,
  file.path(results_tables, "taxonomy_silva.csv"),
  row.names = TRUE
)


# ==========================================================
# 13. PROFUNDIDAD DE SECUENCIACIÓN
# ==========================================================

sequencing_depth <- data.frame(
  sample = rownames(seqtab.nochim),
  SIBO = metadata$SIBO,
  hypothyroidism = metadata$hypothyroidism,
  reads = rowSums(seqtab.nochim)
)

write.csv(
  sequencing_depth,
  file.path(results_tables, "sequencing_depth_by_sample.csv"),
  row.names = FALSE
)

depth_summary <- do.call(
  rbind,
  lapply(
    split(
      sequencing_depth$reads,
      sequencing_depth$SIBO
    ),
    function(x) {
      data.frame(
        median = median(x),
        Q1 = unname(quantile(x, 0.25)),
        Q3 = unname(quantile(x, 0.75)),
        min = min(x),
        max = max(x)
      )
    }
  )
)

depth_summary$SIBO <- rownames(depth_summary)
rownames(depth_summary) <- NULL

write.csv(
  depth_summary,
  file.path(
    results_tables,
    "sequencing_depth_summary_by_SIBO.csv"
  ),
  row.names = FALSE
)


# ==========================================================
# 14. OBJETO PHYLOSEQ
# ==========================================================

OTU <- otu_table(
  seqtab.nochim,
  taxa_are_rows = FALSE
)

TAX <- tax_table(
  taxa
)

META <- sample_data(
  metadata
)

physeq <- phyloseq(
  OTU,
  TAX,
  META
)

saveRDS(
  physeq,
  file.path(results_rds, "physeq_ASV.rds")
)

cat(
  "\nObjeto phyloseq: ",
  nsamples(physeq),
  " muestras / ",
  ntaxa(physeq),
  " ASVs\n",
  sep = ""
)


# ==========================================================
# 15. DIVERSIDAD ALFA
# ==========================================================

observed_asvs <- rowSums(seqtab.nochim > 0)

shannon <- vegan::diversity(
  seqtab.nochim,
  index = "shannon"
)

simpson <- vegan::diversity(
  seqtab.nochim,
  index = "simpson"
)

alpha_div <- data.frame(
  sample = rownames(seqtab.nochim),
  SIBO = metadata$SIBO,
  hypothyroidism = metadata$hypothyroidism,
  Observed_ASVs = observed_asvs,
  Shannon = shannon,
  Simpson = simpson
)

write.csv(
  alpha_div,
  file.path(results_tables, "alpha_diversity_by_sample.csv"),
  row.names = FALSE
)

alpha_summary_fun <- function(x) {
  paste0(
    round(median(x), 2),
    " (",
    round(unname(quantile(x, 0.25)), 2),
    "–",
    round(unname(quantile(x, 0.75)), 2),
    ")"
  )
}

alpha_summary <- aggregate(
  cbind(
    Observed_ASVs,
    Shannon,
    Simpson
  ) ~ SIBO,
  data = alpha_div,
  FUN = alpha_summary_fun
)

wilcox_observed <- wilcox.test(
  Observed_ASVs ~ SIBO,
  data = alpha_div
)

wilcox_shannon <- wilcox.test(
  Shannon ~ SIBO,
  data = alpha_div
)

wilcox_simpson <- wilcox.test(
  Simpson ~ SIBO,
  data = alpha_div
)

alpha_pvalues <- data.frame(
  Metric = c(
    "Observed_ASVs",
    "Shannon",
    "Simpson"
  ),
  p_value = c(
    wilcox_observed$p.value,
    wilcox_shannon$p.value,
    wilcox_simpson$p.value
  )
)

write.csv(
  alpha_summary,
  file.path(results_tables, "alpha_diversity_median_IQR.csv"),
  row.names = FALSE
)

write.csv(
  alpha_pvalues,
  file.path(results_tables, "alpha_diversity_pvalues.csv"),
  row.names = FALSE
)


# ==========================================================
# 15.1. FIGURA DEFINITIVA DE DIVERSIDAD ALFA
# ==========================================================

alpha_long <- rbind(
  data.frame(
    SIBO = alpha_div$SIBO,
    Metric = "Riqueza observada",
    Value = alpha_div$Observed_ASVs
  ),
  data.frame(
    SIBO = alpha_div$SIBO,
    Metric = "Shannon",
    Value = alpha_div$Shannon
  ),
  data.frame(
    SIBO = alpha_div$SIBO,
    Metric = "Simpson",
    Value = alpha_div$Simpson
  )
)

alpha_long$Metric <- factor(
  alpha_long$Metric,
  levels = c(
    "Riqueza observada",
    "Shannon",
    "Simpson"
  )
)

alpha_long$SIBO <- factor(
  alpha_long$SIBO,
  levels = c("SIBO_NEG", "SIBO_POS"),
  labels = c("SIBO−", "SIBO+")
)

p_labels <- data.frame(
  Metric = factor(
    c(
      "Riqueza observada",
      "Shannon",
      "Simpson"
    ),
    levels = c(
      "Riqueza observada",
      "Shannon",
      "Simpson"
    )
  ),
  label = paste0(
    "p = ",
    formatC(
      alpha_pvalues$p_value,
      format = "f",
      digits = 3,
      decimal.mark = ","
    )
  ),
  x = 1.5,
  y = c(
    max(alpha_div$Observed_ASVs, na.rm = TRUE) * 1.08,
    max(alpha_div$Shannon, na.rm = TRUE) * 1.08,
    max(alpha_div$Simpson, na.rm = TRUE) * 1.08
  )
)

p_alpha <- ggplot(
  alpha_long,
  aes(x = SIBO, y = Value)
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    size = 1.5,
    alpha = 0.65
  ) +
  facet_wrap(
    ~ Metric,
    scales = "free_y",
    nrow = 1
  ) +
  geom_text(
    data = p_labels,
    aes(
      x = x,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    size = 4.5
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    strip.text = element_text(
      face = "bold",
      size = 14
    ),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 11),
    panel.spacing = grid::unit(1.5, "lines")
  )

ggsave(
  filename = file.path(
    results_figures,
    "alpha_diversity.png"
  ),
  plot = p_alpha,
  width = 11,
  height = 4.8,
  dpi = 300
)


# ==========================================================
# 16. DIVERSIDAD BETA
#     BRAY-CURTIS SOBRE ABUNDANCIAS RELATIVAS
# ==========================================================

asv_matrix <- as.matrix(seqtab.nochim)

storage.mode(asv_matrix) <- "numeric"

asv_relative <- asv_matrix /
  rowSums(asv_matrix)

bray_dist <- vegan::vegdist(
  asv_relative,
  method = "bray"
)


# ==========================================================
# 16.1. PCoA
# ==========================================================

pcoa_calc <- cmdscale(
  bray_dist,
  k = 2,
  eig = TRUE
)

positive_eig <- pcoa_calc$eig[
  pcoa_calc$eig > 0
]

PC1_percent <- round(
  100 * pcoa_calc$eig[1] /
    sum(positive_eig),
  1
)

PC2_percent <- round(
  100 * pcoa_calc$eig[2] /
    sum(positive_eig),
  1
)

pcoa_df <- data.frame(
  Sample = rownames(asv_relative),
  SIBO = metadata$SIBO,
  hypothyroidism = metadata$hypothyroidism,
  PC1 = pcoa_calc$points[, 1],
  PC2 = pcoa_calc$points[, 2]
)

write.csv(
  pcoa_df,
  file.path(
    results_tables,
    "pcoa_bray_curtis_coordinates.csv"
  ),
  row.names = FALSE
)

pcoa_df$SIBO_plot <- factor(
  pcoa_df$SIBO,
  levels = c("SIBO_NEG", "SIBO_POS"),
  labels = c("SIBO−", "SIBO+")
)

p_pcoa <- ggplot(
  pcoa_df,
  aes(
    x = PC1,
    y = PC2,
    colour = SIBO_plot,
    fill = SIBO_plot
  )
) +
  stat_ellipse(
    geom = "polygon",
    alpha = 0.12,
    linewidth = 1,
    type = "t"
  ) +
  geom_point(
    size = 2.5,
    alpha = 0.8
  ) +
  labs(
    x = paste0("PC1 (", PC1_percent, " %)"),
    y = paste0("PC2 (", PC2_percent, " %)"),
    colour = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "bottom",
    axis.title = element_text(
      size = 14,
      face = "bold"
    )
  )

ggsave(
  filename = file.path(
    results_figures,
    "PCoA_BrayCurtis.png"
  ),
  plot = p_pcoa,
  width = 8.5,
  height = 6.5,
  dpi = 300
)


# ==========================================================
# 16.2. PERMDISP
# ==========================================================

sibo_group <- factor(
  metadata$SIBO
)

dispersion_sibo <- vegan::betadisper(
  bray_dist,
  sibo_group
)

set.seed(123)

permdisp_sibo <- vegan::permutest(
  dispersion_sibo,
  permutations = 999
)

permdisp_table <- as.data.frame(
  permdisp_sibo$tab
)

write.csv(
  permdisp_table,
  file.path(results_tables, "permdisp_SIBO.csv"),
  row.names = TRUE
)

capture.output(
  permdisp_sibo,
  file = file.path(
    results_tables,
    "permdisp_SIBO.txt"
  )
)


# ==========================================================
# 16.3. PERMANOVA NO AJUSTADA
# ==========================================================

set.seed(123)

permanova_unadjusted <- vegan::adonis2(
  bray_dist ~ SIBO,
  data = metadata,
  permutations = 999
)

write.csv(
  as.data.frame(permanova_unadjusted),
  file.path(
    results_tables,
    "permanova_unadjusted.csv"
  ),
  row.names = TRUE
)

capture.output(
  permanova_unadjusted,
  file = file.path(
    results_tables,
    "permanova_unadjusted.txt"
  )
)


# ==========================================================
# 16.4. PERMANOVA AJUSTADA POR HIPOTIROIDISMO
# ==========================================================

set.seed(123)

permanova_adjusted <- vegan::adonis2(
  bray_dist ~ SIBO + hypothyroidism,
  data = metadata,
  permutations = 999,
  by = "margin"
)

write.csv(
  as.data.frame(permanova_adjusted),
  file.path(
    results_tables,
    "permanova_adjusted_hypothyroidism.csv"
  ),
  row.names = TRUE
)

capture.output(
  permanova_adjusted,
  file = file.path(
    results_tables,
    "permanova_adjusted_hypothyroidism.txt"
  )
)


# ==========================================================
# 17. COMPOSICIÓN TAXONÓMICA A NIVEL DE GÉNERO
# ==========================================================

physeq_genus <- tax_glom(
  physeq,
  taxrank = "Genus",
  NArm = FALSE
)

physeq_genus_rel <- transform_sample_counts(
  physeq_genus,
  function(x) x / sum(x)
)

saveRDS(
  physeq_genus,
  file.path(results_rds, "physeq_genus_counts.rds")
)

saveRDS(
  physeq_genus_rel,
  file.path(results_rds, "physeq_genus_relative.rds")
)

otu_genus <- as(
  otu_table(physeq_genus_rel),
  "matrix"
)

if (taxa_are_rows(physeq_genus_rel)) {
  otu_genus <- t(otu_genus)
}

tax_genus <- as(
  tax_table(physeq_genus_rel),
  "matrix"
)

genus_names <- tax_genus[, "Genus"]

genus_names[
  is.na(genus_names) |
    genus_names == ""
] <- "Unassigned"

stopifnot(
  identical(
    colnames(otu_genus),
    rownames(tax_genus)
  )
)

mean_global <- colMeans(otu_genus)

genus_summary <- data.frame(
  Taxon_ID = names(mean_global),
  Genus = genus_names[names(mean_global)],
  Mean_global = as.numeric(mean_global),
  stringsAsFactors = FALSE
)

genus_summary <- genus_summary[
  order(
    genus_summary$Mean_global,
    decreasing = TRUE
  ),
]

top15 <- head(
  genus_summary,
  15
)

meta_phy <- data.frame(
  sample_data(physeq_genus_rel)
)

meta_phy <- meta_phy[
  rownames(otu_genus),
  ,
  drop = FALSE
]

neg <- meta_phy$SIBO == "SIBO_NEG"
pos <- meta_phy$SIBO == "SIBO_POS"

top15_ids <- top15$Taxon_ID

mean_neg <- colMeans(
  otu_genus[
    neg,
    top15_ids,
    drop = FALSE
  ]
)

mean_pos <- colMeans(
  otu_genus[
    pos,
    top15_ids,
    drop = FALSE
  ]
)

top15_SIBO <- data.frame(
  Genus = top15$Genus,
  SIBO_NEG =
    round(mean_neg[top15_ids] * 100, 2),
  SIBO_POS =
    round(mean_pos[top15_ids] * 100, 2),
  Media_global =
    round(top15$Mean_global * 100, 2),
  row.names = NULL
)

write.csv(
  top15_SIBO,
  file.path(
    results_tables,
    "top15_genus_abundance_by_SIBO.csv"
  ),
  row.names = FALSE
)


# ==========================================================
# 17.1. DOTPLOT TOP 15
# ==========================================================

plot_taxa <- rbind(
  data.frame(
    Genus = top15_SIBO$Genus,
    Group = "SIBO−",
    Abundance = top15_SIBO$SIBO_NEG
  ),
  data.frame(
    Genus = top15_SIBO$Genus,
    Group = "SIBO+",
    Abundance = top15_SIBO$SIBO_POS
  )
)

plot_taxa$Genus <- factor(
  plot_taxa$Genus,
  levels = rev(top15_SIBO$Genus)
)

p_taxa_dot <- ggplot(
  plot_taxa,
  aes(
    x = Abundance,
    y = Genus,
    shape = Group
  )
) +
  geom_line(
    aes(group = Genus),
    linewidth = 0.7
  ) +
  geom_point(
    size = 3.2
  ) +
  labs(
    x = "Abundancia relativa media (%)",
    y = NULL,
    shape = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(
      face = "italic",
      size = 11
    ),
    axis.title.x = element_text(
      size = 13,
      face = "bold"
    )
  )

ggsave(
  filename = file.path(
    results_figures,
    "taxonomic_composition_top15_dotplot.png"
  ),
  plot = p_taxa_dot,
  width = 9,
  height = 7,
  dpi = 300
)


# ==========================================================
# 18. ANCOM-BC2
#     FILTRADO POR PREVALENCIA
# ==========================================================

otu_genus_counts <- as(
  otu_table(physeq_genus),
  "matrix"
)

if (!taxa_are_rows(physeq_genus)) {
  otu_genus_counts <- t(otu_genus_counts)
}

prevalence <- rowSums(
  otu_genus_counts > 0
) / ncol(otu_genus_counts)

prevalence_summary <- data.frame(
  Threshold = c(
    ">=5%",
    ">=10%",
    ">=20%",
    ">=25%"
  ),
  Taxa = c(
    sum(prevalence >= 0.05),
    sum(prevalence >= 0.10),
    sum(prevalence >= 0.20),
    sum(prevalence >= 0.25)
  )
)

write.csv(
  prevalence_summary,
  file.path(
    results_tables,
    "ANCOMBC2_prevalence_summary.csv"
  ),
  row.names = FALSE
)


# ==========================================================
# 18.1. ANCOM-BC2 NO AJUSTADO
# ==========================================================

sample_data(physeq_genus)$SIBO <- relevel(
  factor(sample_data(physeq_genus)$SIBO),
  ref = "SIBO_NEG"
)

set.seed(123)

ancom_sibo <- ancombc2(
  data = physeq_genus,
  tax_level = "Genus",
  fix_formula = "SIBO",
  p_adj_method = "BH",
  prv_cut = 0.10,
  lib_cut = 0,
  group = "SIBO",
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = 0.05,
  global = FALSE,
  pairwise = FALSE,
  dunnet = FALSE,
  trend = FALSE,
  iter_control = list(
    tol = 1e-5,
    max_iter = 100,
    verbose = TRUE
  ),
  em_control = list(
    tol = 1e-5,
    max_iter = 100
  ),
  lme_control = NULL,
  mdfdr_control = list(
    fwer_ctrl_method = "holm",
    B = 100
  ),
  verbose = TRUE
)

res_ancom <- ancom_sibo$res

res_sibo <- data.frame(
  Genus = res_ancom$taxon,
  logFC = res_ancom$lfc_SIBOSIBO_POS,
  SE = res_ancom$se_SIBOSIBO_POS,
  W = res_ancom$W_SIBOSIBO_POS,
  p_value = res_ancom$p_SIBOSIBO_POS,
  FDR = res_ancom$q_SIBOSIBO_POS,
  Differential =
    res_ancom$diff_SIBOSIBO_POS,
  Passed_SS =
    res_ancom$passed_ss_SIBOSIBO_POS,
  Robust =
    res_ancom$diff_robust_SIBOSIBO_POS,
  stringsAsFactors = FALSE
)

res_sibo <- res_sibo[
  order(res_sibo$FDR),
]

write.csv(
  res_sibo,
  file.path(
    results_tables,
    "ANCOMBC2_SIBO_unadjusted.csv"
  ),
  row.names = FALSE
)


# ==========================================================
# 18.2. ANCOM-BC2 AJUSTADO POR HIPOTIROIDISMO
# ==========================================================

sample_data(physeq_genus)$SIBO <- relevel(
  factor(sample_data(physeq_genus)$SIBO),
  ref = "SIBO_NEG"
)

sample_data(physeq_genus)$hypothyroidism <- relevel(
  factor(
    sample_data(physeq_genus)$hypothyroidism
  ),
  ref = "Hypo_NEG"
)

set.seed(123)

ancom_sibo_adjusted <- ancombc2(
  data = physeq_genus,
  tax_level = "Genus",
  fix_formula = "SIBO + hypothyroidism",
  p_adj_method = "BH",
  prv_cut = 0.10,
  lib_cut = 0,
  group = "SIBO",
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = 0.05,
  global = FALSE,
  pairwise = FALSE,
  dunnet = FALSE,
  trend = FALSE,
  iter_control = list(
    tol = 1e-5,
    max_iter = 100,
    verbose = TRUE
  ),
  em_control = list(
    tol = 1e-5,
    max_iter = 100
  ),
  lme_control = NULL,
  mdfdr_control = list(
    fwer_ctrl_method = "holm",
    B = 100
  ),
  verbose = TRUE
)

res_adj <- ancom_sibo_adjusted$res

res_sibo_adj <- data.frame(
  Genus = res_adj$taxon,
  logFC = res_adj$lfc_SIBOSIBO_POS,
  SE = res_adj$se_SIBOSIBO_POS,
  W = res_adj$W_SIBOSIBO_POS,
  p_value = res_adj$p_SIBOSIBO_POS,
  FDR = res_adj$q_SIBOSIBO_POS,
  Differential =
    res_adj$diff_SIBOSIBO_POS,
  Passed_SS =
    res_adj$passed_ss_SIBOSIBO_POS,
  Robust =
    res_adj$diff_robust_SIBOSIBO_POS,
  stringsAsFactors = FALSE
)

res_sibo_adj <- res_sibo_adj[
  order(res_sibo_adj$FDR),
]

write.csv(
  res_sibo_adj,
  file.path(
    results_tables,
    "ANCOMBC2_SIBO_adjusted_hypothyroidism.csv"
  ),
  row.names = FALSE
)


# ==========================================================
# 18.3. COMPARACIÓN DE MODELOS ANCOM-BC2
# ==========================================================

sig_unadj <- res_sibo$Genus[
  res_sibo$FDR < 0.05
]

sig_adj <- res_sibo_adj$Genus[
  res_sibo_adj$FDR < 0.05
]

comparison_ancom <- merge(
  res_sibo[
    ,
    c("Genus", "logFC", "FDR")
  ],
  res_sibo_adj[
    ,
    c("Genus", "logFC", "FDR")
  ],
  by = "Genus",
  suffixes = c(
    "_unadjusted",
    "_adjusted"
  )
)

comparison_sig <- comparison_ancom[
  comparison_ancom$FDR_unadjusted < 0.05 |
    comparison_ancom$FDR_adjusted < 0.05,
]

comparison_sig <- comparison_sig[
  order(comparison_sig$FDR_adjusted),
]

write.csv(
  comparison_sig,
  file.path(
    results_tables,
    "ANCOMBC2_comparison_unadjusted_adjusted.csv"
  ),
  row.names = FALSE
)


# ==========================================================
# 18.4. TABLA FINAL ANCOM-BC2
#      MODELO AJUSTADO POR HIPOTIROIDISMO
# ==========================================================

tabla_ancom <- res_sibo_adj[
  res_sibo_adj$FDR < 0.05,
  c(
    "Genus",
    "logFC",
    "SE",
    "FDR"
  )
]

tabla_ancom$IC95_inf <-
  tabla_ancom$logFC -
  1.96 * tabla_ancom$SE

tabla_ancom$IC95_sup <-
  tabla_ancom$logFC +
  1.96 * tabla_ancom$SE

tabla_ancom <- tabla_ancom[
  order(tabla_ancom$FDR),
]

tabla_final <- data.frame(
  Genero = tabla_ancom$Genus,
  logFC = round(tabla_ancom$logFC, 2),
  IC95 = paste0(
    round(tabla_ancom$IC95_inf, 2),
    " a ",
    round(tabla_ancom$IC95_sup, 2)
  ),
  FDR = format(
    tabla_ancom$FDR,
    scientific = TRUE,
    digits = 3
  ),
  check.names = FALSE
)

write.csv(
  tabla_final,
  file.path(
    results_tables,
    "ANCOMBC2_differential_genera_adjusted.csv"
  ),
  row.names = FALSE
)


# ==========================================================
# 18.5. FOREST PLOT ANCOM-BC2 AJUSTADO
# ==========================================================

forest_data <- res_sibo_adj[
  res_sibo_adj$FDR < 0.05,
]

forest_data$CI_low <-
  forest_data$logFC -
  1.96 * forest_data$SE

forest_data$CI_high <-
  forest_data$logFC +
  1.96 * forest_data$SE

forest_data <- forest_data[
  order(forest_data$logFC),
]

forest_data$Genus <- factor(
  forest_data$Genus,
  levels = forest_data$Genus
)

p_forest <- ggplot(
  forest_data,
  aes(
    x = logFC,
    y = Genus
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_errorbar(
    aes(
      xmin = CI_low,
      xmax = CI_high
    ),
    orientation = "y",
    width = 0.18,
    linewidth = 0.7
  ) +
  geom_point(
    size = 3
  ) +
  labs(
    x = "logFC (SIBO+ vs SIBO−)",
    y = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.y =
      element_text(face = "italic"),
    axis.title.x =
      element_text(face = "bold"),
    plot.margin =
      margin(10, 20, 10, 10)
  )

ggsave(
  filename = file.path(
    results_figures,
    "ANCOMBC2_forest_SIBO_adjusted.png"
  ),
  plot = p_forest,
  width = 9,
  height = 7,
  units = "in",
  dpi = 600
)


# ==========================================================
# 19. RESUMEN DE RESULTADOS EN CONSOLA
# ==========================================================

cat("\n============================================\n")
cat("RESUMEN DEL PIPELINE\n")
cat("============================================\n")

cat(
  "Muestras finales: ",
  nrow(seqtab.nochim),
  "\n",
  sep = ""
)

cat(
  "ASVs finales: ",
  ncol(seqtab.nochim),
  "\n",
  sep = ""
)

cat("\nDistribución clínica:\n")
print(subgroup_counts)

cat("\nDiversidad alfa - p-valores Wilcoxon:\n")
print(alpha_pvalues)

cat("\nPERMDISP:\n")
print(permdisp_sibo)

cat("\nPERMANOVA no ajustada:\n")
print(permanova_unadjusted)

cat("\nPERMANOVA ajustada por hipotiroidismo:\n")
print(permanova_adjusted)

cat(
  "\nANCOM-BC2 significativos sin ajustar (FDR < 0.05): ",
  sum(res_sibo$FDR < 0.05, na.rm = TRUE),
  "\n",
  sep = ""
)

cat(
  "ANCOM-BC2 significativos ajustados (FDR < 0.05): ",
  sum(res_sibo_adj$FDR < 0.05, na.rm = TRUE),
  "\n",
  sep = ""
)

cat(
  "ANCOM-BC2 robustos ajustados: ",
  sum(res_sibo_adj$Robust, na.rm = TRUE),
  "\n",
  sep = ""
)

cat("\nSignificativos que dejan de serlo al ajustar:\n")
print(setdiff(sig_unadj, sig_adj))

cat("\nNuevos significativos tras el ajuste:\n")
print(setdiff(sig_adj, sig_unadj))

cat("\nSignificativos comunes a ambos modelos:\n")
print(intersect(sig_unadj, sig_adj))


# ==========================================================
# 20. REPRODUCIBILIDAD
# ==========================================================

capture.output(
  sessionInfo(),
  file = file.path(
    project_dir,
    "sessionInfo.txt"
  )
)

package_versions <- data.frame(
  Package = required_packages,
  Version = vapply(
    required_packages,
    function(pkg) {
      as.character(packageVersion(pkg))
    },
    character(1)
  )
)

write.csv(
  package_versions,
  file.path(
    results_tables,
    "package_versions.csv"
  ),
  row.names = FALSE
)

cat("\n============================================\n")
cat("PIPELINE FINALIZADO\n")
cat("Resultados guardados en:\n")
cat(results_dir, "\n")
cat("============================================\n")
