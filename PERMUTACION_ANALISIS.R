
# --------------------------------------------------

# RNA-seq Analysis in R

# Date: 07-09-2026

# Author: Melannie Cedeño Valencia


# --------------------------------------------------

#1 Set work directory

PERM_DIR <- "C:/Users/mcede/Documents/TFM_2026/R_GENE_EXPRESSION/PERMUTACION_ANALISIS"

SALMON_DIR <- "C:/Users/mcede/Documents/TFM_2026/salmon_results"

REFERENCE_DIR <- "C:/Users/mcede/Documents/TFM_2026/reference"

setwd(PERM_DIR)

##2 Llamado de librerias

library(tximport)
library(DESeq2)
library(edgeR)
library(rtracklayer)
library(ggplot2)
library(tidyverse)
library(pheatmap)
library(dplyr)
library(openxlsx)
library(ggrepel)
library(VennDiagram)
library(grid)

#########################

##CARGA DE ANOTACION

#########################
# 3 Archivo GTF

gtf_file <- file.path(
  REFERENCE_DIR,
  "gencode.v49.annotation.gtf.gz"
)

gtf <- import(gtf_file)

tx2gene_p <- unique(
  data.frame(
    transcript_id = mcols(gtf)$transcript_id,
    gene_id = mcols(gtf)$gene_id,
    gene_name = mcols(gtf)$gene_name,
    gene_type = mcols(gtf)$gene_type
  )
)

tx2gene_p <- na.omit(tx2gene_p)

gene_annotation_p <- unique(
  tx2gene_p[, c("gene_id", "gene_name", "gene_type")]
)

cat(
  "Numero de transcritos en tx2gene:",
  nrow(tx2gene_p),
  "\n"
)


head(tx2gene_p)

#########################

#Tabla de anotaciones

#########################

gene_annotation_p <- unique(
  tx2gene_p[, c("gene_id", "gene_name", "gene_type")]
)

head(gene_annotation_p)

#########################

#4 DEFINICION DE MUESTRAS REALES

#########################
samples_p <- data.frame(
  sample = c(
    "SRR6006895",
    "SRR6006896",
    "SRR6006898",
    "SRR6006900",
    "SRR6006904",
    "SRR6006906"
  ),
  condition_real = c(
    "Control",
    "Control",
    "Control",
    "JIA",
    "JIA",
    "JIA"
  )
)
############################
# 5. ARCHIVOS QUANT.SF DE SALMON
############################

files_p <- file.path(
  SALMON_DIR,
  c(
    "SRR6006895_quant",
    "SRR6006896_quant",
    "SRR6006898_quant",
    "SRR6006900_quant",
    "SRR6006904_quant",
    "SRR6006906_quant"
  ),
  "quant.sf"
)

names(files_p) <- samples_p$sample


############################
# 6. IMPORTACION CON TXIMPORT
############################

txi_p <- tximport(
  files_p,
  type = "salmon",
  tx2gene = tx2gene_p[, c("transcript_id", "gene_id")]
)

dim(txi_p$counts)

############################
# 7. ORDENAR Y VERIFICAR METADATOS
############################

samples_p <- samples_p[
  match(colnames(txi_p$counts), samples_p$sample),
]

print(samples_p)


############################
# 8. CREAR PERMUTACIONES 3 vs 3
############################

# Vector con los nombres de las muestras
muestras_p <- samples_p$sample

# Muestras reales Control y JIA
control_real_p <- samples_p$sample[
  samples_p$condition_real == "Control"
]

jia_real_p <- samples_p$sample[
  samples_p$condition_real == "JIA"
]

# Crear todas las combinaciones posibles de 3 muestras
combinaciones_p <- combn(
  muestras_p,
  3,
  simplify = FALSE
)

# Función para crear una clave única por comparación
# Esto evita contar dos veces Grupo_A vs Grupo_B y Grupo_B vs Grupo_A
crear_clave_p <- function(grupoA_p) {
  
  grupoB_p <- setdiff(muestras_p, grupoA_p)
  
  claveA_p <- paste(sort(grupoA_p), collapse = "|")
  claveB_p <- paste(sort(grupoB_p), collapse = "|")
  
  clave_final_p <- paste(
    sort(c(claveA_p, claveB_p)),
    collapse = "__VS__"
  )
  
  return(clave_final_p)
}

# Crear claves para todas las combinaciones
claves_p <- sapply(
  combinaciones_p,
  crear_clave_p
)

# Eliminar comparaciones duplicadas
combinaciones_unicas_p <- combinaciones_p[
  !duplicated(claves_p)
]

cat("Número de combinaciones únicas 3 vs 3:", 
    length(combinaciones_unicas_p), "\n")


##################################################
# 9. EXCLUIR LA COMPARACION REAL CONTROL vs JIA
##################################################

es_comparacion_real_p <- function(grupoA_p) {
  
  setequal(grupoA_p, control_real_p) |
    setequal(grupoA_p, jia_real_p)
}

permutaciones_p <- combinaciones_unicas_p[
  !sapply(combinaciones_unicas_p, es_comparacion_real_p)
]

cat("Número de permutaciones artificiales:", 
    length(permutaciones_p), "\n")


###########################################
# 10. GUARDAR TABLA DE PERMUTACIONES
###########################################

tabla_permutaciones_p <- data.frame(
  Permutacion = paste0("P", seq_along(permutaciones_p)),
  Grupo_A = sapply(
    permutaciones_p,
    function(x) paste(x, collapse = ", ")
  ),
  Grupo_B = sapply(
    permutaciones_p,
    function(x) paste(setdiff(muestras_p, x), collapse = ", ")
  )
)

print(tabla_permutaciones_p)

write.xlsx(
  tabla_permutaciones_p,
  file.path(PERM_DIR, "tabla_permutaciones_generadas.xlsx"),
  rowNames = FALSE
)


########################################################
# 11. FUNCION PARA ANALISIS DESEQ2 EN CADA PERMUTACION
########################################################

analisis_deseq2_p <- function(grupoA_p, nombre_perm_p) {
  
  grupoB_p <- setdiff(muestras_p, grupoA_p)
  
  # Crear metadatos permutados
  metadata_perm_p <- data.frame(
    row.names = muestras_p,
    condition_perm = ifelse(
      muestras_p %in% grupoA_p,
      "Grupo_A",
      "Grupo_B"
    )
  )
  
  # Grupo_B será la referencia
  metadata_perm_p$condition_perm <- factor(
    metadata_perm_p$condition_perm,
    levels = c("Grupo_B", "Grupo_A")
  )
  
  # Asegurar que el orden de metadatos coincida con txi_p$counts
  metadata_perm_p <- metadata_perm_p[
    colnames(txi_p$counts),
    ,
    drop = FALSE
  ]
  
  # Crear objeto DESeq2
  dds_perm_p <- DESeqDataSetFromTximport(
    txi_p,
    colData = metadata_perm_p,
    design = ~ condition_perm
  )
  
  # Filtrado igual al análisis principal
  keep_perm_p <- rowSums(counts(dds_perm_p) >= 10) >= 3
  
  dds_perm_p <- dds_perm_p[
    keep_perm_p,
  ]
  
  # Ejecutar DESeq2
  dds_perm_p <- DESeq(
    dds_perm_p,
    quiet = TRUE
  )
  
  # Extraer resultados Grupo_A vs Grupo_B
  res_deseq_perm_p <- results(
    dds_perm_p,
    contrast = c("condition_perm", "Grupo_A", "Grupo_B")
  )
  
  res_deseq_perm_p <- as.data.frame(res_deseq_perm_p)
  res_deseq_perm_p$gene_id <- rownames(res_deseq_perm_p)
  
  # Añadir anotación génica
  res_deseq_perm_p <- merge(
    res_deseq_perm_p,
    gene_annotation_p,
    by = "gene_id",
    all.x = TRUE
  )
  
  # Ordenar por padj
  res_deseq_perm_p <- res_deseq_perm_p[
    order(res_deseq_perm_p$padj),
  ]
  
  # Genes significativos en la permutación
  sig_deseq_perm_p <- subset(
    res_deseq_perm_p,
    !is.na(padj) &
      padj < 0.05 &
      abs(log2FoldChange) > 1
  )
  
  # Agregar información de la permutación
  if (nrow(sig_deseq_perm_p) > 0) {
    
    sig_deseq_perm_p$Permutacion <- nombre_perm_p
    sig_deseq_perm_p$Grupo_A <- paste(grupoA_p, collapse = ", ")
    sig_deseq_perm_p$Grupo_B <- paste(grupoB_p, collapse = ", ")
    
  }
  
  # Contar genes up y down
  n_up_p <- sum(
    sig_deseq_perm_p$log2FoldChange > 1,
    na.rm = TRUE
  )
  
  n_down_p <- sum(
    sig_deseq_perm_p$log2FoldChange < -1,
    na.rm = TRUE
  )
  
  return(
    list(
      resultados = res_deseq_perm_p,
      significativos = sig_deseq_perm_p,
      n_significativos = nrow(sig_deseq_perm_p),
      n_up = n_up_p,
      n_down = n_down_p
    )
  )
}


#######################################################
# 12. FUNCION PARA ANALISIS EDGER EN CADA PERMUTACION
#######################################################

analisis_edger_p <- function(grupoA_p, nombre_perm_p) {
  
  grupoB_p <- setdiff(muestras_p, grupoA_p)
  
  # Crear grupos permutados
  group_perm_p <- ifelse(
    muestras_p %in% grupoA_p,
    "Grupo_A",
    "Grupo_B"
  )
  
  # Grupo_B será la referencia
  group_perm_p <- factor(
    group_perm_p,
    levels = c("Grupo_B", "Grupo_A")
  )
  
  # Matriz de conteos
  counts_edger_perm_p <- txi_p$counts
  
  # Crear objeto DGEList
  dge_perm_p <- DGEList(
    counts = counts_edger_perm_p,
    group = group_perm_p
  )
  
  # Filtrado igual al análisis principal
  keep_edger_perm_p <- rowSums(dge_perm_p$counts >= 10) >= 3
  
  dge_perm_p <- dge_perm_p[
    keep_edger_perm_p,
    ,
    keep.lib.sizes = FALSE
  ]
  
  # Normalización TMM
  dge_perm_p <- calcNormFactors(
    dge_perm_p,
    method = "TMM"
  )
  
  # Matriz de diseño
  design_perm_p <- model.matrix(
    ~ group_perm_p
  )
  
  # Estimación de dispersión
  dge_perm_p <- estimateDisp(
    dge_perm_p,
    design_perm_p
  )
  
  # Ajuste del modelo GLM
  fit_perm_p <- glmFit(
    dge_perm_p,
    design_perm_p
  )
  
  # Prueba de razón de verosimilitud
  lrt_perm_p <- glmLRT(
    fit_perm_p,
    coef = 2
  )
  
  # Extraer resultados
  res_edger_perm_p <- as.data.frame(
    topTags(
      lrt_perm_p,
      n = Inf
    )
  )
  
  res_edger_perm_p$gene_id <- rownames(res_edger_perm_p)
  
  # Añadir anotación génica
  res_edger_perm_p <- merge(
    res_edger_perm_p,
    gene_annotation_p,
    by = "gene_id",
    all.x = TRUE
  )
  
  # Ordenar por FDR
  res_edger_perm_p <- res_edger_perm_p[
    order(res_edger_perm_p$FDR),
  ]
  
  # Genes significativos en la permutación
  sig_edger_perm_p <- subset(
    res_edger_perm_p,
    !is.na(FDR) &
      FDR < 0.05 &
      abs(logFC) > 1
  )
  
  # Agregar información de la permutación
  if (nrow(sig_edger_perm_p) > 0) {
    
    sig_edger_perm_p$Permutacion <- nombre_perm_p
    sig_edger_perm_p$Grupo_A <- paste(grupoA_p, collapse = ", ")
    sig_edger_perm_p$Grupo_B <- paste(grupoB_p, collapse = ", ")
    
  }
  
  # Contar genes up y down
  n_up_p <- sum(
    sig_edger_perm_p$logFC > 1,
    na.rm = TRUE
  )
  
  n_down_p <- sum(
    sig_edger_perm_p$logFC < -1,
    na.rm = TRUE
  )
  
  return(
    list(
      resultados = res_edger_perm_p,
      significativos = sig_edger_perm_p,
      n_significativos = nrow(sig_edger_perm_p),
      n_up = n_up_p,
      n_down = n_down_p
    )
  )
}


###########################################
# 13. EJECUTAR TODAS LAS PERMUTACIONES
##########################################

resumen_permutaciones_p <- data.frame()

sig_deseq_todas_p <- data.frame()
sig_edger_todas_p <- data.frame()

for (i in seq_along(permutaciones_p)) {
  
  nombre_perm_p <- paste0("P", i)
  
  grupoA_p <- permutaciones_p[[i]]
  grupoB_p <- setdiff(muestras_p, grupoA_p)
  
  cat("\n====================================\n")
  cat("Ejecutando permutación:", nombre_perm_p, "\n")
  cat("Grupo A:", paste(grupoA_p, collapse = ", "), "\n")
  cat("Grupo B:", paste(grupoB_p, collapse = ", "), "\n")
  cat("====================================\n")
  
  # Ejecutar DESeq2
  deseq_perm_p <- analisis_deseq2_p(
    grupoA_p,
    nombre_perm_p
  )
  
  # Ejecutar edgeR
  edger_perm_p <- analisis_edger_p(
    grupoA_p,
    nombre_perm_p
  )
  
  # Crear fila resumen
  resumen_temp_p <- data.frame(
    Permutacion = nombre_perm_p,
    Grupo_A = paste(grupoA_p, collapse = ", "),
    Grupo_B = paste(grupoB_p, collapse = ", "),
    
    FP_empiricos_DESeq2 = deseq_perm_p$n_significativos,
    Up_DESeq2 = deseq_perm_p$n_up,
    Down_DESeq2 = deseq_perm_p$n_down,
    
    FP_empiricos_edgeR = edger_perm_p$n_significativos,
    Up_edgeR = edger_perm_p$n_up,
    Down_edgeR = edger_perm_p$n_down
  )
  
  resumen_permutaciones_p <- rbind(
    resumen_permutaciones_p,
    resumen_temp_p
  )
  
  # Guardar genes significativos DESeq2
  if (nrow(deseq_perm_p$significativos) > 0) {
    
    sig_deseq_todas_p <- rbind(
      sig_deseq_todas_p,
      deseq_perm_p$significativos
    )
    
  }
  
  # Guardar genes significativos edgeR
  if (nrow(edger_perm_p$significativos) > 0) {
    
    sig_edger_todas_p <- rbind(
      sig_edger_todas_p,
      edger_perm_p$significativos
    )
    
  }
}


############################
# 14. RESUMEN FINAL
############################

resumen_final_p <- data.frame(
  Herramienta = c("DESeq2", "edgeR"),
  
  Promedio_FP_empiricos = c(
    mean(resumen_permutaciones_p$FP_empiricos_DESeq2),
    mean(resumen_permutaciones_p$FP_empiricos_edgeR)
  ),
  
  Mediana_FP_empiricos = c(
    median(resumen_permutaciones_p$FP_empiricos_DESeq2),
    median(resumen_permutaciones_p$FP_empiricos_edgeR)
  ),
  
  Minimo_FP_empiricos = c(
    min(resumen_permutaciones_p$FP_empiricos_DESeq2),
    min(resumen_permutaciones_p$FP_empiricos_edgeR)
  ),
  
  Maximo_FP_empiricos = c(
    max(resumen_permutaciones_p$FP_empiricos_DESeq2),
    max(resumen_permutaciones_p$FP_empiricos_edgeR)
  ),
  
  Total_FP_empiricos = c(
    sum(resumen_permutaciones_p$FP_empiricos_DESeq2),
    sum(resumen_permutaciones_p$FP_empiricos_edgeR)
  )
)

print(resumen_permutaciones_p)
print(resumen_final_p)


############################
# 15. EXPORTAR RESULTADOS
############################

write.xlsx(
  resumen_permutaciones_p,
  file.path(PERM_DIR, "resumen_permutaciones_falsos_positivos.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  resumen_final_p,
  file.path(PERM_DIR, "resumen_final_falsos_positivos_empiricos.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  sig_deseq_todas_p,
  file.path(PERM_DIR, "genes_significativos_permutaciones_DESeq2.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  sig_edger_todas_p,
  file.path(PERM_DIR, "genes_significativos_permutaciones_edgeR.xlsx"),
  rowNames = FALSE
)


###########################################
# 16. GRAFICO COMPARATIVO BARRA Y PUNTO
###########################################

grafico_perm_p <- ggplot(
  resumen_permutaciones_p,
  aes(x = Permutacion)
) +
  geom_col(
    aes(y = FP_empiricos_DESeq2),
    fill="#F4A261",
    alpha = 0.7
  ) +
  geom_point(
    aes(y = FP_empiricos_edgeR),
    color = "#E76F51",
    size = 2.5
  ) +
  theme_bw() +
  labs(
    title = "Falsos positivos empíricos por permutación",
    subtitle = "Barras: DESeq2 | Puntos: edgeR",
    x = "Permutación",
    y = "Número de genes significativos"
  )

print(grafico_perm_p)

ggsave(
  filename = file.path(PERM_DIR, "grafico_falsos_positivos_permutaciones.png"),
  plot = grafico_perm_p,
  width = 8,
  height = 6,
  dpi = 300
)

###########################################
# GRAFICO COMPARATIVO DE BARRAS
###########################################

grafico_barras_p <- ggplot(
  resumen_largo_p,
  aes(
    x = Permutacion,
    y = Falsos_positivos,
    fill = Herramienta
  )
) +
  geom_col(
    position = "dodge",
    alpha = 0.8
  )  +
  scale_fill_manual(
    values = c(
      "DESeq2" = "#F4A261",
      "edgeR" = "#E76F51"
    )
  ) +
  theme_bw() +
  labs(
    title = "Comparación de falsos positivos empíricos por permutación",
    x = "Permutación",
    y = "Número de genes significativos",
    fill = "Herramienta"
  )

print(grafico_barras_p)

ggsave(
  filename = file.path(PERM_DIR, "grafico_falsos_positivos_permutaciones2.png"),
  plot = grafico_barras_p,
  width = 8,
  height = 6,
  dpi = 300
)

#########################
#GUARDAR ENTORNO DE R
save.image()
#########################
