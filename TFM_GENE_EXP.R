# --------------------------------------------------

# RNA-seq Analysis in R

# Date: 14-06-2026

# Author: Melannie Cedeño Valencia

# --------------------------------------------------



# Description:

# This script performs RNA-seq data analysis in R, 

# including metadata processing and raw count matrix import.


# --------------------------------------------------



# Set work directory

RESULTS_DIR <- "C:/Users/mcede/Documents/TFM_2026/R_GENE_EXPRESSION"
SALMON_DIR <- "C:/Users/mcede/Documents/TFM_2026/salmon_results"

REFERENCE_DIR <- "C:/Users/mcede/Documents/TFM_2026/reference"

# Definir directorio de trabajo para guardar salidas
setwd(RESULTS_DIR)



##Llamado de librerias

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
library(EnhancedVolcano)


#########################

##CARGA DE ANOTACION

#########################
# Archivo GTF
gtf_file <- file.path(
  REFERENCE_DIR,
  "gencode.v49.annotation.gtf.gz"
)

cat("Importando GTF...\n")

gtf <- import(gtf_file)


tx2gene <- unique(
  data.frame(
    transcript_id = mcols(gtf)$transcript_id,
    gene_id = mcols(gtf)$gene_id,
    gene_name = mcols(gtf)$gene_name,
    gene_type = mcols(gtf)$gene_type
  )
)

tx2gene <- na.omit(tx2gene)

cat(
  "Numero de transcritos en tx2gene:",
  nrow(tx2gene),
  "\n"
)


head(tx2gene)

#########################

#Tabla de anotaciones

#########################

gene_annotation <- unique(
  tx2gene[, c("gene_id", "gene_name","gene_type")]
)


head(gene_annotation)


write.csv(
  gene_annotation,
  file.path(
    RESULTS_DIR,
    "gene_annotation_gencode_v49.csv"
  ),
  row.names = FALSE
)


#########################

#METADATOS

#########################

samples <- data.frame(
  sample = c(
    "SRR6006895",
    "SRR6006896",
    "SRR6006898",
    "SRR6006900",
    "SRR6006904",
    "SRR6006906"
  ),
  condition = c(
    "Control",
    "Control",
    "Control",
    "JIA",
    "JIA",
    "JIA"
  )
)

metadata <- data.frame(
  row.names = samples$sample,
  condition = factor(
    samples$condition,
    levels = c("Control", "JIA")
  )
)

write.csv(
  metadata,
  file.path(
    RESULTS_DIR,
    "metadata_samples.csv"
  )
)

print(metadata)

#########################

##ARCHIVOS SALMON

#########################

files <- file.path(
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

names(files) <- samples$sample

print(file.exists(files))

list.files(SALMON_DIR)

#########################

#VERIFICACION DE QUANT.SF

#########################

quant_test <- read.delim(files[1])

cat(
  "Numero de transcritos cuantificados:",
  nrow(quant_test),
  "\n"
)

head(quant_test)

#########################

#IMPORTACION SALMON

#########################

cat("Importando resultados Salmon...\n")

txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene[, c("transcript_id","gene_id","gene_type")]
)

cat(
  "Dimensiones matriz inicial de tximport:",
  dim(txi$counts),
  "\n"
)

head(txi$counts)

#########################

#CREACION OBJETO DESEQ

#########################

dds <- DESeqDataSetFromTximport(
  txi,
  colData = metadata,
  design = ~ condition
)

#########################
# FILTRADO DE BAJA EXPRESION
#########################

cat(
  "Genes antes del filtrado:",
  nrow(dds),
  "\n"
)

keep <- rowSums(counts(dds) >= 10) >= 3

cat(
  "Genes retenidos despues del filtrado:",
  sum(keep),
  "\n"
)

dds <- dds[keep, ]


cat(
  "Dimensiones finales del objeto dds:",
  dim(dds),
  "\n"
)

#########################
# ANALISIS DESEQ2
#########################

dds <- DESeq(dds)

names(mcols(dds))
############################
# GRAFICO DE DISPERSION
############################


png(
  filename = "DESeq2_dispersion_plot.png",
  width = 900,
  height = 700
)

plotDispEsts(dds)

dev.off()

#########################
# CONTEOS NORMALIZADOS
#########################

normalized_counts <- counts(
  dds,
  normalized = TRUE
)

write.csv(
  normalized_counts,
  file.path(
    RESULTS_DIR,
    "DESeq2_normalized_counts.csv"
  )
)

cat(
  "\nDimensiones matriz normalizada:",
  dim(normalized_counts),
  "\n"
)

#########################
# RESULTADOS DESEQ2
#########################

res_deseq <- results(
  dds,
  contrast = c(
    "condition",
    "JIA",
    "Control"
  )
)

# Convertir a data.frame
res_deseq <- as.data.frame(res_deseq)

# Agregar gene_id y gene_type
res_deseq$gene_id <- rownames(res_deseq)


# Añadir nombres de genes
res_deseq <- merge(
  res_deseq,
  gene_annotation,
  by = "gene_id",
  all.x = TRUE
)
dim(res_deseq)

str(res_deseq)

cat(
  "\nGenes dentro de resdeseq:",
  nrow(res_deseq),
  "\n"
)

# Eliminar genes sin estadísticos calculados
res_deseq <- res_deseq[
  !is.na(res_deseq$padj),
]

cat(
  "\nGenes después de eliminar genes sin estadísticos:",
  nrow(res_deseq),
  "\n"
)

# Ordenar por significancia
res_deseq <- res_deseq[
  order(res_deseq$padj),
]


write.xlsx(
  res_deseq,
  "DESeq2_results.xlsx",
  rowNames = FALSE
)

#########################
# GENES SIGNIFICATIVOS
#########################

sig_deseq <- subset(
  res_deseq,
  padj < 0.05 &
    abs(log2FoldChange) > 1
)
sig_deseq <- sig_deseq[
  order(sig_deseq$padj),
]

write.xlsx(
  sig_deseq,
  "DESeq2_significant_genes.xlsx",
  rowNames = FALSE
)

###############################################
# RESUMEN DESEQ2 - TOP 20 MAS SIGNIFICATIVOS
###############################################

cat(
  "\nGenes analizados:",
  nrow(res_deseq),
  "\n"
)

cat(
  "Genes significativos (padj < 0.05 y |log2FC| > 1):",
  nrow(sig_deseq),
  "\n"
)

cat(
  "\nTop 20 genes más significativos:\n"
)

sigdeseq_top20<-head(
    sig_deseq[
      ,
      c(
        "gene_name",
        "gene_id",
        "gene_type",
        "log2FoldChange",
        "padj"
      )
    ],
    20
  )

print(sigdeseq_top20)

#########################
# EXPORTAR TABLA TOP 20
#########################

write.xlsx(
  sigdeseq_top20,
  file.path(
    RESULTS_DIR,
    "SigDeseq2_top20_genes.xlsx"
  ),
  rowNames = FALSE
)


##########################################
# RESUMEN DESEQ2 - TOP 50
##########################################

cat(
  "\nTop 50 genes más significativos:\n"
)

sigdeseq_top50 <- head(
    sig_deseq[
      ,
      c(
        "gene_name",
        "gene_id",
        "gene_type",
        "log2FoldChange",
        "padj"
      )
    ],
    50
  )

print(sigdeseq_top50)

#########################
# EXPORTAR TABLA TOP 50
#########################

write.xlsx(
  sigdeseq_top50,
  file.path(
    RESULTS_DIR,
    "SigDESeq2_top50_genes.xlsx"
  ),
  rowNames = FALSE
)


#########################
# RESUMEN Sig_DESEQ2
#########################

cat(
  "\nGenes analizados:",
  nrow(res_deseq),
  "\n"
)

cat(
  "Genes significativos (padj < 0.05 y |log2FC| > 1):",
  nrow(sig_deseq),
  "\n"
)


##############################
# GENES UPREGULATED EN JIA
##############################

sig_upregulated <- subset(
  sig_deseq,
  padj < 0.05 &
    log2FoldChange > 1
)

head(
  sig_upregulated[
    ,
    c(
      "gene_name",
      "gene_id",
      "gene_type",
      "log2FoldChange",
      "padj"
    )
  ],
  20
)

write.xlsx(
  sig_upregulated,
  file.path(
    RESULTS_DIR,
    "Sig_DESeq2_upregulated_JIA.xlsx"
  ),
  rowNames = FALSE
)

cat(
  "Genes upregulated en JIA:",
  nrow(sig_upregulated),
  "\n"
)

##############################
# GENES DOWNREGULATED EN JIA
##############################

sig_downregulated <- subset(
  sig_deseq,
  padj < 0.05 &
    log2FoldChange < -1
)

head(
  sig_downregulated[
    ,
    c(
      "gene_name",
      "gene_id",
      "gene_type",
      "log2FoldChange",
      "padj"
    )
  ],
  20
)

write.xlsx(
  sig_downregulated,
  file.path(
    RESULTS_DIR,
    "Sig_DESeq2_downregulated_JIA.xlsx"
  ),
  rowNames = FALSE
)

cat(
  "Genes downregulated en JIA:",
  nrow(sig_downregulated),
  "\n"
)



###############################
# Genes presentes en JIA
# y no detectados en Control
# usando umbral de 10 conteos
###############################

# Matriz de conteos crudos
counts_raw <- counts(dds)

# Definir muestras
control_samples <- c("SRR6006895", "SRR6006896", "SRR6006898")
jia_samples <- c("SRR6006900", "SRR6006904", "SRR6006906")

# Verificar que las muestras existan en la matriz de conteos
all(control_samples %in% colnames(counts_raw))
all(jia_samples %in% colnames(counts_raw))

# Genes con >=10 conteos en al menos 2 muestras JIA
# y <10 conteos en todas las muestras Control
genes_jia_only_10 <- rownames(counts_raw)[
  rowSums(counts_raw[, jia_samples, drop = FALSE] >= 10) >= 2 &
    rowSums(counts_raw[, control_samples, drop = FALSE] >= 10) == 0
]

# Filtrar con estadística DESeq2
genes_jia_only_10_significativos <- res_deseq[
  res_deseq$gene_id %in% genes_jia_only_10 &
    res_deseq$padj < 0.05 &
    res_deseq$log2FoldChange > 0 &
    res_deseq$stat > 0,
]

# Ordenar por padj
genes_jia_only_10_significativos <- genes_jia_only_10_significativos[
  order(genes_jia_only_10_significativos$padj),
]

# Extraer conteos de esos genes
counts_jia_only_10 <- counts_raw[
  genes_jia_only_10_significativos$gene_id,
  ,
  drop = FALSE
]

# Convertir conteos a data.frame
counts_jia_only_10_df <- as.data.frame(counts_jia_only_10)

# Agregar gene_id a la tabla de conteos
counts_jia_only_10_df$gene_id <- rownames(counts_jia_only_10_df)

# Unir resultados DESeq2 + conteos
genes_jia_only_10_final <- merge(
  genes_jia_only_10_significativos,
  counts_jia_only_10_df,
  by = "gene_id",
  all.x = TRUE
)

# Ordenar nuevamente por padj después del merge
genes_jia_only_10_final <- genes_jia_only_10_final[
  order(genes_jia_only_10_final$padj),
]

# Ver dimensiones y primeros resultados
dim(genes_jia_only_10_final)
head(genes_jia_only_10_final)

# Exportar resultados
write.xlsx(
  genes_jia_only_10_final,
  file = "genes_presentes_JIA_no_control_umbral10.xlsx",
  rowNames = FALSE
)

#########################

# PCA

#########################

vsd <- vst(dds)

pca_data <- plotPCA(
  vsd,
  intgroup = "condition",
  returnData = TRUE
)

percentVar <- round(
  100 * attr(pca_data, "percentVar")
)

pca_plot <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    color = condition
  )
) +
  geom_point(size = 4) +
  xlab(
    paste0(
      "PC1: ",
      percentVar[1],
      "%"
    )
  ) +
  ylab(
    paste0(
      "PC2: ",
      percentVar[2],
      "%"
    )
  ) +
  theme_bw()

ggsave(
  filename = file.path(
    RESULTS_DIR,
    "PCA_samples.png"
  ),
  plot = pca_plot,
  width = 8,
  height = 6
)


#########################

# HEATMAP DE MUESTRAS

#########################

sampleDistMatrix <- dist(
  t(
    assay(vsd)
  )
)

pheatmap(
  as.matrix(sampleDistMatrix),
  filename = file.path(
    RESULTS_DIR,
    "Sample_distance_heatmap.png"
  )
)


##################################

# VOLCANO DESEQ2 - con res_deseq

#################################


# PREPARAR DATOS PARA VOLCANO


volcano_df <- res_deseq

# Eliminar genes sin valores necesarios
volcano_df <- volcano_df[
  !is.na(volcano_df$log2FoldChange) &
    !is.na(volcano_df$padj),
]

# Evitar problemas si algún padj es 0
volcano_df$padj[volcano_df$padj == 0] <- min(
  volcano_df$padj[volcano_df$padj > 0],
  na.rm = TRUE
)

# Crear etiquetas: usar gene_name si existe, si no usar gene_id
volcano_df$label_gene <- ifelse(
  is.na(volcano_df$gene_name) | volcano_df$gene_name == "",
  volcano_df$gene_id,
  volcano_df$gene_name
)

############################
# GENES A ETIQUETAR
############################

genes_label <- sig_deseq$gene_name

genes_label <- genes_label[
  !is.na(genes_label) &
    genes_label != ""
]

genes_label <- head(genes_label, 20)

############################
# COLORES PERSONALIZADOS
############################

keyvals <- ifelse(
  volcano_df$log2FoldChange >= 1 & volcano_df$padj < 0.05,
  "red",
  ifelse(
    volcano_df$log2FoldChange <= -1 & volcano_df$padj < 0.05,
    "forestgreen",
    "grey70"
  )
)

names(keyvals)[keyvals == "red"] <- "Upregulated en JIA"
names(keyvals)[keyvals == "forestgreen"] <- "Downregulated en JIA"
names(keyvals)[keyvals == "grey70"] <- "No significativo"

############################
# CREAR VOLCANO PLOT
############################

volcano_enhanced2 <- EnhancedVolcano(
  volcano_df,
  lab = volcano_df$label_gene,
  x = "log2FoldChange",
  y = "padj",
  title = "Volcano plot DESeq2",
  subtitle = "JIA vs Control",
  caption = "Criterios: padj < 0.05 y |log2FoldChange| > 1",
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 2.0,
  labSize = 3.5,
  selectLab = genes_label,
  xlab = "log2FoldChange",
  ylab = "-log10(padj)",
  legendPosition = "right",
  legendLabSize = 10,
  legendIconSize = 4,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colCustom = keyvals
)

print(volcano_enhanced2)

############################
# EXPORTAR PNG
############################

ggsave(
  filename = "Volcano_Enhanced_DESeq2_JIA_vs_Control_2_2.png",
  plot = volcano_enhanced2,
  width = 8,
  height = 7,
  dpi = 300
)


###############################################################
#RESUMEN Deseq2
###############################################################


cat("\n RESUMEN DESeq2")


cat("\nGenes totales analizados:", nrow(res_deseq), "\n")

cat(
  "Genes con padj no NA:",
  sum(!is.na(res_deseq$padj)),
  "\n"
)

cat(
  "Genes significativos DESeq2 (padj < 0.05 y |log2FC| > 1):",
  nrow(sig_deseq),
  "\n"
)

cat(
  "Genes upregulated DESeq2 en JIA:",
  nrow(sig_upregulated),
  "\n"
)

cat(
  "Genes downregulated DESeq2 en JIA:",
  nrow(sig_downregulated),
  "\n"
)


###############################################################
# ANALISIS DE EXPRESION DIFERENCIAL CON edgeR
# JIA vs CONTROL
# A partir de resultados Salmon importados con tximport
###############################################################

###############################
# DIRECTORIO RESULTADOS edgeR
###############################

OUTPUT_DIR_2 <- "C:/Users/mcede/Documents/TFM_2026/R_GENE_EXPRESSION/edgeR"


if (!dir.exists(OUTPUT_DIR_2)) {
  dir.create(OUTPUT_DIR_2, recursive = TRUE)
}


###############################################################
# MATRIZ DE CONTEOS PARA edgeR
###############################################################

counts_edger <- txi$counts

cat(
  "Dimensiones matriz inicial para edgeR:",
  dim(counts_edger),
  "\n"
)

###############################################################
# DEFINIR GRUPOS
###############################################################

group <- factor(
  samples$condition,
  levels = c("Control", "JIA")
)

cat(
  "\nGrupos experimentales:\n"
)

print(group)

###############################################################
# CREAR OBJETO DGEList
###############################################################

dge <- DGEList(
  counts = counts_edger,
  group = group
)

cat(
  "\nDimensiones iniciales del objeto edgeR:\n"
)

print(dim(dge))

###############################################################
# FILTRADO DE BAJA EXPRESION
###############################################################

keep_edger <- rowSums(dge$counts >= 10) >= 3

cat(
  "\nGenes antes del filtrado edgeR:",
  nrow(dge),
  "\n"
)

cat(
  "Genes retenidos despues del filtrado edgeR:",
  sum(keep_edger),
  "\n"
)

cat(
  "Genes retenidos despues del filtrado Deseq2:",
  sum(keep),
  "\n"
)


dge <- dge[
  keep_edger,
  ,
  keep.lib.sizes = FALSE
]

cat(
  "Dimensiones despues del filtrado edgeR:",
  dim(dge),
  "\n"
)

###############################################################
# NORMALIZACION TMM
###############################################################

dge <- calcNormFactors(
  dge,
  method = "TMM"
)

cat(
  "\nFactores de normalizacion TMM:\n"
)

print(dge$samples)

###############################################################
# GUARDAR MATRIZ CPM NORMALIZADA
###############################################################

cpm_edger <- cpm(
  dge,
  normalized.lib.sizes = TRUE
)

write.xlsx(
  as.data.frame(cpm_edger),
  file.path(
    OUTPUT_DIR_2,
    "edgeR_CPM_normalized_counts.xlsx"
  ),
  rowNames = TRUE
)

###############################################################
# MATRIZ DE DISEÑO
###############################################################

design <- model.matrix(
  ~ group
)

cat(
  "\nMatriz de diseño edgeR:\n"
)

print(design)

###############################################################
# ESTIMACION DE DISPERSION
###############################################################

dge <- estimateDisp(
  dge,
  design
)




###############################################################
# GRAFICO BCV
###############################################################

png(
  filename = file.path(
    OUTPUT_DIR_2,
    "edgeR_BCV_plot.png"
  ),
  width = 1200,
  height = 900,
  res = 150
)

plotBCV(dge)

dev.off()
###############################################################
# AJUSTE DEL MODELO GLM
# Metodo menos conservador: Likelihood Ratio Test
###############################################################

fit <- glmFit(
  dge,
  design
)
###############################################################
# TEST DE EXPRESION DIFERENCIAL
# coef = 2 corresponde a groupJIA
# Comparacion: JIA vs Control
###############################################################

lrt <- glmLRT(
  fit,
  coef = 2
)

###############################################################
# RESULTADOS edgeR
###############################################################

res_edger <- as.data.frame(
  topTags(
    lrt,
    n = Inf
  )
)

# Agregar gene_id
res_edger$gene_id <- rownames(res_edger)

# Añadir gene_name y gene_type desde la anotación
res_edger <- merge(
  res_edger,
  gene_annotation,
  by = "gene_id",
  all.x = TRUE
)

# Ordenar por FDR
res_edger <- res_edger[
  order(res_edger$FDR),
]

cat(
  "\nGenes dentro de res_edger:",
  nrow(res_edger),
  "\n"
)

str(res_edger)
###############################################################
# GUARDAR RESULTADOS COMPLETOS edgeR
###############################################################

write.xlsx(
  res_edger,
  file.path(
    OUTPUT_DIR_2,
    "edgeR_results.xlsx"
  ),
  rowNames = FALSE
)

###############################################################
# GENES SIGNIFICATIVOS edgeR
###############################################################

sig_edger <- subset(
  res_edger,
  FDR < 0.05 &
    abs(logFC) > 1
)

sig_edger <- sig_edger[
  order(sig_edger$FDR),
]

cat(
  "\nGenes significativos edgeR (FDR < 0.05 y |logFC| > 1):",
  nrow(sig_edger),
  "\n"
)

write.xlsx(
  sig_edger,
  file.path(
    OUTPUT_DIR_2,
    "edgeR_significant_genes.xlsx"
  ),
  rowNames = FALSE
)


###############################################################
# GENES UPREGULATED EN JIA - edgeR
###############################################################

sig_edger_upregulated_JIA <- subset(
  res_edger,
  FDR < 0.05 &
    logFC > 1
)

sig_edger_upregulated_JIA <- sig_edger_upregulated_JIA[
  order(sig_edger_upregulated_JIA$FDR),
]

cat(
  "Genes upregulated en JIA por edgeR:",
  nrow(sig_edger_upregulated_JIA),
  "\n"
)

write.xlsx(
  sig_edger_upregulated_JIA,
  file.path(
    OUTPUT_DIR_2,
    "edgeR_significant_upregulated_JIA.xlsx"
  ),
  rowNames = FALSE
)

###############################################################
# GENES DOWNREGULATED EN JIA - edgeR
###############################################################

sig_edger_downregulated_JIA <- subset(
  res_edger,
  FDR < 0.05 &
    logFC < -1
)

sig_edger_downregulated_JIA <- sig_edger_downregulated_JIA[
  order(sig_edger_downregulated_JIA$FDR),
]

cat(
  "Genes downregulated en JIA por edgeR:",
  nrow(sig_edger_downregulated_JIA),
  "\n"
)

write.xlsx(
  sig_edger_downregulated_JIA,
  file.path(
    OUTPUT_DIR_2,
    "edgeR_significant_downregulated_JIA.xlsx"
  ),
  rowNames = FALSE
)

###############################################################
# TOP 20 GENES edgeR
###############################################################

edger_top20 <- head(
  res_edger[
    ,
    c(
      "gene_name",
      "gene_id",
      "gene_type",
      "logFC",
      "logCPM",
      "LR",
      "PValue",
      "FDR"
    )
  ],
  20
)

cat(
  "\nTop 20 genes más significativos por edgeR LRT:\n"
)

print(edger_top20)

write.xlsx(
  edger_top20,
  file.path(
    OUTPUT_DIR_2,
    "edgeR_top20_genes.xlsx"
  ),
  rowNames = FALSE
)

###############################################################
# VOLCANO PLOT edgeR
###############################################################

volcano_edger <- res_edger

############################
# GENES A ETIQUETAR edgeR
############################

genes_label_edger <- sig_edger$gene_name

genes_label_edger <- genes_label_edger[
  !is.na(genes_label_edger) &
    genes_label_edger != "" &
    genes_label_edger != sig_edger$gene_id
]

genes_label_edger <- head(  genes_label_edger,  20)

print(genes_label_edger)

###############################################################
# CREAR ETIQUETAS
###############################################################

volcano_edger$label_gene <- ifelse(
  is.na(volcano_edger$gene_name) |
    volcano_edger$gene_name == "" |
    volcano_edger$gene_name == volcano_edger$gene_id,
  volcano_edger$gene_id,
  volcano_edger$gene_name
)
###############################################################
# COLORES PERSONALIZADOS
###############################################################

keyvals <- ifelse(
  volcano_edger$logFC >= 1 & volcano_edger$FDR < 0.05,
  "orange",
  ifelse(
    volcano_edger$logFC <= -1 & volcano_edger$FDR < 0.05,
    "darkgreen",
    "grey70"
  )
)

names(keyvals)[keyvals == "orange"] <- "Upregulated en JIA"
names(keyvals)[keyvals == "darkgreen"] <- "Downregulated en JIA"
names(keyvals)[keyvals == "grey70"] <- "No significativo"

###############################################################
# CREAR VOLCANO PLOT
###############################################################


volcano_enhanced_edger <- EnhancedVolcano(
  volcano_edger,
  lab = volcano_edger$label_gene,
  selectLab = genes_label_edger,
  x = "logFC",
  y = "FDR",
  title = "Volcano plot edgeR",
  subtitle = "JIA vs Control",
  caption = "Criterios: FDR < 0.05 y |logFC| > 1",
  pCutoff = 0.05,
  FCcutoff = 1,
  colCustom = keyvals,
  pointSize = 2.0,
  labSize = 3.5,
  xlab = "logFC",
  ylab = "-log10 FDR",
  legendPosition = "right",
  drawConnectors = TRUE
)

print(volcano_enhanced_edger)

###############################################################
# GUARDAR VOLCANO PLOT
###############################################################

ggsave(
  filename = file.path(
    OUTPUT_DIR_2,
    "edgeR_VolcanoPlot_EnhancedVolcano.png"
  ),
  plot = volcano_enhanced_edger,
  width = 9,
  height = 7,
  dpi = 300
)

###############################################################
# RESUMEN FINAL edgeR
###############################################################

cat(
  "\nResumen edgeR:\n"
)

cat(
  "Genes analizados por edgeR:",
  nrow(res_edger),
  "\n"
)

cat(
  "Genes significativos edgeR (FDR < 0.05 y |logFC| > 1):",
  nrow(sig_edger),
  "\n"
)

cat(
  "Genes upregulated en JIA por edgeR:",
  nrow(sig_edger_upregulated_JIA),
  "\n"
)

cat(
  "Genes downregulated en JIA por edgeR:",
  nrow(sig_edger_downregulated_JIA),
  "\n"
)


###############################################################
# COMPARACION DE LOG2FOLDCHANGE ENTRE LOS DOS METODOS
#EVALUAR LA CORRELACIÓN DE LOS DATOS
###############################################################


comparison <- merge(
  res_deseq[, c("gene_id", "log2FoldChange", "padj")],
  res_edger[, c("gene_id", "logFC", "FDR")],
  by = "gene_id"
)

cor(
  comparison$log2FoldChange,
  comparison$logFC,
  use = "complete.obs"
)

#########################
# GUARDAR ENTORNO DE R
#########################

save.image(
  file.path(
    RESULTS_DIR,
    "all_data.RData"
  )
)

cat(
  "\nAnalisis completado. Todos los archivos generados por R fueron guardados en:",
  RESULTS_DIR,
  "\n"
)



