#!/bin/bash

### Control de calidad
# Carpeta donde se encuentran los archivos
FASTQ_DIR="SRA_downloads"

# Carpeta de salida
OUTDIR="../../QC"

# Número de hilos
THREADS=6

echo "Iniciando análisis de calidad..."

# Analizar todos los fastq.gz

for FILE in ${FASTQ_DIR}/*.fastq.gz
do
    echo "Analizando: $FILE"

    fastqc "$FILE" \
        --threads $THREADS \
        --outdir "$OUTDIR"

done

echo "Análisis FastQC completado."

# Generar reporte MultiQC

echo "Generando reporte MultiQC..."

multiqc "$OUTDIR" -o "$OUTDIR"

echo "Proceso finalizado."
