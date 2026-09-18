#!/bin/bash

# Carpeta donde están los FASTQ
FASTQ_DIR="."

# Carpeta de salida para trimming
TRIM_DIR="../../TRIMMED"

# Crear carpeta de salida
mkdir -p "$TRIM_DIR"

# Número de hilos
THREADS=6

echo "Iniciando trimming de lecturas paired-end..."

# Buscar todos los archivos R1

for R1 in ${FASTQ_DIR}/*_1.fastq.gz
do
    [ -e "$R1" ] || continue

    BASE=$(basename "$R1" _1.fastq.gz)
    R2="${FASTQ_DIR}/${BASE}_2.fastq.gz"

    echo "Procesando muestra: $BASE"
    echo "R1: $R1"
    echo "R2: $R2"

    if [[ ! -f "$R2" ]]; then
        echo "ERROR: falta R2 para $BASE"
        continue
    fi

    trim_galore \
        --paired \
        --quality 20 \
        --length 30 \
        --trim-n \
        --cores $THREADS \
        --output_dir "$TRIM_DIR" \
        "$R1" "$R2"

done

echo "Trimming completado."


###
