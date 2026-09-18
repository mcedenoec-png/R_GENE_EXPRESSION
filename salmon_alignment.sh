#!/bin/bash

# ==========================================
# RNA-seq Cuantificación con SALMON
# ==========================================

# ----------- CONFIGURACIÓN ----------------

# Carpeta con reads trimmed
READS_DIR="$HOME/GSE103501/TRIMMED"

# Transcriptoma humano - referencia
TRANSCRIPTOME="$HOME/GSE103501/REFERENCE/gencode.v49.transcripts.fa.gz"

# Carpeta donde se guardará el índice
INDEX_DIR="$HOME/GSE103501/ALIGNMENT/gencode_index"

# Carpeta de resultados
OUTPUT_DIR="$HOME/GSE103501/ALIGNMENT/salmon_results"

# Número de threads
THREADS=4

# Tamaño del k-mer
KMER=31


# ------------------------------------------

echo "======================================="
echo " INICIANDO PIPELINE DE SALMON "
echo "======================================="

# Crear carpeta de resultados
mkdir -p “$INDEX_DIR”
mkdir -p "$OUTPUT_DIR"

# ==========================================
# PASO 1: Se crea el índice
# ==========================================

if [ ! -f "$INDEX_DIR/versionInfo.json" ]; then

    echo ""
    echo "Creando índice de Salmon..."
    echo "Transcriptoma: $TRANSCRIPTOME"
    echo "Índice: $INDEX_DIR"

    salmon index \
        -t "$TRANSCRIPTOME" \
        -i "$INDEX_DIR" \
        -k "$KMER" \
        -p "$THREADS"

    echo "Índice creado correctamente."

else

    echo ""
    echo "El índice ya existe. Saltando creación."

fi

# ==========================================
# PASO 2: Se cueantifica las muestras
# ==========================================

echo ""
echo "Iniciando cuantificación de muestras..."

for R1 in ${READS_DIR}/*_1_val_1.fq.gz
do

    # Obtener nombre de muestra
    SAMPLE=$(basename $R1 _1_val_1.fq.gz)

    # Definir archivo R2
    R2=${READS_DIR}/${SAMPLE}_2_val_2.fq.gz

    # Verificar existencia de R2
    if [ ! -f "$R2" ]; then
        echo "ERROR: No se encontró R2 para $SAMPLE"
        continue
    fi

    echo ""
    echo "---------------------------------------"
    echo "Procesando muestra: $SAMPLE"
    echo "R1: $R1"
    echo "R2: $R2"

    # Se ejecuta Salmon
    salmon quant \
        -i $INDEX_DIR \
        -l A \
        -1 $R1 \
        -2 $R2 \
        -p $THREADS \
        --validateMappings \
        --gcBias \
        --seqBias \
        -o ${OUTPUT_DIR}/${SAMPLE}_quant

    echo "Muestra completada: $SAMPLE"

done

echo ""
echo "======================================="
echo " CUANTIFICACIÓN FINALIZADA "
echo "======================================="
echo "Resultados guardados en:"
echo "$OUTPUT_DIR"
