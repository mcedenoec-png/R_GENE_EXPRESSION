#!/bin/bash

# Archivo con IDs SRR
LISTA="SRR_Acc_List.txt"

# Carpeta de salida
mkdir -p SRA_downloads

# Descargar cada SRR
while read SRR; do
    echo "Descargando $SRR ..."
    
    # Descargar archivo .sra
    prefetch "$SRR" --output-directory SRA_downloads

    # Convertir a FASTQ (paired-end)
    fasterq-dump "SRA_downloads/$SRR/$SRR.sra" \
        --split-files \
        --threads 4 \
        --temp ~/tmp_sra \
        --outdir SRA_downloads

done < "$LISTA"

echo "Descarga y conversión completadas."
