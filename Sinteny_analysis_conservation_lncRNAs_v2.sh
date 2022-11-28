#!/bin/bash

# Busqueda de conservacion por sintenia en lncRNAs
# PARTE 1.- Generacion de bloques sintenicos por especie.
# USAGE= ./Sinteny_analysis_conservation_lncRNAs_v1.sh At_team_lncRNAs_intergenicos_GENES_IDs

# AUTHOR: Evelia Lorena Coss Navarrete
# AIM: Analisis de conservacion por sintenia en lncRNAs intergenicos (lincRNAs).
# https://github.com/EveliaCoss

mkdir At_Sinteny_block_genes
mkdir At_Genomic_segment_species
mkdir At_Neighbors_At_and_species_pcoding_10kb_BED

for i in  $(cat $1);
#for i in  $1; # pruebas
do
echo $i
echo "--------------"

#---Arabidopsis thaliana----

# 1.- Busqueda de las coordenadas del lncRNA
grep $i At_team_lncRNAs_clean_v2.bed > At_lncRNA_contador.bed

# 2.- Localizar las proteinas up and downstream de cada lncRNA
bedtools window -a At_lncRNA_contador.bed -b Araport11_pcoding_transcripts.bed -w 10000 > At_lncRNA_neighbors_pcoding_10kb.bed 

#(-w te da up y downstream)

# 3.- Extraer los nombres de los genes codificantes vecinos al lncRNA
awk 'BEGIN {FS=OFS="\t"} {print $16}' At_lncRNA_neighbors_pcoding_10kb.bed | sort -V | uniq > At_lncRNA_neighbors_pcoding_10kb_IDs

# 4.- Obtener y ordendar las coordenadas de los genes vecinos
grep -Ff At_lncRNA_neighbors_pcoding_10kb_IDs Araport11_pcoding_transcripts.bed | sort -k1,1 -k2,2n > At_lncRNAs_neighbors_pcoding_10kb_coord.bed

# 5.- Unir archivos de coordenadas y ordenarlos
cat At_lncRNA_contador.bed At_lncRNAs_neighbors_pcoding_10kb_coord.bed | sort -k1,1 -k2,2n  > ./At_Sinteny_block_genes/At_lncRNA_${i}_syntenic_block.bed

# 6.- Archivo condensado de las coordenadas para cargar en R
awk 'BEGIN {FS=OFS="\t"} {$4=substr($4, 1, length($4)-2); print}' ./At_Sinteny_block_genes/At_lncRNA_${i}_syntenic_block.bed > ./At_Sinteny_block_genes/At_lncRNA_${i}_syntenic_isoforms_condense.bed

sort -k4,4 ./At_Sinteny_block_genes/At_lncRNA_${i}_syntenic_isoforms_condense.bed | groupBy -g 1,4 -c 4,2,3 -o count,min,max | awk -v OFS='\t' '{print $1, $4, $5, $2, $3}' | sed 's/Chr/Arabidopsis_thaliana_chr/g' > ./At_Sinteny_block_genes/At_lncRNA_${i}_syntenic_block_GENES.bed

# 7.- Generar archivo para R
awk 'BEGIN {FS=OFS="\t"} {print $1, $4, $2, $3}' ./At_Sinteny_block_genes/At_lncRNA_${i}_syntenic_block_GENES.bed > ./At_Sinteny_block_genes/At_lncRNA_${i}_syntenic_block_GENES_sorted.bed

# Mover archivos
# > Bloques sintenicos por especie, coordenadas. 
# Se puede usar para el UCSC con referencia en A. thaliana
mv ./At_Sinteny_block_genes/At_lncRNA_${i}_syntenic_block.bed ./At_Neighbors_At_and_species_pcoding_10kb_BED

#--- Otras especies-----

# Busqueda de las coordenadas del lncRNA
# Resultado de hallifover
# r=("Brassica_napus" "Brassica_oleracea" "Brassica_rapa" "Eutrema_parvulum" "Eutrema_salsugineum" "Arabis_alpina" "Boechera_stricta" "Camelina_sativa" "Capsella_rubella" "Arabidopsis_lyrata" "Arabidopsis_halleri" "Aethionema_arabicum")

# Especies de las que tengo orthofinder
r=("Brassica_oleracea" "Brassica_rapa" "Thellungiella_parvula" "Eutrema_salsugineum" "Arabis_alpina" "Camelina_sativa" "Capsella_rubella" "Arabidopsis_lyrata" "Arabidopsis_halleri" "Aethionema_arabicum" "Boechera_stricta")

  for specie in "${r[@]}"; do
    # 1.- Obtener las coordenadas del lncRNA en cada especie
    grep $i ./HalLiftover_lncRNAs_results/${specie}_At_team_lncRNAs.bed > ./At_Sinteny_block_genes/At_${i}_lncRNA_${specie}_transcripts.bed 
    # 2.- Obtener las coordenadas de los genes codificantes vecinos
    grep -Ff At_lncRNA_neighbors_pcoding_10kb_IDs ../OrthoFinder_files/Pcoding_Orthologos_Position_species/${specie}_At_protein_genes_ORTHO.bed > ./At_Sinteny_block_genes/At_${i}_${specie}_neighbors_pcoding_10kb.bed
    # ---En que segmentos de cada genoma se encuentran?---
    # Ubicacion del lncRNA
    grep -w $i ./At_Sinteny_block_genes/At_${i}_lncRNA_${specie}_transcripts.bed  | awk '{print $1}' | sort -V | uniq > ./At_Sinteny_block_genes/At_${i}_${specie}_segment_genomic

    # Ubicacion de los genes codificantes
    for pcoding in $(cat At_lncRNA_neighbors_pcoding_10kb_IDs); do
      grep -w $pcoding ./At_Sinteny_block_genes/At_${i}_${specie}_neighbors_pcoding_10kb.bed | awk '{print $1}' | sort -V | uniq > ./At_Sinteny_block_genes/At_${i}_${specie}_${pcoding}_segment_genomic
    done
    #----------
    # 3.- Unir archivos de coordenadas y ordenarlos
      cat ./At_Sinteny_block_genes/At_${i}_lncRNA_${specie}_transcripts.bed ./At_Sinteny_block_genes/At_${i}_${specie}_neighbors_pcoding_10kb.bed | sort -k1,1 -k2,2n > ./At_Sinteny_block_genes/At_${i}_${specie}_syntenic_block.bed
    # 4.- Archivo condensado de las coordenadas para cargar en R
      awk 'BEGIN {FS=OFS="\t"} {$4=substr($4, 1, length($4)-2); print}' ./At_Sinteny_block_genes/At_${i}_${specie}_syntenic_block.bed > ./At_Sinteny_block_genes/At_${i}_${specie}_syntenic_isoforms_condense.bed
    # 5.- Seleccionar por segmentos, donde se encuentra el lncRNA
    for lsegment in $(cat ./At_Sinteny_block_genes/At_${i}_${specie}_segment_genomic); do
      # Seleccion por cada segmento
      awk -v SEGMENT="$lsegment" 'BEGIN {FS=OFS="\t"} {if($1==SEGMENT) print $0}' ./At_Sinteny_block_genes/At_${i}_${specie}_syntenic_isoforms_condense.bed > ./At_Sinteny_block_genes/At_${i}_${specie}_syntenic_selected_${lsegment}.bed
      # Condensar los archivos
      sort -k4,4 ./At_Sinteny_block_genes/At_${i}_${specie}_syntenic_selected_${lsegment}.bed | groupBy -g 1,4 -c 4,2,3 -o count,min,max | awk -v OFS='\t' '{print $1, $4, $5, $2, $3}' > ./At_Sinteny_block_genes/At_lncRNA_${i}_${specie}_syntenic_block_GENES_${lsegment}.bed
      # Archivo para R de cada segmento
      awk 'BEGIN {FS=OFS="\t"} {print $1, $4, $2, $3}' ./At_Sinteny_block_genes/At_lncRNA_${i}_${specie}_syntenic_block_GENES_${lsegment}.bed > ./At_Sinteny_block_genes/At_lncRNA_${i}_${specie}_syntenic_block_GENES_${lsegment}_sorted.bed
      # Renombrar la columna uno con la especie y el segmento genomico
      # Archivo para cargar en R
      name="${specie}_${lsegment}"
      echo "Specie_segmento: $name"
      sed "s|$lsegment|$name|g" ./At_Sinteny_block_genes/At_lncRNA_${i}_${specie}_syntenic_block_GENES_${lsegment}_sorted.bed > ./At_Sinteny_block_genes/At_lncRNA_${i}_${specie}_syntenic_block_GENES_${lsegment}_rename.bed
    done
  done

# Mover archivos
# > Segmentos genomicos de los genes codificantes y el lncRNA de interes, por especie
mv ./At_Sinteny_block_genes/*_segment_genomic ./At_Genomic_segment_species
# > Bloques sintenicos por especie, coordenadas. 
# Se puede usar para el UCSC con referencia en A. thaliana
mv ./At_Sinteny_block_genes/*_syntenic_block.bed ./At_Neighbors_At_and_species_pcoding_10kb_BED

done

#rm -rf At_Genomic_segment_species At_Neighbors_At_and_species_pcoding_10kb_BED  At_Sinteny_block_genes

