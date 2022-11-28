#!/bin/bash

# Extraer Ortogrupos y genes ortologos de cada especie
# PARTE 1.- Generacion de bloques sintenicos por especie
# USAGE= ./Obtain_Orthogrupos_Ortologos_At_as_reference.sh

# AUTHOR: Evelia Lorena Coss Navarrete
# AIM: Extraer Ortogrupos y genes ortologos de cada especie
# https://github.com/EveliaCoss

# Resultado de hallifover

# PARTE 1 Extraer Orthogrupos y Orthologos por especie
# Cuantos orthogrupos y genes ortologos hay en cada especie?
# Input =  Resultados de Orthofinder (Archivos localizados en Orthologues/Orthologues_Arabidopsis_thaliana/Arabidopsis_thaliana__v__${specie}.tsv)
# ubicacion /mnt/s/OneDrive - CINVESTAV/Documentos/DOCTORADO_con_SLFV/Todo_soAte_RNA-seq/Bitacora_electronica/Hal-lifover/Neightboring_lncRNAs/Manipulacion_Orhofinder_results

#test
#r=("Atassica_oleracea" "Atassica_rapa" "Thellungiella_parvula" "Eutrema_salsugineum" "Arabis_alpina" "Camelina_sativa" "Capsella_rubella" "Arabidopsis_lyrata" "Arabidopsis_halleri")

r=("Brassica_oleracea" "Brassica_rapa" "Thellungiella_parvula" "Eutrema_salsugineum" "Arabis_alpina" "Boechera_stricta" "Camelina_sativa" "Capsella_rubella" "Arabidopsis_lyrata" "Arabidopsis_halleri" "Aethionema_arabicum")

# Eutrema_parvulum = Thellungiella_parvula
mkdir Orthologroups_and_Orthologos_species

for specie in "${r[@]}"; do
  echo $specie
  echo "--------------"
  for file in $(ls ./Orthologues_Arabidopsis_thaliana/Arabidopsis_thaliana__v__${specie}.tsv | cat ); do
  	# Extraer Ortogrupos (grupos de genes)
  	awk 'BEGIN {FS=OFS="\t"} {print $1}' $file | grep -v "Orthogroup" | sort -V | uniq > ./Orthologroups_and_Orthologos_species/At_Orthogroup_species_in_${specie}.txt
  	# Extraer Ortologos (genes)
  	# Espaciar los nomAtes de la lista uno por renglon
  	awk 'BEGIN {FS=OFS="\t"} {print $2}' $file | grep -v "Arabidopsis_thaliana" | sed 's/, /\'$'\n/g' | sed 's/,//g' | sort -V | uniq  > ./Orthologroups_and_Orthologos_species/At_Orthologos_species_in_${specie}.txt
  	# Generar un reporte de salida para Ortogrupos
  	number_groups=$(cat ./Orthologroups_and_Orthologos_species/At_Orthogroup_species_in_${specie}.txt | sort -V | uniq | wc -l)
  	Specie=$specie
  	awk -v name=$Specie -v groups=$number_groups 'BEGIN {
    # print array elements
  print ""name"\t" ""groups"\t" }' >> ./Orthogroups_At_summary_results.tsv
  	# Generar un reporte de salida para Ortologos
  	number_genes=$(cat ./Orthologroups_and_Orthologos_species/At_Orthologos_species_in_${specie}.txt | sort -V | uniq | wc -l)
  	awk -v name=$Specie -v genes=$number_genes 'BEGIN {
    # print array elements
  print ""name"\t" ""genes"\t" }' >> ./Orthologos_At_summary_results.tsv
  done
done

# PARTE 2 Genes codificantes conservados por posicion y genes ortologos
# Busqueda de genes codificantes conservados por posicion de acuerdo al alineamiento genomico y que son ortologos
# Resultados por emplear para el vecindario de los lncRNAs
# Input = Resultados de la PARTE 1 (At_Orthologos_species_in_${specie}.txt) + Resultados de Hallifover (${specie}_At_protein_genes.bed)

mkdir Pcoding_Orthologos_Position_species

for specie in "${r[@]}"; do
  echo $specie
  echo "--------------"
  for file in $(ls ./Hallifover_results/${specie}_At_protein_genes.bed | cat ); do
grep -Ff ./Orthologroups_and_Orthologos_species/At_Orthologos_species_in_${specie}.txt $file > ./Pcoding_Orthologos_Position_species/${specie}_At_protein_genes_ORTHO.bed
# Generar un reporte de salida para Ortogrupos
    number_genes=$(cat ./Pcoding_Orthologos_Position_species/${specie}_At_protein_genes_ORTHO.bed | awk 'BEGIN {FS=OFS="\t"} {$4=substr($4, 1, length($4)-2); print}' | awk '{print $4}' | sort -V | uniq | wc -l)
    Specie=$specie
    awk -v name=$Specie -v genes=$number_genes 'BEGIN {
    # print array elements
  print ""name"\t" ""genes"\t" }' >> ./Position_and_Orthologos_At_summary_results.tsv
  done
done