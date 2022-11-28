#!/bin/bash

# Busqueda de conservacion por sintenia en lncRNAs
# PARTE 2.- Generacion de matrices
# USAGE= ./Sinteny_conservation_lncRNAs_make_matriz_count_p2.sh At_team_lncRNAs_intergenicos_GENES_IDs

# AUTHOR: Evelia Lorena Coss Navarrete
# AIM 1: Generar un resumen sobre los lcnRNAs conservados por sintenia en las especies analizadas.
# AIM 2: Analisis de filostratigrafia de los lncRNAs en cada filostratum.
# https://github.com/EveliaCoss

for i in  $(cat $1);
  #for i in  $1;
  do
  echo $i
  echo "--------------"

# PARTE 1.- Obtener la informacion de Arabidopsis thaliana y de las otras especies en forma de matriz
# lncRNAs conservados por sintenia, ubicacion genomica y numero de genes codificantes que lo flanquean.

    for file in $(ls ./Sinteny_block_genes/At_lncRNA_${i}_syntenic_block_GENES_sorted.bed | cat );
    do
    echo $file
    segment=$(awk '{print $1}' $file | sort -V | uniq)
    number_genes=$(awk '{print $2}' $file | grep -v $i | wc -l)
   #Tabla = Nombre del lncRNA ($1), Especie ($2), segmento genomico ($3) y Numero de genes vecinos ($4)
       awk -v lncRNA=$i -v name=Arabidopsis_thaliana -v segment=$segment -v genes=$number_genes 'BEGIN {
    # print array elements
    print  ""lncRNA"\t" ""name"\t" ""segment"\t" ""genes"\t" }' >> At_lncRNA_summary_localization.tsv
  done

  # Obtener informacion sobre las otras especies
r=("Brassica_oleracea" "Brassica_rapa" "Thellungiella_parvula" "Eutrema_salsugineum" "Arabis_alpina" "Camelina_sativa" "Capsella_rubella" "Arabidopsis_lyrata" "Arabidopsis_halleri" "Aethionema_arabicum" "Boechera_stricta")

 for specie in "${r[@]}"; do
   for lsegment in $(cat ./Genomic_segment_species/At_${i}_${specie}_segment_genomic); do
     for file in $(ls ./Sinteny_block_genes/At_lncRNA_${i}_${specie}_syntenic_block_GENES_${lsegment}_rename.bed | cat ); do
      segment=$(awk '{print $1}' $file | sort -V | uniq)
     number_genes=$(awk '{print $2}' $file | grep -v $i | wc -l)
          #Tabla = Nombre del lncRNA ($1), Especie ($2), segmento genomico ($3) y Numero de genes vecinos ($4)
      awk -v lncRNA=$i -v name=$specie -v segment=$segment -v genes=$number_genes 'BEGIN {
      # print array elements
      print  ""lncRNA"\t" ""name"\t" ""segment"\t" ""genes"\t" }' >> At_lncRNA_species_summary_localization.tsv
     done
    done
  done
 cat At_lncRNA_summary_localization.tsv At_lncRNA_species_summary_localization.tsv | sort -V > At_lncRNA_species_summary_sinteny_analysis.tsv

# USAGE= ./Sinteny_conservation_lncRNAs_make_matriz_count_p2.sh  At_lncRNA_species_summary_sinteny_analysis.tsv

# PARTE 2.- Cuantos lncRNAs se encuentran conservados por sintenia en cada especie?
#r=("Brassica_oleracea" "Brassica_rapa" "Thellungiella_parvula" "Eutrema_salsugineum" "Arabis_alpina" "Camelina_sativa" "Capsella_rubella" "Arabidopsis_lyrata" "Arabidopsis_halleri" "Aethionema_arabicum" "Boechera_stricta")

  for specie in "${r[@]}"; do
   echo $specie
   echo "--------------"
  file=$(ls At_lncRNA_species_summary_sinteny_analysis.tsv)
    # Busqueda global por especies.
    # cuantos lncRNAs se conservan por sintenia en las otras especies? (con respecto a A. thaliana)
   number_genes=$(grep $specie $file | awk '{print $1}' | sort -V | uniq | wc -l)
    # Generar tabla
   awk -v name=$specie -v genes=$number_genes 'BEGIN {
    # print array elements
    print ""name"\t" ""genes"\t" }' >> At_lncRNA_species_summary_Sinteny_results_GLOBAL.tsv

    # lncRNAs conservados por sintenia en cada especie 
   grep $specie $file | awk '{print $1}' | sort -V | uniq > At_lncRNAs_Sinteny_results_${specie}_all_lncRNAs.tsv
 done
done

# USAGE= ./Sinteny_conservation_lncRNAs_make_matriz_count_p2.sh

# Comparacion nodo 5 = At, Al y Ah
#  Specie1=Arabidopsis_lyrata
#  Specie2=Arabidopsis_halleri
#  number_genes=$(cat At_lncRNAs_Sinteny_results_${Specie1}_all_lncRNAs.tsv At_lncRNAs_Sinteny_results_${Specie2}_all_lncRNAs.tsv | sort -V | uniq | wc -l)
#  node=node5
#  awk -v name=$node -v genes=$number_genes 'BEGIN {
    # print array elements
#  print ""name"\t" ""genes"\t" }' >> At_lncRNA_species_Phylostrata_results.tsv

# Comparacion nodo 4 = Cs y Cr
#  Specie3=Camelina_sativa
#  Specie4=Capsella_rubella
#  number_genes=$(cat At_lncRNAs_Sinteny_results_${Specie3}_all_lncRNAs.tsv At_lncRNAs_Sinteny_results_${Specie4}_all_lncRNAs.tsv | sort -V | uniq | wc -l)
#  node=node4
#  awk -v name=$node -v genes=$number_genes 'BEGIN {
    # print array elements
#  print ""name"\t" ""genes"\t" }' >> At_lncRNA_species_Phylostrata_results.tsv

  # Comparacion nodo 3 = Bs
#  Specie5=Boechera_stricta
#  number_genes=$(cat At_lncRNAs_Sinteny_results_${Specie5}_all_lncRNAs.tsv | sort -V | uniq | wc -l)
#  node=node3
#  awk -v name=$node -v genes=$number_genes 'BEGIN {
    # print array elements
#  print ""name"\t" ""genes"\t" }' >> At_lncRNA_species_Phylostrata_results.tsv

# Comparacion nodo 2 =  Aa, Es, Tp, Br y Bo
#  Specie6=Arabis_alpina
#  Specie7=Brassica_oleracea
#  Specie8=Brassica_rapa
#  Specie9=Thellungiella_parvula
#  Specie10=Eutrema_salsugineum
#  number_genes=$(cat At_lncRNAs_Sinteny_results_${Specie6}_all_lncRNAs.tsv At_lncRNAs_Sinteny_results_${Specie7}_all_lncRNAs.tsv \
#  At_lncRNAs_Sinteny_results_${Specie8}_all_lncRNAs.tsv At_lncRNAs_Sinteny_results_${Specie9}_all_lncRNAs.tsv \
#  At_lncRNAs_Sinteny_results_${Specie10}_all_lncRNAs.tsv | sort -V | uniq | wc -l)
#  node=node2
#  awk -v name=$node -v genes=$number_genes 'BEGIN {
    # print array elements
#  print ""name"\t" ""genes"\t" }' >> At_lncRNA_species_Phylostrata_results.tsv

  # Comparacion nodo 1 = Aa
#  Specie11=Aethionema_arabicum
#  number_genes=$(cat At_lncRNAs_Sinteny_results_${Specie11}_all_lncRNAs.tsv | sort -V | uniq | wc -l)
#  node=node3
#  awk -v name=$node -v genes=$number_genes 'BEGIN {
#    # print array elements
#  print ""name"\t" ""genes"\t" }' >> At_lncRNA_species_Phylostrata_results.tsv

