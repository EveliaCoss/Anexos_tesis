#!/bin/bash

# Obtener BED de los lncRNAs conservados por sintenia
# PARTE 1.- Extraer los BEDs de los lncRNAs conservados por sintenia
# USAGE= ./Ext_Sinteny_BED_files.sh
# sed -i 's/\r//' Ext_Sinteny_BED_files.sh

# AUTHOR: Evelia Lorena Coss Navarrete
# https://github.com/EveliaCoss

mkdir Sinteny_BEDs

r=("Brassica_oleracea" "Brassica_rapa" "Thellungiella_parvula" "Eutrema_salsugineum" "Arabis_alpina" "Camelina_sativa" "Capsella_rubella" "Arabidopsis_lyrata" "Arabidopsis_halleri" "Aethionema_arabicum" "Boechera_stricta")

 for specie in "${r[@]}"; do
	grep -Ff At_lncRNAs_Sinteny_results_${specie}_all_lncRNAs.tsv ./HalLiftover_lncRNAs_results/${specie}_At_team_lncRNAs.bed > ./Sinteny_BEDs/${specie}.bed
done

