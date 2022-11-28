#!/bin/bash

# Obtener los motivos conservados entre las diferentes especies
# USAGE= ./lncRNAs_FASTA_filter_v2.sh At_team_lncRNAs_intergenicos_GENES_IDs
# sed -i 's/\r//' lncRNAs_FASTA_filter_v2.sh

# AUTHOR: Evelia Lorena Coss Navarrete
# https://github.com/EveliaCoss
mkdir FASTA_specie_specific
mkdir FASTA_filter

for i in  $(cat $1);
# AT1G03106
#for i in  $1;
do
echo $i
echo "--------------"

# PARTE 1.- Seleccionar por tamanos
r=("Brassica_oleracea" "Brassica_rapa" "Thellungiella_parvula" "Eutrema_salsugineum" "Arabis_alpina" "Camelina_sativa" "Capsella_rubella" "Arabidopsis_lyrata" "Arabidopsis_halleri" "Arabidopsis_thaliana" "Aethionema_arabicum" "Boechera_stricta")
#r=("Arabidopsis_lyrata" "Arabidopsis_halleri" "Arabidopsis_thaliana")

	for specie in "${r[@]}"; do
		# Seleccionar la primera isoforma de cada gen
		file=$(ls ./FASTA_files/${i}* | sort -V | head -1)
		# Analizar los tamanos de cada secuencia con infoseq
		infoseq $file -only -name -type -length | sort -k 3n > At_transcripts_length_contador
		# Secuencias mayores a 10 nt y extraer los IDs
		awk '$3 >= 10 {print}' At_transcripts_length_contador | awk '{print $1}' | tail -n +2 | sort -V  > At_transcripts_IDs_species_contador
		# Eliminar species con secuencias pequenas
		awk 'BEGIN{OFS="\t"} {if(/^>/){s=""; h=$1}else{s=s$1}; seqs[h]=s;} END{for (h in seqs) print h, seqs[h]} ' $file | grep -Ff At_transcripts_IDs_species_contador | sed 's/\t/\n/' > ./FASTA_filter/${i}.fa
	done

# Mover los FASTA de especie especifica a otra carpeta
   # Solo se encuentran en Arabidopsis thaliana
      number_species=$(grep ">" ./FASTA_filter/${i}.fa | wc -l)
      # 
    if [[ "$number_species" == 1 ]] ; then
      echo "Solo se encuentra en Arabidopsis_thaliana"
      mv ./FASTA_filter/${i}.fa ./FASTA_specie_specific/
    else
      echo "Se encuentra en otras especies"
    fi
done
