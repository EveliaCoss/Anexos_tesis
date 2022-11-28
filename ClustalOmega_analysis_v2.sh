#!/bin/bash

# PARTE 1.- Porcentaje de identidad con Clustal Omega
# USAGE= ./ClustalOmega_analysis_v2.sh
# sed -i 's/\r//' ClustalOmega_analysis_v2.sh

# AUTHOR: Evelia Lorena Coss Navarrete
# https://github.com/EveliaCoss

#./AtlincRNA_allspecies_FASTA_files/*_AtlincRNAvsSpTranscript.fasta
#./AtlincRNA_allspecies_FASTA_files/*_AtlincRNAvsSplincRNA.fasta

mkdir ./AtlincRNA_SplincRNA_ClustalOmega_analysis

#------------------
# Obtener el PID de los AtlincRNA conservados e identificados en las otras especies (SplincRNAs) 
#------------------

if [[ -f ./AtlincRNA_SplincRNA_ClustalOmega_analysis/AT1G04087.1_PID_AtlincRNAvsSpTranscript.distmat ]] 
then
  echo "Analisis con ClustalOmega -AtlincRNA-SplincRNA -check"
else

echo  "--->----<-----"
echo "-----> Analisis con ClustalOmega -AtlincRNA-SplincRNA <-----"

file=$(ls ./AtlincRNA_allspecies_FASTA_files/*_AtlincRNAvsSplincRNA.fasta)
	for fi in $(echo $file)
	#for i in  $(cBo $1);
	#Bo3G09922
	#for i in  $1;
	do
	echo $fi
	echo "--------------"
	name=$(echo $fi | sed 's/_AtlincRNAvsSplincRNA.fasta//;s/AtlincRNA_allspecies_FASTA_files//;s/\.//;s/\///g') # renombrar
	clustalo --infile  $fi --threads 8 -t DNA --full --percent-id --distmat-out=./AtlincRNA_SplincRNA_ClustalOmega_analysis/${name}_PID_AtlincRNAvsSplincRNA.distmat --output-order=tree-order --outfmt=clu -o ./AtlincRNA_SplincRNA_ClustalOmega_analysis/${name}_clustalo_alignment_AtlincRNAvsSplincRNA.txt --force

# [--infile]              Input
# [-t]                    Tipo de secuencia
# [--full --percent-id --distmBo-out=] Parametros para el porcentaje de identidad
# [--output-order=]       Orden de los archivos
# [--outfmt=]   Tipo de salida del alineamiento, clustal
# [--force]    Sobreescribir salida
# [-o]         salida

	done
fi

#------------------
# Obtener el PID de los AtlincRNA conservados e identificados por posicion (alineamineto) en las otras especies (SpTranscript) 
#------------------


file=$(ls ./AtlincRNA_allspeciesOther_ALLFASTA_files/*_AtlincRNAvsSpTranscript.fasta)
	for fi in $(echo $file)
	#for i in  $(cBo $1);
	#Bo3G09922
	#for i in  $1;
	do
	echo $fi
	echo "--------------"
	name=$(echo $fi | sed 's/_AtlincRNAvsSpTranscript.fasta//;s/AtlincRNA_allspecies_FASTA_files//;s/\.//;s/\///g') # renombrar
	clustalo --infile  $fi --threads 8 -t DNA --full --percent-id --distmat-out=./AtlincRNA_SplincRNA_ClustalOmega_analysis/${name}_PID_AtlincRNAvsSpTranscript.distmat --output-order=tree-order --outfmt=clu -o ./AtlincRNA_SplincRNA_ClustalOmega_analysis/${name}_clustalo_alignment_AtlincRNAvsSpTranscript.txt --force

# [--infile]              Input
# [-t]                    Tipo de secuencia
# [--full --percent-id --distmBo-out=] Parametros para el porcentaje de identidad
# [--output-order=]       Orden de los archivos
# [--outfmt=]   Tipo de salida del alineamiento, clustal
# [--force]    Sobreescribir salida
# [-o]         salida

	done