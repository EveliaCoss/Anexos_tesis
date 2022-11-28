# 03_Conservacion por Sintenia en lncRNAs

## Construccion de los bloques sintenicos

El script `Sinteny_analysis_conservation_lncRNAs_v1.sh` emplea los archivos:

```
Carpeta
 |At_team_lncRNAs_clean_v2.bed           # BED de los lncRNAs
 |Araport11_pcoding_transcripts.bed      # BED de los genes codificantes
 |At_team_lncRNAs_clasificados_IDs_GENES # IDs de lncRNAs intergenicos (lincRNAs)
 |-HalLiftover_lncRNAs_results/          # BED de los lncRNAs conservados por posicion, resultado de HalLiftover
   |-${specie}_At_team_lncRNAs.bed


# Ademas emplea los resultados de los genes codificantes conservados por posicion y ortologia.
 ../OrthoFinder_files/Pcoding_Orthologos_Position_species

 # USAGE:
./Sinteny_analysis_conservation_lncRNAs_v1.sh At_team_lncRNAs_intergenicos_GENES_IDs

```

Ademas, el script emplea el programa `Bedtools window`, por lo que es necesario que instalen bedtools y lo incorporen en el PATH en caso de que lo corran en su computadora. En mazorka el modulo es `bedtools2/2.25.0`.

Genera de salida los siguientes archivos:

* `At_lncRNA_contador.bed`: Coordenadas del lncRNA analizado en ese momento (cambia por cada ciclo).
* `At_lncRNA_neighbors_pcoding_10kb.bed`: Coordenadas de los genes vecinos al lncRNA (cambia por cada ciclo).
* `At_lncRNA_neighbors_pcoding_10kb_IDs`: IDs de los genes vecinos al lncRNA (cambia por cada ciclo).
* `At_lncRNAs_neighbors_pcoding_10kb_coord.bed`: Coordenadas del lncRNA y de los genes vecinos (cambia por cada ciclo).

Dentro de la carpeta `Sinteny_block_genes` se encuentran los archivos:

* `At_lncRNA_${i}_syntenic_isoforms_condense.bed`: Bloques del vecindario genetico, en el geneID se quitaron las isoformas.
* `At_lncRNA_${i}_syntenic_block_GENES.bed`: Archivo con cinco columnas, contiene: el numero del cromosoma en *A. thaliana* ($1), coordenada inicial ($2) y final ($3), nombre de los genes ($4), asi como el numero de fragmentos encontrados para cada gen ($5).
* `At_lncRNA_${i}_syntenic_block_GENES_sorted.bed`: Archivo de cuatro columnas, contiene: el numero del cromosoma en *A. thaliana* ($1), nombre de los genes ($2), coordenada inicial ($3) y final ($4). Este archivo puede emplearse para ser visualizado en `R` con `ggeenes`.
* `At_${i}_lncRNA_${specie}_transcripts.bed `: Coordenadas del lncRNA en cada especie analizada.
* `At_${i}_${specie}_neighbors_pcoding_10kb.bed`: Coordenadas de los genes codificantes vecinos al lncRNA para cada especie.
* `At_${i}_${specie}_segment_genomic`: Segmentos en donde se localiza el lncRNA en cada especie.
* `At_${i}_${specie}_syntenic_isoforms_condense.bed`: Bloques del vecindario genetico, en el geneID se quitaron las isoformas para cada especie.
* `At_${i}_${specie}_syntenic_selected_${lsegment}.bed`: Coordenadas donde se localiza al lncRNA con los genes vecinos en cada especie.
* `At_${i}_${specie}_syntenic_isoforms_condense.bed`: Bloques del vecindario genetico, en el geneID se quitaron las isoformas en cada especie y lncRNA.
* `At_lncRNA_${i}_${specie}_syntenic_block_GENES_${lsegment}.bed`: Archivo con cinco columnas, contiene: el numero del cromosoma en cada especie ($1), coordenada inicial ($2) y final ($3), nombre de los genes ($4), asi como el numero de fragmentos encontrados para cada gen ($5).
* `At_lncRNA_${i}_${specie}_syntenic_block_GENES_${lsegment}_sorted.bed`: chivo de cuatro columnas, contiene: el numero del cromosoma en *A. thaliana* ($1), nombre de los genes ($2), coordenada inicial ($3) y final ($4).
* `At_lncRNA_${i}_${specie}_syntenic_block_GENES_${lsegment}_rename.bed`: Archivo de cuatro columnas, contiene: el numero del cromosoma en cada especie ($1), nombre de los genes ($2), coordenada inicial ($3) y final ($4). Se renombro el cromosoma con la especie y el segmento en donde se localizo el lncRNA con sus genes vecinos. Este archivo puede emplearse para ser visualizado en `R` con `ggeenes`.

Dentro de la carpeta `Neighbors_At_and_species_pcoding_10kb_BED` se encuentran los archivos:

* `At_lncRNA_${i}_syntenic_block.bed`: Bloques del vecindario genetico para cada lncRNA en *A. thaliana*.
* `At_${i}_${specie}_syntenic_block.bed`: Coordenadas de los genes codificantes y el lncRNA en cada segment genomico para cada especie y lncRNA.

Dentro de la carpeta `Genomic_segment_species` se encuentran los archivos:
* `At_${i}_${specie}_${pcoding}_segment_genomic`: Coordenadas de los genes codificantes en cada segmento genomico.


## Interpretacion de los resultados de la conservacion por sintenia

El script `Sinteny_conservation_lncRNAs_make_matriz_count_p2.sh` genera los siguientes archivos:

* `At_lncRNA_summary_localization.tsv`: Tabla con los lncRNAs en *A. thaliana* y el numero de genes codificantes que lo flanquean en la ultima columna.
*  `At_lncRNA_species_summary_localization.tsv`: Tabla de los lncRNAs localizados por sintenia en cada especie y el numero de genes codificantes que lo flanquean en la ultima columna.
* `At_lncRNA_species_summary_sinteny_analysis.tsv`: Tabla con los resultados en conjunto de los lncRNAs conservados por sintenia en todas las especies.
* `At_lncRNAs_Sinteny_results_${specie}_all_lncRNAs.tsv`: IDs se los lncRNAs conservados por sintenia, cada archivo es una especie analizada.

genera una tabla de resumen que indica el numero de lncRNAs conservados por sintenia en cada especie. 
