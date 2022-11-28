# 01_Localizacion de genes codificantes y lncRNAs mediante alineamiento genomico

El alineamiento genomico generado por Jose Antonio Corona-Gomez presenta solo la numeracion en el cromosoma (1-5), sin Chr o chr. Por lo que, se debe eliminar de la columna de coordenadas.

## Busqueda de proteinas conservadas por posicion

El archivo `Araport11_pcoding_transcripts_sorted.bed` son las coordenadas BED de los 27,535 genes de proteinas presentes en *Arabidopsis thaliana*.

Primero se elimina el `Chr` presente en el nombre de los cromosomas.

```
sed 's/Chr//g' Araport11_pcoding_transcripts_sorted.bed > Araport11_pcoding_transcripts_sorted_ed.bed 

```

El script `Halliftover_allspecies_proteins.sh` emplea la funcion `halLiftover` para obtener los genes codificantes conservados por posicion en otros genomas.

```
# USAGE
# halLiftover 16WGA-brasicaceas.hal Especie_referencia BED Especie_analisis output

halLiftover 16WGA-brasicaceas.hal Arabidopsis_thaliana Araport11_pcoding_transcripts_sorted_ed.bed $i ./Species_New_181021/${i}_At_protein_genes.bed

# 16WGA-brasicaceas.hal: Archivo comprimido del alineamiento genomico generado previamente con Haltools.
# Especie_referencia: Especie de referencia presente en el alineamiento genomico.
# BED: Archivo BED de los genes codificantes de la especie de referencia.
# Especie_analisis: Especie con la que se quiere comparar para obtener las coordenadas.
# output: Archivo de salida de las coordenadas del gen codificante en la especie analizada.
```

Las especies analizadas fueron las mismas que en el alineamiento senaladas en el script como `r`.

## Busqueda de lncRNAs conservados por posicion

Modificar el archivo BED.

```
sed 's/Chr//g' At_team_lncRNAs_clean_v2.bed > At_team_lncRNAs_ed.bed 
```

El script para los lcnRNAs fue `Halliftover_sintenia_lncRNAs.sh`.

```
halLiftover 16WGA-brasicaceas.hal Arabidopsis_thaliana At_team_lncRNAs_ed.bed $i ./Species_NewlncRNAs_181021/${i}_At_team_lncRNAs.bed
```

Se analizaron las mismas especies.