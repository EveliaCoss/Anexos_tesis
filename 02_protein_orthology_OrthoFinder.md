# 02_Busqueda de proteinas Ortologas con Orthofinder

Se emplearon las siguientes guias:

a) [Step-by-step Orthofinder tutorials](https://davidemms.github.io/menu/tutorials.html).

b) En el tutorial **Step-by-step Orthofinder tutorials** en la seccion [Running an example Orthofinder analysis](https://davidemms.github.io/orthofinder_tutorials/running-an-example-orthofinder-analysis.html). 

Se mencionan los siguientes pasos:

##  1.- Crear un folder con el nombre sin espacios

Ejemplo:

path = /mnt/c/Users/eveli/Documents/Orthofinder_brassicas, de acuerdo a la [guia de instalacion](https://davidemms.github.io/orthofinder_tutorials/downloading-and-running-orthofinder.html).

```
# cd /mnt/c/Users/eveli/Documents/Orthofinder_brassicas
mkdir orthofinder_tutorial
```

## 2.- Descargar los archivos (.pep o .faa) de cada especie 

### Ensembl Plants

Localizar la especie y en Download FASTA, seleccionar la carpeta *pep* o entrar en el siguiente enlace: 
http://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/release-50/fasta/

Ejemplo:

```
http://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/release-50/fasta/arabidopsis_thaliana/pep/Arabidopsis_thaliana.TAIR10.pep.all.fa.gz
```

Especies descargadas de esta plataforma: *Arabidopsis thaliana*, *Arabidopsis halleri*, *Arabidopsis lyrata*, *Brassica rapa*, *Brassica oleracea*, *Camelina sativa*, *Arabis alpina* y *Eutrema salsugineum*.

### Phythozome v13

Crear una cuenta de usuario, iniciar sesion para poder usar el Portal de Phytome: 
https://genome.jgi.doe.gov/portal/pages/dynamicOrganismDownload.jsf?organism=Phytozome

Seleccionar la carpeta **PhytozomeV13**, buscar la especie, seleccionar la version del genoma y la carpeta que dice *annotation*. Posteriormente, dar click sobre el archivo 
*protein_primaryTranscriptOnly.fa.gz* y descargar el archivo.

Ejemplo:

```
PhytozomeV13/BrapaFPsc/v1.3/annotation/Phytozome: BrapaFPsc_277_v1.3.protein_primaryTranscriptOnly.fa.gz 
```

En el mismo tutorial de OrthoFinder te mencionan como descargarlo de phytozome, en el [link](https://davidemms.github.io/orthofinder_tutorials/getting-input-data.html).

### Tao, et al. (2018)

Ademas se emplearon los mismos archivo que empleo [Tao, et al (2018)](https://github.com/zhaotao1987/SynNet-Pipeline) para [Syn-MRL](https://www.nature.com/articles/s41467-021-23665-0) de las especies *Aethionema arabicum* y *Tarenaya hassleriana*. Al no localizarse en las otras bases de datos.

Orthofinder puede emplear las extensiones: `".fa”, “.faa”, “.fasta”, “.fas” o “.pep"`. Las secuencias deben encontrarse en aminoacidos.

## 3.- Descomprimir los archivos GZ:

```{bash, eval=F}
gunzip *.gz
```

## 4.- Secuencia primaria

Orthofinder emplea en su analisis solo la **secuencia primaria o una sola isoforma de cada gen** (primary_transcript), por lo que, mediante el script **primary_transcript.py**, se puede extraer estas secuencias para cada especie. Excepto para los archivos que fueron obtenidos de **Phytozome**, ni los archivos de Tao, como son: "Sparvula_574_v2.2.transcript_primaryTranscriptOnly.fa", "Crubella_474_v1.1.protein_primaryTranscriptOnly.fa", y "Bstricta_278_v1.2.cds_primaryTranscriptOnly.fa".

El script **primary_transcript.py** se encuentra dentro de la carpeta de OrthoFinder, por lo que, antes de correrlo en mazorka, es importante procesar las secuencias con este script.

```
# descagar
wget https://github.com/davidemms/OrthoFinder/releases/latest/download/OrthoFinder.tar.gz


# extraer
tar xzvf OrthoFinder.tar.gz
cd OrthoFinder/

# mi script
for f in *fa ; do python /mnt/c/Users/eveli/Documents/Orthofinder_brassicas/orthofinder_tutorial/OrthoFinder/tools/primary_transcript.py $f ; done

```

El script genera una carpeta de salida llamada **primary_transcripts**.

## 5.- Acotar los nombres de los archivos

Renombrar los archivos solo con el nombre de la especie.

```
mv Arabidopsis_thaliana.TAIR10.pep.all.fa Arabidopsis_thaliana.fa 
```

## 6.- Subirlo el archivo y correr el script de OrthoFinder en mazorka

El pipeline normal de Orthofinder es:

1) DIAMOND or MMseqs2 (recommended, although BLAST+ can be used instead)

2) The MCL graph clustering algorithm

3) FastME (The appropriate version for your system, e.g. 'fastme-2.1.5-linux64', should be renamed `fastme`, see instructions below.)

Y con el MSA workflow:

4) Multiple sequence alignment program: MAFFT recommended

5) Tree inference program: FastTree* recommended.

Para subir el archivo se emplea: `rsync -av ./primary_transcripts usuario@mazorka.langebio.cinvestav.mx:"/LUSTRE/usuario/usuario/results/Orthofinder/"`.

El script `orthofinder_pp.sh` muestra las siguientes opciones:

```
orthofinder -f primary_transcripts/ -t 8 -a 8 -S blast -I 1.7 -M msa -A mafft -T fasttree 

#1. -f  solo le dice que lo que sigue es el directorio con tus fastas ordenados por especie.
#2. -S le dice que método para analizar secuencias usar, puede ser blast, diamond o mmseq2, el recomendado es diamond porque es eficiente  y preciso, blast es un poco más preciso pero va a tardar.
#3. El parámetro -I (MCL inflation parameter) puedes ajustarlo, pero 1.7 funciona bien para hacer pruebas, no es recomendable que sea mayor a 2.0. MCL inflation parameter [Default = 1.5]
#En mazorka agregar -t para los threads.
#-T <opt>: Tree inference program opt=fasttree,raxml,iqtree,... user-extendable (requires '-M msa') [Default = fasttree]
# -a number_of_orthofinder_threads = los cuales en plantas se calculan:

# 0.2 Gb como requerimiento minimo en cada threads
# El requerimiento minimo es de 64 GB.
# -I <int>        MCL inflation parameter [Default = 1.5]
```

## 7.- Analisis de los resultados

En la parte de **Exploring Orthofinder'results** de la [guia](https://davidemms.github.io/orthofinder_tutorials/exploring-orthofinders-results.html), en **Quality Control** contamos con:

**OrthoFinder assigned 357029 genes (92.8% of total) to 32524 orthogroups**. Fifty percent of all genes were in orthogroups with 13 or more genes (G50 was 13) and were contained in the largest 8843 orthogroups (O50 was 8843). There were 10924 orthogroups with all species present and 306 of these consisted entirely of single-copy genes.

Para mas informacion se puede observar el archivo **Comparative_Genomics_Statistics/Statistics_Overall.tsv**.

That’s pretty good, in general it’s nice to see at least 80% of your genes assigned to orthogroups. Fewer than this means that you are probably missing orthology relationships that actually exist for some of the remaining genes, poor species sampling is the most likely cause for this. Let’s also check the percentages on a per species basis. There’s a tab-delimited file called **Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv**.

Results Files: Orthogroups Directory (deprecated)

The orthogroups in Phylogenetic_Hierarchical_Orthogroups/ should be used instead. They are identifed using rooted genes trees and are 12%-20% more accurate.

- **Orthogroups.tsv** (deprecated) is a tab separated text file. Each row contains the genes belonging to a single orthogroup. The genes from each orthogroup are organized into columns, one per species. The orthogroups in Phylogenetic_Hierarchical_Orthogroups/N0.tsv should be used instead.

- **Orthogroups_UnassignedGenes.tsv** is a tab separated text file that is identical in format to Orthogroups.csv but contains all of the genes that were not assigned to any orthogroup.

- **Orthogroups.txt** (legacy format) is a second file containing the orthogroups described in the Orthogroups.tsv file but using the OrthoMCL output format.

- **Orthogroups.GeneCount.tsv** is a tab separated text file that is identical in format to Orthogroups.csv but contains counts of the number of genes for each species in each orthogroup.

- **Orthogroups_SingleCopyOrthologues.txt** is a list of orthogroups that contain exactly one gene per species i.e. they contain one-to-one orthologues. They are ideally suited to between-species comparisons and to species tree inference.

```{bash, eval=F}
# count the total number of groups
wc -l Orthogroups.GeneCount.tsv # 32525 Orthogroups.GeneCount.tsv

# count the number of groups with one gene per strain
awk '{if (($2==1)&&($3==1)&&($4==1)&&($5==1)) print}' Orthogroups.GeneCount.tsv | wc -l
8198
```

Especies empleadas en el analisis

```
Species used: 
0: Arabidopsis_halleri.fa
1: Arabidopsis_lyrata.fa
2: Arabidopsis_thaliana.fa
3: Arabis_alpina.fa
4: Brassica_oleracea.fa
5: Brassica_rapa.fa
6: Camelina_sativa.fa
7: Capsella_rubella.fa
8: Eutrema_salsugineum.fa
9: Thellungiella_parvula.fa
10: Aethionema_arabicum.pep
11: Boechera_stricta.fa
12: Tarenaya_hassleriana.pep
```

## 8.- Extraer los genes ortologos y los ortogrupos

El script `Obtain_Orthogrupos_Ortologos_At_as_reference.sh` extraer el numero de genes ortologos y ortogrupos obtenidos con OrthoFinder, para su empleo es necesario contar con dos carpetas: 
* `Hallifover_results`: contiene los resultados de las proteinas conservadas por posicion en el alineamiento genomico en formato BED. 
* `Orthologues_Arabidopsis_thaliana`: contiene los resultados de OrthoFinder en la identificacion de los genes ortologos y a los ortogrupos a los cuales pertenecen para cada especie en formato TSV. Normalmente esta carpeta se localiza en: `Orthologues/Orthologues_Arabidopsis_thaliana/`.

```

OrthoFinder_files
 |-Hallifover_results
 |-Orthologues_Arabidopsis_thaliana

# Hallifover_results/Aethionema_arabicum_At_protein_genes.bed

head Aethionema_arabicum_At_protein_genes.bed
KE150568.1      81573   81612   AT1G01020.2     0       -       81573   81612   0,0,0   1       39      0
KE151742.1      252555  252597  AT1G01020.2     0       -       252555  252597  0,0,0   2       2,37    0,5
KE150568.1      81613   81650   AT1G01020.2     0       -       81613   81650   0,0,0   1       37      0
KE151742.1      252598  252635  AT1G01020.2     0       -       252598  252635  0,0,0   1       37      0
KE150568.1      81744   81828   AT1G01020.2     0       -       81744   81828   0,0,0   1       84      0
KE151742.1      252742  252815  AT1G01020.2     0       -       252742  252815  0,0,0   1       73      0
KE150568.1      81828   81839   AT1G01020.2     0       -       81828   81839   0,0,0   1       11      0
KE151742.1      252815  252823  AT1G01020.2     0       -       252815  252823  0,0,0   1       8       0
KE150568.1      81840   81841   AT1G01020.2     0       -       81840   81841   0,0,0   1       1       0
KE151742.1      252875  252876  AT1G01020.2     0       -       252875  252876  0,0,0   1       1       0

# Orthologues_Arabidopsis_thaliana/Arabidopsis_thaliana__v__Aethionema_arabicum.tsv

head Arabidopsis_thaliana__v__Aethionema_arabicum.tsv
Orthogroup      Arabidopsis_thaliana    Aethionema_arabicum
OG0000000       AT2G13770       AA172G00003, AA110G00005, AA705G00001, AA64G00013, AA771G00002, AA50G00053
OG0000000       AT5G35490       AA168G00040, AA98G00019
OG0000001       AT5G59930       AA97G00041, AA8G00367, AA8G00366
OG0000001       AT3G48400       AA97G00009
OG0000001       AT3G46810       AA16G00021, AA16G00020, AA16G00019, AA16G00001, AA16G00018
OG0000001       AT2G19660, AT2G04680, AT5G37620, AT2G19650      AA44G00395, AA28G00008
OG0000001       AT5G26190, AT5G43030, AT5G43040 AA53G01007
OG0000001       AT1G50190, AT1G55700, AT2G23100, AT5G42280, AT1G55380, AT1G55420, AT1G55430, AT1G55390, AT1G55440      AA76G00006, AA22G00021, AA54G00111, AA87G00273, AA54G00327, AA15G00172
OG0000001       AT2G02680, AT2G02620, AT2G02690, AT2G02700, AT2G02630, AT2G02610, AT2G02640     AA96G00139, AA96G00166, AA102G00156
```

El script se ejecuta fuera de ambas carpetas:

```
./Obtain_Orthogrupos_Ortologos_At_as_reference.sh
```

Se obtienen dos nuevas carpetas y dos archivos de resumen:

* `Pcoding_Orthologos_Position_species`: Coordenadas de los genes codificantes conservados por posicion y que son ortologos con respecto a *Arabidopsis thaliana*.
* `Orthologroups_and_Orthologos_species`: Contiene los IDs de cada especie para ortologos y ortogrupos.
* `Orthogroups_At_summary_results.tsv`: Resumen del numero de ortogrupos detectados en cada especie.
* `Orthologos_At_summary_results.tsv`: Resumen del numero de numero de genes ortologos detectados en cada especie.
* `Position_and_Orthologos_At_summary_results.tsv`: Resumen del numero de genes codificantes conservados por posicion y ortologia.

```
OrthoFinder_files
 |-Hallifover_results
 |-Orthologues_Arabidopsis_thaliana
 |-Pcoding_Orthologos_Position_species
 |-Orthologroups_and_Orthologos_species
 |Orthogroups_At_summary_results.tsv
 |Orthologos_At_summary_results.tsv
 |Position_and_Orthologos_At_summary_results.tsv

# Numero de genes ortologos encontrados en cada especie
cat Orthologos_At_summary_results.tsv
Brassica_oleracea       19961
Brassica_rapa   19753
Thellungiella_parvula   19679
Eutrema_salsugineum     20205
Arabis_alpina   14559
Boechera_stricta        21391
Camelina_sativa 22286
Capsella_rubella        21914
Arabidopsis_lyrata      22562
Arabidopsis_halleri     21908
Aethionema_arabicum     16191

# Numero de ortogrupos encontrados en cada especie
cat Orthogroups_At_summary_results.tsv
Brassica_oleracea       16832
Brassica_rapa   16715
Thellungiella_parvula   16641
Eutrema_salsugineum     16946
Arabis_alpina   12692
Boechera_stricta        17514
Camelina_sativa 17938
Capsella_rubella        17691
Arabidopsis_lyrata      17813
Arabidopsis_halleri     17390
Aethionema_arabicum     14359

# Numero de genes codificantes conservados por posicion y ortologia
cat Position_and_Orthologos_At_summary_results.tsv
Brassica_oleracea       19113
Brassica_rapa   19107
Thellungiella_parvula   19273
Eutrema_salsugineum     19552
Arabis_alpina   14092
Boechera_stricta        21026
Camelina_sativa 21752
Capsella_rubella        21316
Arabidopsis_lyrata      22197
Arabidopsis_halleri     21125
Aethionema_arabicum     15893

```
