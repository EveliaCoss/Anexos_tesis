# 03_Analisis con MEME suite (motivos)

## Analisis de microsecuencia de todos los transcritos

#### Eliminar los archivos vacios de MEME_suite

Los archivos meme.txt de cada carpeta (lncRNA) fueron analizados y se movieron a la carpeta (*./MEME_suite/MEME_filter/*)

```
# detectar archivos txt vacios y moverlos, lo malo es que como tienen el mismo nombre se sobreesciben y las carpetas quedan vacias
find ./MEME_suite/ -empty -name "*.txt" -exec mv {} ./MEME_filter/ \;
# mover las carpetas vacias
find ./MEME_suite/ -empty -type d -exec mv {} ./MEME_filter/ \;

# ejemplo = https://stackoverflow.com/questions/23254002/how-can-i-check-if-a-file-is-empty
# https://stackoverflow.com/questions/9964823/how-to-check-if-a-file-is-empty-in-bash/9964890
```

### Cargar datos en R

Cargar JAVA en R.

```
# install.packages("rJava")
# Despues de descargar JAVA e instalarlo en esta ubicacion lo pude llamar.
# Sys.setenv(JAVA_HOME="C:/Program Files/Java/jre1.8.0_311/")
library(rJava)

if (!require("pacman")) install.packages("pacman")
pacman::p_load_gh(
    "trinker/qdapDictionaries",
    "trinker/qdapRegex",
    "trinker/qdapTools",
    "trinker/qdap"
)


# Cargar los siguientes paquetes
library(memes)
library(magrittr)
library(GenomicRanges)
library(rJava)
library(qdap)
library(dplyr)
```

### Tabla del numero de motivos presentes en cada lncRNA por especie (v2)

Nueva matriz con la seleccion de Evalue 1-e3 y que al menos el motivo se encuentre en dos especies para hablar de microsecuencia.

```
library(memes)
library(magrittr)
library(GenomicRanges)
library(rJava)
library(qdap)
library(dplyr)

# Declarar las variables de entrada y salida
lncRNA_number_of_motif_db <-NULL
#lncRNA_number_of_motif_matrix <-NULL

lncRNA_number_of_motif_matrix <- data.frame(rep(0,12))
row.names(lncRNA_number_of_motif_matrix) <- c("Aethionema_arabicum", "Arabidopsis_halleri", "Arabidopsis_lyrata" , "Arabidopsis_thaliana",  "Arabis_alpina", "Boechera_stricta", "Brassica_oleracea", "Brassica_rapa", "Camelina_sativa", "Capsella_rubella", "Eutrema_parvulum", "Eutrema_salsugineum")

species <- c("Aethionema_arabicum", "Arabidopsis_halleri", "Arabidopsis_lyrata" , "Arabidopsis_thaliana",  "Arabis_alpina", "Boechera_stricta", "Brassica_oleracea", "Brassica_rapa", "Camelina_sativa", "Capsella_rubella", "Eutrema_parvulum", "Eutrema_salsugineum")

# Declarar variables para generar la matriz de entrada
# # Matriz de los motivos reportados en el lncRNA
lncRNA_motif_db <- NULL
lncRNA_motif_selected <- NULL
lncRNA_motif_matrix <-NULL

lncRNA_matrix_db <-data.frame(Motif = "",
                             species = "")

# Tabla del numero de motivos presentes en cada especie en cada lncRNA

for (i in 1:length(lincRNAs_names)) {
#for (i in 1:5) {
  # Cargar los resultados obtenidos con MEME para cada lncRNA.
  lncRNA_meme_path <-paste(complete_dir[i],"meme.txt", sep="/")
  lncRNA_meme_txt <-importMeme(lncRNA_meme_path, parse_genomic_coord = FALSE, combined_sites = FALSE) 
  
  # Matriz de los motivos reportados en el lncRNA
  # Generar Matriz de las especies en donde se encuentra cada motivo
  # Con seleccion de p value <=1e-3
  for (j in 1:5) { 
  lncRNA_motif_db <- lncRNA_meme_txt$sites_hits[[j]]
  lncRNA_motif_selected <- lncRNA_motif_db %>% filter(pvalue <= 1e-3) %>% select(sequence)
  
  # Eliminar motivos donde solo hay una especie, por cada lncRNA
    if (length(lncRNA_motif_selected$sequence) >= 2){
      lncRNA_motif_matrix <-data.frame(Motif = c(rep(paste("Motif", j, sep = ""), length(lncRNA_motif_selected$sequence))), 
                             species = c(lncRNA_motif_selected$sequence))
      lncRNA_matrix_db <- rbind(lncRNA_matrix_db, lncRNA_motif_matrix)
    }
  }

  # Eliminar motivos donde solo hay una especie, por cada lncRNA
      if (length(lncRNA_matrix_db$species) >= 2){
    
        # Eliminar la primera final (helper row)
        lncRNA_matrix_db <- lncRNA_matrix_db[-1,]
        # Tabla de acumulacion del numero de motivos presentes en cada especie en cada lncRNA
      lncRNA_number_of_motif_db <- wfm(lncRNA_matrix_db[, 1], lncRNA_matrix_db[, 2], lower.case = FALSE)
      lncRNA_number_of_motif_db <-as.data.frame(lncRNA_number_of_motif_db)
      colnames(lncRNA_number_of_motif_db)<- lincRNAs_names[i]
    
      # combinar archivos dependiendo de la fila (row)
      name <-lincRNAs_names[i]
      lncRNA_number_of_motif_matrix[name] <- ifelse(rownames(lncRNA_number_of_motif_matrix) %in% rownames(lncRNA_number_of_motif_db), lncRNA_number_of_motif_db[,1], 0)
      
      # Reiniciar variable
      lncRNA_matrix_db <-data.frame(Motif = "",
                             species = "")
    
      }
}

lncRNA_number_of_motif_matrix <- lncRNA_number_of_motif_matrix[,-1]
#lncRNA_matrix_db
head(lncRNA_number_of_motif_matrix)
#lncRNA_number_of_motif_db

# Respaldar el archivo
lncRNA_number_of_motif_matrix_evalue <- lncRNA_number_of_motif_matrix
head(lncRNA_number_of_motif_matrix_evalue)
```

Extraer el nombre de los AtlncRNas conservados

```
AtlincRNAs_conserved_by_synteny_seq_names <- colnames(lncRNA_number_of_motif_matrix_evalue)
head(AtlincRNAs_conserved_by_synteny_seq_names)
```

### Tabla del numero de motivos presentes en cada lncRNA por especie (v2)

Nueva matriz con la seleccion de Evalue 1-e3 y que al menos el motivo se encuentre en dos especies para hablar de microsecuencia.

```
library(memes)
library(magrittr)
library(GenomicRanges)
library(rJava)
library(qdap)
library(dplyr)

# Declarar las variables de entrada y salida
lncRNA_number_of_motif_db <-NULL
#lncRNA_number_of_motif_matrix <-NULL

lncRNA_number_of_motif_matrix <- data.frame(rep(0,12))
row.names(lncRNA_number_of_motif_matrix) <- c("Aethionema_arabicum", "Arabidopsis_halleri", "Arabidopsis_lyrata" , "Arabidopsis_thaliana",  "Arabis_alpina", "Boechera_stricta", "Brassica_oleracea", "Brassica_rapa", "Camelina_sativa", "Capsella_rubella", "Eutrema_parvulum", "Eutrema_salsugineum")

species <- c("Aethionema_arabicum", "Arabidopsis_halleri", "Arabidopsis_lyrata" , "Arabidopsis_thaliana",  "Arabis_alpina", "Boechera_stricta", "Brassica_oleracea", "Brassica_rapa", "Camelina_sativa", "Capsella_rubella", "Eutrema_parvulum", "Eutrema_salsugineum")

# Declarar variables para generar la matriz de entrada
# # Matriz de los motivos reportados en el lncRNA
lncRNA_motif_db <- NULL
lncRNA_motif_selected <- NULL
lncRNA_motif_matrix <-NULL

lncRNA_matrix_db <-data.frame(Motif = "",
                             species = "")

# Tabla del numero de motivos presentes en cada especie en cada lncRNA

for (i in 1:length(lincRNAs_names)) {
#for (i in 1:5) {
  # Cargar los resultados obtenidos con MEME para cada lncRNA.
  lncRNA_meme_path <-paste(complete_dir[i],"meme.txt", sep="/")
  lncRNA_meme_txt <-importMeme(lncRNA_meme_path, parse_genomic_coord = FALSE, combined_sites = FALSE) 
  
  # Matriz de los motivos reportados en el lncRNA
  # Generar Matriz de las especies en donde se encuentra cada motivo
  # Con seleccion de p value <=1e-3
  for (j in 1:5) { 
  lncRNA_motif_db <- lncRNA_meme_txt$sites_hits[[j]]
  lncRNA_motif_selected <- lncRNA_motif_db %>% filter(pvalue <= 1e-3) %>% select(sequence)
  
  # Eliminar motivos donde solo hay una especie, por cada lncRNA
    if (length(lncRNA_motif_selected$sequence) >= 2){
      lncRNA_motif_matrix <-data.frame(Motif = c(rep(paste("Motif", j, sep = ""), length(lncRNA_motif_selected$sequence))), 
                             species = c(lncRNA_motif_selected$sequence))
      lncRNA_matrix_db <- rbind(lncRNA_matrix_db, lncRNA_motif_matrix)
    }
  }

  # Eliminar motivos donde solo hay una especie, por cada lncRNA
      if (length(lncRNA_matrix_db$species) >= 2){
    
        # Eliminar la primera final (helper row)
        lncRNA_matrix_db <- lncRNA_matrix_db[-1,]
        # Tabla de acumulacion del numero de motivos presentes en cada especie en cada lncRNA
      lncRNA_number_of_motif_db <- wfm(lncRNA_matrix_db[, 1], lncRNA_matrix_db[, 2], lower.case = FALSE)
      lncRNA_number_of_motif_db <-as.data.frame(lncRNA_number_of_motif_db)
      colnames(lncRNA_number_of_motif_db)<- lincRNAs_names[i]
    
      # combinar archivos dependiendo de la fila (row)
      name <-lincRNAs_names[i]
      lncRNA_number_of_motif_matrix[name] <- ifelse(rownames(lncRNA_number_of_motif_matrix) %in% rownames(lncRNA_number_of_motif_db), lncRNA_number_of_motif_db[,1], 0)
      
      # Reiniciar variable
      lncRNA_matrix_db <-data.frame(Motif = "",
                             species = "")
    
      }
}

lncRNA_number_of_motif_matrix <- lncRNA_number_of_motif_matrix[,-1]
#lncRNA_matrix_db
head(lncRNA_number_of_motif_matrix)
#lncRNA_number_of_motif_db

# Respaldar el archivo
lncRNA_number_of_motif_matrix_evalue <- lncRNA_number_of_motif_matrix
head(lncRNA_number_of_motif_matrix_evalue)
```

Extraer el nombre de los AtlncRNas conservados

```
AtlincRNAs_conserved_by_synteny_seq_names <- colnames(lncRNA_number_of_motif_matrix_evalue)
head(AtlincRNAs_conserved_by_synteny_seq_names)
```

Lista de los AtlincRNAs conservados por sintenia y microsecuencia

```
write.table(AtlincRNAs_conserved_by_synteny_seq_names, file = "S:/Repositorio_scripts/Sequence_analysis/Results_R_MEME/AtlincRNAs_conserved_by_synteny_seq_names.txt", row.names = F) #(11/Julio/2022)
