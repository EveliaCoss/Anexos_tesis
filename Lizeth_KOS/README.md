# Visualizacion de genes conservados por sintenia en R

```
library(dplyr)
all_AtlincRNA_CSS_IDs %>% filter(AtlincRNAs %in% c("AT2G04405", "AT3G61198", "AT4G00975", "AT5G01775"))
```

Los lncRNAs AT2G04405, AT3G61198 y AT5G01775 se conservan por sintenia y microhomologia a traves de las especies.

Cargar los archivos:

```
# ----Cargar archivos a R----
# Arabidopsis thaliana coordendadas
#getwd()
AtlincRNA_AT2G044051_coord <- read.table("D:/Repositorio_scripts/Synteny_analysis/At_Sinteny_block_genes/At_lncRNA_AT2G04405_syntenic_block_GENES_sorted.bed", header = FALSE)
colnames(AtlincRNA_AT2G044051_coord) <- c("molecule", "gene", "start", "end")

# Arabidopsis lyrata
#Al_chr3_IPS1_lncRNA_coord <- read.table("D:/Mazorka/results/Phylostratographic/Sinteny_block_genes/At_lncRNA_AT3G09922_Arabidopsis_lyrata_syntenic_block_GENES_3_rename.bed",header = FALSE)
AllincRNA_AT2G04405_coord <- read.table("D:/Repositorio_scripts/Synteny_analysis/At_Sinteny_block_genes/At_lncRNA_AT2G04405_Arabidopsis_lyrata_syntenic_block_GENES_5_sorted.bed",header = FALSE)
AllincRNA_AT2G04405_coord$V1 <- "Arabidopsis_lyrata_chr5"
colnames(AllincRNA_AT2G04405_coord) <- c("molecule", "gene", "start", "end")

# Arabidopsis halleri
AhlincRNA_AT2G04405_coord <- read.table("D:/Repositorio_scripts/Synteny_analysis/At_Sinteny_block_genes/At_lncRNA_AT2G04405_Arabidopsis_halleri_syntenic_block_GENES_Scaffold5103_rename.bed",header = FALSE)
colnames(AhlincRNA_AT2G04405_coord) <- c("molecule", "gene", "start", "end")

# Camelina sativa
CslincRNA_AT2G04405_coord <- read.table("D:/Repositorio_scripts/Synteny_analysis/At_Sinteny_block_genes/At_lncRNA_AT2G04405_Camelina_sativa_syntenic_block_GENES_CM002730.1_rename.bed",header = FALSE)
colnames(CslincRNA_AT2G04405_coord) <- c("molecule", "gene", "start", "end")

# tiene tres duplicados
#At_lncRNA_AT3G09922_Camelina_sativa_syntenic_block_GENES_CM002729.1_rename.bed
#At_lncRNA_AT3G09922_Camelina_sativa_syntenic_block_GENES_CM002743.1_rename.bed

#At_lncRNA_AT3G09922_Camelina_sativa_syntenic_block_GENES_CM002747.1_rename.bed

# capsella rubella
CrlincRNA_AT2G04405_coord <- read.table("D:/Repositorio_scripts/Synteny_analysis/At_Sinteny_block_genes/At_lncRNA_AT2G04405_Capsella_rubella_syntenic_block_GENES_scaffold_5_rename.bed",header = FALSE)
colnames(CrlincRNA_AT2G04405_coord) <- c("molecule", "gene", "start", "end")

# Eutrema_salsugineum
EslincRNA_AT2G04405_coord <- read.table("D:/Repositorio_scripts/Synteny_analysis/At_Sinteny_block_genes/At_lncRNA_AT2G04405_Eutrema_salsugineum_syntenic_block_GENES_NW_006256547.1_sorted.bed",header = FALSE)
colnames(EslincRNA_AT2G04405_coord) <- c("molecule", "gene", "start", "end")

#At_lncRNA_AT3G09922_Eutrema_salsugineum_syntenic_block_GENES_NW_006256829.1_rename.bed
#At_lncRNA_AT3G09922_Eutrema_salsugineum_syntenic_block_GENES_NW_006256885.1_rename.bed

# Thellungiella parvula
TplincRNA_AT2G04405_coord <- read.table("D:/Repositorio_scripts/Synteny_analysis/At_Sinteny_block_genes/At_lncRNA_AT2G04405_Thellungiella_parvula_syntenic_block_GENES_CM001188.1_rename.bed",header = FALSE)
colnames(TplincRNA_AT2G04405_coord) <- c("molecule", "gene", "start", "end")

#At_lncRNA_AT3G09922_Eutrema_parvulum_syntenic_block_GENES_CM001189.1_rename.bed
#At_lncRNA_AT3G09922_Eutrema_parvulum_syntenic_block_GENES_CM001192.1_rename.bed

# Brassica rapa coordendadas
# A05
BrlincRNA_AT2G04405_coord <- read.table("D:/Repositorio_scripts/Synteny_analysis/At_Sinteny_block_genes/At_lncRNA_AT2G04405_Brassica_rapa_syntenic_block_GENES_A09_rename.bed",header = FALSE)
colnames(BrlincRNA_AT2G04405_coord) <- c("molecule", "gene", "start", "end")

#At_lncRNA_AT3G09922_Brassica_rapa_syntenic_block_GENES_A01_rename.bed
#At_lncRNA_AT3G09922_Brassica_rapa_syntenic_block_GENES_A02_rename.bed
#At_lncRNA_AT3G09922_Brassica_rapa_syntenic_block_GENES_A03_rename.bed
#At_lncRNA_AT3G09922_Brassica_rapa_syntenic_block_GENES_A05_rename.bed
#At_lncRNA_AT3G09922_Brassica_rapa_syntenic_block_GENES_A10_rename.bed

# Brassica oleracea
# C5
BolincRNA_AT2G04405_coord <- read.table("D:/Repositorio_scripts/Synteny_analysis/At_Sinteny_block_genes/At_lncRNA_AT2G04405_Brassica_oleracea_syntenic_block_GENES_C9_rename.bed",header = FALSE)
colnames(BolincRNA_AT2G04405_coord) <- c("molecule", "gene", "start", "end")

#At_lncRNA_AT3G09922_Brassica_oleracea_syntenic_block_GENES_C1_rename.bed
#At_lncRNA_AT3G09922_Brassica_oleracea_syntenic_block_GENES_C2_rename.bed
#At_lncRNA_AT3G09922_Brassica_oleracea_syntenic_block_GENES_C3_rename.bed
#At_lncRNA_AT3G09922_Brassica_oleracea_syntenic_block_GENES_C5_rename.bed
#At_lncRNA_AT3G09922_Brassica_oleracea_syntenic_block_GENES_C7_rename.bed
#At_lncRNA_AT3G09922_Brassica_oleracea_syntenic_block_GENES_C9_rename.bed

#Aethionema_arabicum
AealincRNA_AT2G04405_IPS1coord <- read.table("D:/Repositorio_scripts/Synteny_analysis/At_Sinteny_block_genes/At_lncRNA_AT2G04405_Aethionema_arabicum_syntenic_block_GENES_KE151134.1_rename.bed",header = FALSE)
colnames(AealincRNA_AT2G04405_IPS1coord) <- c("molecule", "gene", "start", "end")

#At_lncRNA_AT3G09922_Aethionema_arabicum_syntenic_block_GENES_KE151514.1_rename.bed
#At_lncRNA_AT3G09922_Aethionema_arabicum_syntenic_block_GENES_KE151550.1_rename.bed

#Boechera_stricta
BslincRN_AT2G04405_coord <- read.table("D:/Repositorio_scripts/Synteny_analysis/At_Sinteny_block_genes/At_lncRNA_AT2G04405_Boechera_stricta_syntenic_block_GENES_MLHT01000028.1_rename.bed",header = FALSE)
colnames(BslincRN_AT2G04405_coord) <- c("molecule", "gene", "start", "end")

#At_lncRNA_AT3G09922_Boechera_stricta_syntenic_block_GENES_MLHT01001405.1_rename.bed
```

Ubicar la direccion de las secuencias:

```
# genes codificantes
At_THESIS_pcoding_annotated_bed %>% filter(GeneID %in% c( "AT2G02870", "AT2G02880", "AT2G02890", "AT2G02910", "AT2G02930", "AT2G02950", "AT2G02955"))

# lncRNAs
At_THESIS_lncRNAs_annotated_bed %>% filter(GeneID == "AT2G04405")
```

Anexar informacion y modificar dataframe

```
# ----Agregar columnas de orientacion de los genes----
library(dplyr)

# Arabidopsis thaliana
# lncRNA AT2G04405
# strand AT2G02870, AT2G02880, AT2G02890, AT2G02910, AT2G02930, AT2G02950, AT2G02955, AT2G04405
# orientation forward = 1  and rever = 1
strand <- data.frame(c("forward", "forward","reverse", "reverse", "reverse", "reverse","forward", "forward")) # At, Al y cr
orientation <- data.frame(c(0,0,1,1,1,1,0,0))

#data
geneInfo<- data.frame(
  gene = c("AT2G02870", "AT2G02880", "AT2G02890", "AT2G02910", "AT2G02930", "AT2G02950", "AT2G02955", "AT2G04405"),
  strand = c("forward", "forward","reverse", "reverse", "reverse", "reverse","forward", "forward"),
  orientation = c(0,0,1,1,1,1,0,0)
)

# A. thaliana
AtlincRNA_AT2G044051_coord <- left_join(AtlincRNA_AT2G044051_coord, geneInfo, by="gene")

# Arabidopsis lyrata
AllincRNA_AT2G04405_coord <- left_join(AllincRNA_AT2G04405_coord, geneInfo, by="gene")

# Arabidopsis halleri
# AT3G09922 AT3G09925
AhlincRNA_AT2G04405_coord <- left_join(AhlincRNA_AT2G04405_coord, geneInfo, by="gene")

# Camelina sativa
CslincRNA_AT2G04405_coord <- left_join(CslincRNA_AT2G04405_coord, geneInfo, by="gene")

# capsella rubella
CrlincRNA_AT2G04405_coord <- left_join(CrlincRNA_AT2G04405_coord, geneInfo, by="gene")

# Thellungiella parvula
TplincRNA_AT2G04405_coord <- left_join(TplincRNA_AT2G04405_coord, geneInfo, by="gene")
 
# Eutrema_salsugineum
EslincRNA_AT2G04405_coord <- left_join(EslincRNA_AT2G04405_coord, geneInfo, by="gene")

# Brassica rapa
BrlincRNA_AT2G04405_coord <- left_join(BrlincRNA_AT2G04405_coord, geneInfo, by="gene")

# Brassica oleracea
BolincRNA_AT2G04405_coord <- left_join(BolincRNA_AT2G04405_coord, geneInfo, by="gene")

# Aethionema_arabicum
AealincRNA_AT2G04405_IPS1coord <- left_join(AealincRNA_AT2G04405_IPS1coord, geneInfo, by="gene")

# Boechera_stricta
BslincRN_AT2G04405_coord <-  left_join(BslincRN_AT2G04405_coord, geneInfo, by="gene")

# chequeo
AtlincRNA_AT2G044051_coord # A. thaliana
AllincRNA_AT2G04405_coord # Arabidopsis lyrata
AhlincRNA_AT2G04405_coord # Arabidopsis halleri
CslincRNA_AT2G04405_coord # Camelina sativa
CrlincRNA_AT2G04405_coord # capsella rubella
EslincRNA_AT2G04405_coord # Eutrema_salsugineum
TplincRNA_AT2G04405_coord # Thellungiella parvula
BrlincRNA_AT2G04405_coord # Brassica rapa
BolincRNA_AT2G04405_coord # Brassica oleracea
AealincRNA_AT2G04405_IPS1coord # Aethionema_arabicum
BslincRN_AT2G04405_coord # Boechera_stricta

# unir coordenadas (Archivo corto)
lincRNA_AT2G04405_syntenic_block <- rbind(AtlincRNA_AT2G044051_coord, AllincRNA_AT2G04405_coord, AhlincRNA_AT2G04405_coord, CslincRNA_AT2G04405_coord, CrlincRNA_AT2G04405_coord,EslincRNA_AT2G04405_coord, TplincRNA_AT2G04405_coord,BrlincRNA_AT2G04405_coord, BolincRNA_AT2G04405_coord, AealincRNA_AT2G04405_IPS1coord, BslincRN_AT2G04405_coord  )


# Renombrar molecule para reducir IDs
lincRNA_AT2G04405_syntenic_blockmodified <- mutate(lincRNA_AT2G04405_syntenic_block, molecule = ifelse(molecule == "Arabidopsis_thaliana_chr2", "Arabidopsis_thaliana", ifelse(molecule ==  "Arabidopsis_lyrata_chr5", "Arabidopsis_lyrata", ifelse(molecule == "Arabidopsis_halleri_Scaffold5103", "Arabidopsis_halleri", ifelse(molecule == "Camelina_sativa_CM002730.1", "Camelina_sativa", ifelse(molecule == "Capsella_rubella_scaffold_5", "Capsella_rubella", ifelse(molecule == "NW_006256547.1", "Eutrema_salsugineum", ifelse(molecule == "Eutrema_salsugineum_NW_006256885.1", "Eutrema_salsugineum", ifelse(molecule == "Eutrema_salsugineum_NW_006256885.1", "Eutrema_salsugineum", ifelse(molecule == "Thellungiella_parvula_CM001188.1", "Thellungiella_parvula", ifelse(molecule == "Brassica_oleracea_C9", "Brassica_oleracea",ifelse(molecule == "Aethionema_arabicum_KE151134.1", "Aethionema_arabicum",ifelse(molecule == "Brassica_rapa_A09", "Brassica_rapa", "Boechera_stricta")))))))))))))
                                   
head(lincRNA_AT2G04405_syntenic_blockmodified)
```

Si se quiere reducir tamanos de las secuencias

```
#Reducir el sitio de inicio del gen codificante AT3G09920 para visualizar los otros genes
lincRNA_AT2G04405_syntenic_blockmodified[38,3] <- 2774231 # valor original start 2323033 # Thellungiella parvula
lincRNA_AT2G04405_syntenic_blockmodified[33,4] <- 6352544 # valor original end 6852836 # Eutrema_salsugineum
lincRNA_AT2G04405_syntenic_blockmodified[19,3] <- 3827125 # valor original start AT3G09940	3819691 # camelina sativa
lincRNA_AT2G04405_syntenic_blockmodified[17,4] <- 3816449 # Valor original end AT3G09925	3815390	3823816
lincRNA_AT2G04405_syntenic_blockmodified[18,4] <- 3824523 # valor original end AT3G09930	3817110	3826513
```

Graficar 

```
library(gggenes)
library(ggplot2)
library(cowplot)

# Grafica
# distribucion de los genes en IPS1
PaperVersion_AT2G04405_allspecies_gggenes_plot <- lincRNA_AT2G04405_syntenic_blockmodified %>% mutate(molecule = factor(molecule, levels = c("Arabidopsis_thaliana", "Arabidopsis_lyrata", "Arabidopsis_halleri","Camelina_sativa",  "Capsella_rubella",  "Boechera_stricta", "Eutrema_salsugineum", "Thellungiella_parvula", "Brassica_rapa",  "Brassica_oleracea", "Aethionema_arabicum"))) %>% ggplot(aes(xmin = start, xmax = end, y = molecule, fill = gene, forward = orientation)) +
  geom_gene_arrow() +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  #scale_fill_manual(values=cbPalette2, name ="genes") +
  scale_fill_brewer(palette = "Set3") +
  theme_genes() #+
# Si quieres poner los genes arriba y eliminar los nombres de las especies
  #theme(legend.position="top", legend.text=element_text(size=12, family="sans"), text=element_text(size=12, family="sans"),  # Arial
  #      axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y = element_text(), strip.text.y = element_blank(), legend.title = element_blank()) # remove information in X axis

PaperVersion_AT2G04405_allspecies_gggenes_plot

# Guardar info
save_plot("AT2G04405_CSS_plot.pdf", PaperVersion_AT2G04405_allspecies_gggenes_plot, base_height = 8, base_width = 12, base_asp = 1.1)

# The geom_motif() is defined in ggtree and it is a wrapper layer of the gggenes::geom_gene_arrow().
#ggsave(file="../Conservacion_x_secuencia/Paper_conservation/Figura2C_lincRNAIPS1SINTENY_block.tiff", plot = PaperVersion_lincRNAIPS1_allspecies_gggenes_plot, device= "tiff",width = 20, height = 10)
```






