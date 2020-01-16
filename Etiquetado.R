############################################################################################################################################################
############################################### ETIQUETAR MATRIZ DEL DEEPMAP  19Q2 ###############################################################################
############################################################################################################################################################
setwd("/home/estrella/TFM/Datos/depmap")
#save.image("./20190725_ETIQUETADO.RData")
lista_genes_prot <- read.table ("/home/estrella/TFM/Datos/depmap/mart_export.txt", sep = "\t", header = T)
lista_genes_prot <-as.data.frame(lista_genes_prot)
library(dplyr)
library(reshape2)
library(data.table)
lista_genes_prot2 <- lista_genes_prot%>%filter(lista_genes_prot$Gene.type=="protein_coding") #Tiene 22722 genes que codifican para proteina
depmap<-fread("/home/estrella/TFM/Datos/depmap/CCLE_expression_19Q2csv")
depmap2<- as.data.frame(depmap)
#colnames(depmap2)[grepl("ENSG00000204805",colnames(depmap))] #Para comprobar si estaban

genes_depmap<-data.frame(colnames(depmap2))
colnames(genes_depmap)<-"genes"
div<-colsplit(genes_depmap$genes, "\\s", c("name","ID")) #Dividimos los colnames en ID de Ensembl y name (HUGO), para luego cruzas los ID con la lista del genes_prot.
t1<-data.frame(gsub("\\(", "", div$ID))
t2<-data.frame(gsub("\\)", "", t1$gsub............div.ID.))
colnames(t2) <-"ID" #Los ID de ensembl sin parentesis
div2 <-data.frame(cbind(div[,1], t2))
GenesFiltrados<-as.data.frame(div2[div2$ID%in%lista_genes_prot2$Gene.stable.ID,]) #LISTADO DE GENES QUE CODIFICAN PARA PROTEINA.

nombres_columnas<-c(gsub("\\s\\(", "_", colnames(depmap2))) #Modifico el nombre de las columnas para que tengan nombreHUGO_IDemsembl y poder cruzarlo con la lista.
nombres_columnas2<-c(gsub("\\)", "", nombres_columnas))
colnames(depmap2)<-nombres_columnas2 

GenesFiltrados2<- data.frame(paste(GenesFiltrados$div...1., GenesFiltrados$ID, sep="_"))#Reconvierto el formato para que sea igual al de las columnas del la matriz depmap.
colnames(GenesFiltrados2) <- "gen_ID"

MFiltradaProt<-data.frame(cbind(depmap2$V1, depmap2[,colnames(depmap2) %in% as.character(GenesFiltrados2$gen_ID)])) ####### MATRIZ CON LAS COLUMNAS DE LOS GENES PROTEIN CODING.

library(data.table)
lista_lineas_tej <- fread ("/home/estrella/TFM/Datos/depmap/sample_info_19Q2.csv")
lista_lineas_tej2 <-as.data.frame(lista_lineas_tej[,c(1,3)])
colnames(lista_lineas_tej2)<-c("depmap2.V1", "CCLE_ID")
#save.image("./20190725_Etiquetado.RData")

M19Q2<-merge(lista_lineas_tej2, MFiltradaProt, by="depmap2.V1")
colnames(M19Q2)<-stringr::str_extract(colnames(M19Q2), "^[^_]+(?=_)") ######## MATRIZ ETIQUETADA 19Q2 (nombres de los tej y los Id) y FILTRADA (solo genes que codifican par proteinas con nomenclaruta HUGO)
M19Q2<- M19Q2[-c(1)]
write.table (M19Q2, file = "M19Q2_etiquetada.csv", col.names = TRUE, sep="\t") #El archivo que sale esta mal y no lo puedo abrir después.

######## Filtro de la matriz 19Q2 con los genes del modelo metabólico #####
GenesMetabolismo <- read.table ("./Recon3D_genes.lst")
colnames(GenesMetabolismo)<-'GenesMetabol'
M19Q2_metabol <-cbind(M19Q2[c(1:2)],M19Q2[,colnames(M19Q2) %in% as.character(GenesMetabolismo$GenesMetabol)])
write.table (M19Q2_metabol, file = "M19Q2_metabol.csv", col.names = TRUE, sep="\t")

#####################################################################################################################################################
########################################### ETIQUETAR MATRIZ CERES 19Q2 #############################################################################
#####################################################################################################################################################

# lista_genes_prot <- read.table ("/home/estrella/TFM/Datos/depmap/mart_export.txt", sep = "\t", header = T)
# lista_genes_prot <-as.data.frame(lista_genes_prot)
library(dplyr)
library(reshape2)
library(data.table)
#lista_genes_prot2 <- lista_genes_prot%>%filter(lista_genes_prot$Gene.type=="protein_coding") #Tiene 22722 genes que codifican para proteina
ceres<-fread("/home/estrella/TFM/Datos/depmap/Achilles_gene_effect.csv") #El nombre de los genes es el HUGO ID, con un numero entre parentesis.
ceres2<- as.data.frame(ceres)
colnames(ceres2)<- sub("\\s.*", "", colnames(ceres2))

# genes_ceres<-data.frame(colnames(ceres2))
# colnames(genes_ceres)<-"genes"
# div<-colsplit(genes_ceres$genes, "\\s", c("name","ID")) #Dividimos los colnames en ID de Ensembl y name (HUGO), para luego cruzas los ID con la lista del genes_prot.
# t1<-data.frame(gsub("\\(", "", div$ID))
# t2<-data.frame(gsub("\\)", "", t1$gsub............div.ID.))
# colnames(t2) <-"ID" #Los ID de ensembl sin parentesis
# div2 <-data.frame(cbind(div[,1], t2))
GenesFiltradosCeres<-as.data.frame(div2[div2$ID%in%lista_genes_prot2$Gene.stable.ID,]) #LISTADO DE GENES QUE CODIFICAN PARA PROTEINA.
lista_prot_HUGO<- as.data.frame(colnames(M19Q2)) #Es el mismo listado pero en nomenclatura HUGO.

# nombres_columnas<-c(gsub("\\s\\(", "_", colnames(ceres2))) #Modifico el nombre de las columnas para que tengan nombreHUGO_IDemsembl y poder cruzarlo con la lista.
# nombres_columnas2<-c(gsub("\\)", "", nombres_columnas))
# colnames(ceres2)<-nombres_columnas2 
# GenesFiltradosCeres2<- data.frame(paste(GenesFiltradosCeres$div...1., GenesFiltradosCeres$ID, sep="_"))#Reconvierto el formato para que sea igual al de las columnas del la matriz depmap.
# colnames(GenesFiltradosCeres2) <- "gen_ID"

MFiltradaProtCeres<-data.frame(cbind(ceres2$V1, ceres2[,colnames(ceres2) %in% as.character(lista_prot_HUGO$`colnames(M19Q2)`)])) ####### MATRIZ CON LAS COLUMNAS DE LOS GENES PROTEIN CODING.

lista_lineas_tej <- fread ("/home/estrella/TFM/Datos/depmap/sample_info_19Q2.csv", sep = ",", header = T)
lista_lineas_tej2 <-as.data.frame(lista_lineas_tej[,c(1,3)])
colnames(lista_lineas_tej2)<-c("ceres2.V1", "CCLE_ID")
#save.image("./20190725_genes.RData")

M19Q2_ceres<-merge(lista_lineas_tej2, MFiltradaProtCeres, by="ceres2.V1")
M19Q2_ceres<- M19Q2_ceres[-c(1)]
write.table (M19Q2_ceres, file = "M19Q2_ceres.csv", col.names = TRUE, sep="\t")


######## Filtro de la matriz 19Q2 con los genes del modelo metabólico #####
GenesMetabolismo <- read.table ("./Recon3D_genes.lst")
colnames(GenesMetabolismo)<-'GenesMetabol'
M19Q2_ceres_metabol <-cbind(M19Q2_ceres[c(1)], M19Q2_ceres[,colnames(M19Q2_ceres) %in% as.character(GenesMetabolismo$GenesMetabol)])
write.table (M19Q2_ceres_metabol, file = "M19Q2_ceres_metabol.csv", col.names = TRUE, sep="\t")

############################################################################################################################################################
############################################### ETIQUETAR MATRIZ DEL DEEPMAP  19Q1 ###############################################################################
############################################################################################################################################################
setwd("/home/estrella/TFM/Datos/depmap/dep_Map_19Q1/")
#save.image("./20190725_ETIQUETADO.RData")
lista_genes_prot <- read.table ("/home/estrella/TFM/Datos/depmap/mart_export.txt", sep = "\t", header = T)
lista_genes_prot <-as.data.frame(lista_genes_prot)
library(dplyr)
library(reshape2)
library(data.table)
lista_genes_prot2 <- lista_genes_prot%>%filter(lista_genes_prot$Gene.type=="protein_coding") #Tiene 22722 genes que codifican para proteina
depmap<-fread("/home/estrella/TFM/Datos/depmap/dep_Map_19Q1/CCLE_depMap_19Q1_TPM.csv")
depmap2<- as.data.frame(depmap)
#colnames(depmap2)[grepl("ENSG00000204805",colnames(depmap))] #Para comprobar si estaban

genes_depmap<-data.frame(colnames(depmap2))
colnames(genes_depmap)<-"genes"
div<-colsplit(genes_depmap$genes, "\\s", c("name","ID")) #Dividimos los colnames en ID de Ensembl y name (HUGO), para luego cruzas los ID con la lista del genes_prot.
t1<-data.frame(gsub("\\(", "", div$ID))
t2<-data.frame(gsub("\\)", "", t1$gsub............div.ID.))
colnames(t2) <-"ID" #Los ID de ensembl sin parentesis
div2 <-data.frame(cbind(div[,1], t2))
GenesFiltrados<-as.data.frame(div2[div2$ID%in%lista_genes_prot2$Gene.stable.ID,]) #LISTADO DE GENES QUE CODIFICAN PARA PROTEINA.

nombres_columnas<-c(gsub("\\s\\(", "_", colnames(depmap2))) #Modifico el nombre de las columnas para que tengan nombreHUGO_IDemsembl y poder cruzarlo con la lista.
nombres_columnas2<-c(gsub("\\)", "", nombres_columnas))
colnames(depmap2)<-nombres_columnas2 

GenesFiltrados2<- data.frame(paste(GenesFiltrados$div...1., GenesFiltrados$ID, sep="_"))#Reconvierto el formato para que sea igual al de las columnas del la matriz depmap.
colnames(GenesFiltrados2) <- "gen_ID"

MFiltradaProt<-data.frame(cbind(depmap2$V1, depmap2[,colnames(depmap2) %in% as.character(GenesFiltrados2$gen_ID)])) ####### MATRIZ CON LAS COLUMNAS DE LOS GENES PROTEIN CODING.

lista_lineas_tej <- fread("/home/estrella/TFM/Datos/depmap/dep_Map_19Q1/sample_info.csv", sep = ",", header = T)
lista_lineas_tej2 <-as.data.frame(lista_lineas_tej[,c(2,7)])
colnames(lista_lineas_tej2)<-c("depmap2.V1", "CCLE_ID")
#save.image("./20190725_genes.RData")

M19Q1<-merge(lista_lineas_tej2, MFiltradaProt, by="depmap2.V1")
colnames(M19Q1)<-stringr::str_extract(colnames(M19Q1), "^[^_]+(?=_)") ######## MATRIZ ETIQUETADA 19Q2 (nombres de los tej y los Id) y FILTRADA (solo genes que codifican par proteinas con nomenclaruta HUGO)
M19Q1<- M19Q1[-c(1)]
write.table (M19Q1, file = "M19Q1.csv", col.names = TRUE, sep="\t")

######## Filtro de la matriz 19Q2 con los genes del modelo metabólico #####
GenesMetabolismo <- read.table ("../Recon3D_genes.lst")
colnames(GenesMetabolismo)<-'GenesMetabol'
M19Q1_metabol <-cbind(M19Q1[c(1:2)],M19Q1[,colnames(M19Q1) %in% as.character(GenesMetabolismo$GenesMetabol)])
write.table (M19Q1_metabol, file = "M19Q1_metabol.csv", col.names = TRUE, sep="\t")

#####################################################################################################################################################
########################################### ETIQUETAR MATRIZ CERES 19Q1 #############################################################################
#####################################################################################################################################################
# lista_genes_prot <- read.table ("/home/estrella/TFM/Datos/depmap/mart_export.txt", sep = "\t", header = T)
# lista_genes_prot <-as.data.frame(lista_genes_prot)
library(dplyr)
library(reshape2)
library(data.table)
lista_genes_prot2 <- lista_genes_prot%>%filter(lista_genes_prot$Gene.type=="protein_coding") #Tiene 22722 genes que codifican para proteina
ceres<-fread("/home/estrella/TFM/Datos/depmap/dep_Map_19Q1/gene_effect_corrected.csv") #El nombre de los genes es el HUGO ID, con un numero entre parentesis.
ceres2<- as.data.frame(ceres)
colnames(ceres2)<- sub("\\s.*", "", colnames(ceres2))

# genes_ceres<-data.frame(colnames(ceres2))
# colnames(genes_ceres)<-"genes"
# div<-colsplit(genes_ceres$genes, "\\s", c("name","ID")) #Dividimos los colnames en ID de Ensembl y name (HUGO), para luego cruzas los ID con la lista del genes_prot.
# t1<-data.frame(gsub("\\(", "", div$ID))
# t2<-data.frame(gsub("\\)", "", t1$gsub............div.ID.))
# colnames(t2) <-"ID" #Los ID de ensembl sin parentesis
# div2 <-data.frame(cbind(div[,1], t2))
GenesFiltradosCeres<-as.data.frame(div2[div2$ID%in%lista_genes_prot2$Gene.stable.ID,]) #LISTADO DE GENES QUE CODIFICAN PARA PROTEINA.
lista_prot_HUGO<- as.data.frame(colnames(M19Q1)) #Es el mismo listado pero en nomenclatura HUGO.

# nombres_columnas<-c(gsub("\\s\\(", "_", colnames(ceres2))) #Modifico el nombre de las columnas para que tengan nombreHUGO_IDemsembl y poder cruzarlo con la lista.
# nombres_columnas2<-c(gsub("\\)", "", nombres_columnas))
# colnames(ceres2)<-nombres_columnas2 
# GenesFiltradosCeres2<- data.frame(paste(GenesFiltradosCeres$div...1., GenesFiltradosCeres$ID, sep="_"))#Reconvierto el formato para que sea igual al de las columnas del la matriz depmap.
# colnames(GenesFiltradosCeres2) <- "gen_ID"

MFiltradaProtCeres<-data.frame(cbind(ceres2$V1, ceres2[,colnames(ceres2) %in% as.character(lista_prot_HUGO$`colnames(M19Q1)`)])) ####### MATRIZ CON LAS COLUMNAS DE LOS GENES PROTEIN CODING.

lista_lineas_tej <- fread("/home/estrella/TFM/Datos/depmap/dep_Map_19Q1/sample_info.csv", sep = ",", header = T)
lista_lineas_tej2 <-as.data.frame(lista_lineas_tej[,c(2,7)])
colnames(lista_lineas_tej2)<-c("ceres2.V1", "CCLE_ID")
#save.image("./20190725_genes.RData")

M19Q1_ceres<-merge(lista_lineas_tej2, MFiltradaProtCeres, by="ceres2.V1")
M19Q1_ceres<- M19Q1_ceres[-c(1)]
write.table (M19Q1_ceres, file = "M19Q1_ceres.csv", col.names = TRUE, sep="\t")


######## Filtro de la matriz 19Q2 con los genes del modelo metabólico #####
GenesMetabolismo <- read.table ("../Recon3D_genes.lst")
colnames(GenesMetabolismo)<-'GenesMetabol'
M19Q1_ceres_metabol <-cbind(M19Q1_ceres[c(1)], M19Q1_ceres[,colnames(M19Q1_ceres) %in% as.character(GenesMetabolismo$GenesMetabol)])
write.table (M19Q1_ceres_metabol, file = "M19Q1_ceres_metabol.csv", col.names = TRUE, sep="\t")

save.image("./20190725_Etiquetado.RData")
