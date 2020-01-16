getwd()
setwd("/home/estrella/TFM/Datos/depmap/FiltradoDatosYAnalisisDatosCrudos_(Agrupado por tejido)/")
# library(data.table)
# matrizDatosCrudos<-fread("../RELABELED_CCLE_depMap_19Q1_TPM.csv")
library(data.table)
matrizDatosCrudos<- read.table ('../M19Q2_etiquetada.csv', header = TRUE, sep= '\t')
#Tissue<- matrizDatosCrudos$CCLE
#matrizDatosCrudos<-cbind (Tissue, matrizDatosCrudos[, names(matrizDatosCrudos) != "tissue"])
library(dplyr)
library(reshape2)
div<-colsplit(matrizDatosCrudos$CCLE, "_", c("ID","TISSUE")) #Divides la columna cell_line(que tiene ID_tej) en dos columnas.
MatrizDatos <- cbind(div,matrizDatosCrudos) #Le unes lo que acabas de separar.

#------ AGRUPAR POR TEJIDO-----------------
unique(MatrizDatos$TISSUE)
tablaFreqTej <- data.frame(table(MatrizDatos$TISSUE))

v <- c(subset(tablaFreqTej, Freq > 5)) #Listado de los tejidos que con mas de 5 lineas por tejido.
v <- data.frame(subset(tablaFreqTej, Freq > 5))
MatrizFiltrada <-MatrizDatos[MatrizDatos$TISSUE%in%as.character(v$Var1),]
unique(MatrizFiltrada$TISSUE)

#--------MEDIA Y MEDIANA DE LOS GENES SEGUN TEJIDO-------------
########  MEDIA  ######
listado_tej<- as.character(v$Var1) #Listado de los tejidos que me valen
Matriz_media <- data.frame()
for (tejido in listado_tej) {
  print(tejido)
  mtejido<-MatrizFiltrada[which(MatrizFiltrada$TISSUE==tejido),][,4:length(MatrizFiltrada)]
  media_genes<-c()
  media_genes<-apply(mtejido,2,mean)
  media_genes<-t(as.data.frame(media_genes))
  rownames(media_genes)<-tejido
  Matriz_media <-rbind(Matriz_media, media_genes)
}

########  MEDIANA  ######

Matriz_mediana <- data.frame()
for (tejido in listado_tej) {
  print(tejido)
  mtejido<-MatrizFiltrada[which(MatrizFiltrada$TISSUE==tejido),][,4:length(MatrizFiltrada)]
  mediana_genes<-c()
  mediana_genes<-apply(mtejido,2,median)
  mediana_genes<-t(as.data.frame(mediana_genes))
  rownames(mediana_genes)<-tejido
  Matriz_mediana<-rbind(Matriz_mediana, mediana_genes)
}

#-------ELIMINANDO GENES NO EXPRESADOS (GENES QUE SE EXPRESEN EN MENOS DE UN 10% DE LAS LINEAS)-----------------------------
library(dplyr)
listado_genes<- colnames(MatrizFiltrada)[-(1:3)]
MatrizGenesNoExpresados<-data.frame() #Cuenta el numero de ceros. Numero de veces que no se expresa un gen.

for (tejido in listado_tej){
  print(tejido)
  t1 <-MatrizFiltrada%>%filter(TISSUE==tejido) #Agrupo por tejido
  vacio_tejido<-cbind(tejido, colSums(t1==0)  %>% .[-(1:3)] %>% as.data.frame() %>% t())
  MatrizGenesNoExpresados <-rbind(MatrizGenesNoExpresados,vacio_tejido)
}

colnames(tablaFreqTej) <- c("tejido", "Freq. Tej")
MatrizGenesNoExpresadosDef <- merge(tablaFreqTej, MatrizGenesNoExpresados, by="tejido") #Nos uniría la columna de frecuencias con el numero total de lineas por tejido.
rownames(MatrizGenesNoExpresadosDef)<-MatrizGenesNoExpresadosDef$tejido #Cambiamos el nombre a rownames la columna tejido. Para que solo se queden las filas y columnas con las que queremos operar.
MatrizGenesNoExpresadosDef<-MatrizGenesNoExpresadosDef[-c(1)]
p1<-apply(MatrizGenesNoExpresadosDef,2,as.numeric) %>% as.data.frame() #Al hacer esto te quita los row names por eso hacemos lo siguiente.
rownames(p1)<-rownames(MatrizGenesNoExpresadosDef)
frecuenciasNoExpr <- p1[-c(1)]/p1$`Freq. Tej` 
frecuenciasNoExpr[frecuenciasNoExpr > 0.1] <- NA #Cambio los valores de frecuencia que sean mayores a 0.1 por NA. Si obtengo una proporcion de ceros mayor del 0.9, quiere decir que solo se expresa en el 10% de las lineas y lo considero como no expresado.
MatrizGenesNoExpr <-frecuenciasNoExpr


#------------IDENTIFICAMOS GENES ALTAMENTE EXPRESADOS EN TODOS LOS TEJIDOS--------------
# Definimos los genes altamente expresados como aquellos que estan por encima del cuartil 75%.
#De la matriz que contiene las medias de cada gen por tejido (media_genes2). Hacer de cada columna (gen) el cuartil. Categorizamos cada variable segun el valor de los cuartiles.

library(dplyr)

p75 <-as.data.frame(t(apply(Matriz_media, 2, quantile, prob = 0.75)))
MatrizAltaExpr <- Matriz_media
MatrizAltaExpresionDef <- data.frame(row.names = rownames(Matriz_media))
for (genes in colnames(Matriz_media)){
  print (genes)
  MatrizAltaExpr2<-sweep(MatrizAltaExpr[genes],1, p75[[genes]], ">") #Aparece como "FALSE", no altamente expresado, "TRUE", altamente expresado.
  MatrizAltaExpresionDef<-cbind(MatrizAltaExpresionDef, MatrizAltaExpr2)
}
save.image("./20190801_filtrado.RData")

#------------GENES ALTAMENTE EXPRESADOS EN UN TEJIDO Y EN OTRO NO EXPRESADO--------------------
# Tienes que seleccionar mismo gen, en un tejido altamente expresado y en otro no expresado. Idea es fusionar las dos matrices anteriores y seleccionar.

# MatrizExpresion <- MatrizGenesNoExpr
MEDefinitiv<- MatrizGenesNoExpr
for (tejido in listado_tej) {
  print(tejido)
  #tejido="CENTRAL_NERVOUS_SYSTEM"
  for (genes in colnames(MatrizGenesNoExpr)){
    ifelse (
      is.na(MEDefinitiv[tejido,genes]),
      next,
      ifelse(
        MatrizAltaExpresionDef[tejido,genes] == "TRUE",
        MEDefinitiv[tejido,genes] <- "SI",
        MEDefinitiv[tejido,genes] <- "NO")
    )
  }
}

library(reshape2)
MEDefParaFiltarGenes<-melt(as.matrix(MEDefinitiv))
colnames(MEDefParaFiltarGenes)<-c("Tejido", "Gen", "ValorExpr")

library(plyr)
a1<-count(MEDefParaFiltarGenes, c('Gen','ValorExpr')) #Contamos por gen, cuanto SI, NO y NA hay.
a2<-dcast(a1, Gen ~ ValorExpr) #Volvemos a cambiar la forma de los datos. Le dices que te coja lo anterior y que gen este modulado por vlor expresion.
rownames(a2) <- a2[,1]
a2 <- a2[-c(1)]
a2[is.na(a2)] <- 0

#Creamos una función para que según el número de ceros, si o no, clasifique a los genes.(19 es el numero total de tejidos que tengo)
filtradoGenes <- function(x,y){ 
  #y=numero de tejidos
  GenesNoExpr <- c()
  GenesAlto <-c()
  GenesPocoExpr <- c()
  GenePocoYnoExpr <- c()
  GenesAltosYnoExpr <- c()
  GenesAltoYBajaExpr <- c()
  OtrosGenes <- c()
  #x=e1
  #genes="FUCA2"
  for (genes in (rownames(x))){
    if (x[genes,"NA"] == y) {
      GenesNoExpr <- c(GenesNoExpr, genes)
    }else{
      if (x[genes,"SI"] == y) {
        GenesAlto <- c(GenesAlto, genes)
      }else {
        if (x[genes, "NO"] == y){
          GenesPocoExpr <- c(GenesPocoExpr, genes)
        }else {
          if (x[genes, "NA"] >= 1 & x[genes, "NO"] >= 1 & (x[genes, "NO"] + x[genes, "NA"] == y)){
            GenePocoYnoExpr <- c(GenePocoYnoExpr, genes)
          }else{
            if (x[genes, "NA"] >= 1 & x[genes, "SI"] >= 1 & (x[genes, "SI"] + x[genes, "NA"] == y)){
              GenesAltosYnoExpr <- c(GenesAltosYnoExpr, genes)
            }else{
              if (x[genes, "SI"] >= 1 & x[genes, "NO"] >= 1 & (x[genes, "SI"] + x[genes, "NO"] == y)){
                GenesAltoYBajaExpr <- c(GenesAltoYBajaExpr, genes)
              }else{
                OtrosGenes <- c(OtrosGenes, genes)
              }
            }
          }
        }
      }
    } 
  }
  listaGenes<-list(GenesNoExpr,GenesAlto, GenesPocoExpr, GenePocoYnoExpr, GenesAltoYBajaExpr,GenesAltosYnoExpr,OtrosGenes,y) #Cuando haces una funcion te devuelve la ultima linea por eso, hacemos una lista con los cuatro vectores y le decimos con el return que sea lo que devuelva.
  names(listaGenes)<-c("GenesNoExpr","GenesAlto", "GenesPocoExpr", "GenePocoYnoExpr", "GenesAltoYBajaExpr","GenesAltosYnoExpr","OtrosGenes","NumeroTejidos")
  return(listaGenes)
}
ntejidos <- 23
ClasificacionGENES<-filtradoGenes(a2, ntejidos)

TablaClasificacionGenes <- data.frame(length(ClasificacionGENES$GenesNoExpr),length(ClasificacionGENES$GenesAlto),
                                      length(ClasificacionGENES$GenesPocoExpr), length(ClasificacionGENES$GenePocoYnoExpr),
                                      length(ClasificacionGENES$GenesAltoYBajaExpr), length(ClasificacionGENES$GenesAltosYnoExpr), 
                                      length(ClasificacionGENES$OtrosGenes))

colnames(TablaClasificacionGenes) <-names(ClasificacionGENES[-c(8)])
rownames(TablaClasificacionGenes)<- "NumeroGenes"
TablaClasificacionGenes <-t(TablaClasificacionGenes)
TablaClasificacionGenes <- as.data.frame((TablaClasificacionGenes))

ClasificacionGenes <- c("Genes no expresados en ningun tejido", "Genes con alta expresion en todos los tejidos", 
                        "Genes poco expresados en todos los tejidos", "Genes poco expresados en algun tejido y no expresados en otros", 
                        "Genes con alta expresion en algunos tejidos y baja expresion en otros", "Genes con alta expresion en algunos tejidos y no expresados en otros",
                        "Otros genes")
#barplot(TablaClasificacionGenes$NumeroGenes)
library(ggplot2)
library("RColorBrewer")
display.brewer.all()
#colores<-viridisLite::viridis(7)
colores <-brewer.pal(n = 7, name = "Spectral")
GENES <- rownames(TablaClasificacionGenes)
png("ClasificacionGenes2.png", width = 10, height = 10, units = "in", res = 110)
ggplot(TablaClasificacionGenes, aes(rownames(TablaClasificacionGenes), (NumeroGenes), fill = GENES))+
    geom_bar(stat = "identity") + 
    scale_fill_manual (values=c(colores))+
    # theme_minimal()  +
    ylab("Número de genes en cada grupo")+
    theme(legend.position="top", legend.key.size = unit(0.8, "cm"), #legend.title = element_text(size=10, face="bold"),
      axis.title.x = element_blank(), legend.title=element_blank(), axis.title.y = element_text(size=24, face="bold"),
      legend.background = element_rect(fill="white", size=1.5, linetype="solid"),
      panel.background = element_rect(fill='white', colour='white'), axis.text.x = element_text(angle = 50, hjust = 1, face="bold", size=16),
      axis.text.y=element_text(face="bold", size=14), legend.text=element_text(size=14),
      axis.line = element_line(color = "black", size = 0.5, linetype = "solid"))
dev.off()

save.image("./20190801_filtrado.RData")

############################################################################################################################################
#################################### Viendo si se agrupan por tejidos las muestras  ########################################################
############################################################################################################################################

resumen <- data.frame (apply(Matriz_media, 2, summary))

library(magrittr)
# library(tidyverse)
library(dplyr)
library(ggplot2)

res2<-resumen %>% t() %>% as.data.frame() #traspones la matriz resumen y la pasas a data.frame
hist(res2$Mean)
hist(res2$Mean,breaks = 40)

############# Representación media vs mediana #####################################################################
res3<-res2 %>% arrange(Mean) #Pierdes con arrange el nombre de los rownames
res3a<-res2[order(res2$Mean),]
png("MediavsMediana.png", width = 7, height = 7, units = "in", res = 110)
ggplot(res2, aes(x=Mean,y=Median)) +
  geom_point()+
  theme(axis.title.x =element_text(size=16, face="bold"), axis.title.y = element_text(size=16, face="bold"),
        axis.text.x = element_text(face="bold", size=14),
        axis.text.y=element_text(face="bold", size=14), title =element_text(size=18, face="bold"))+
  labs(title="Media vs mediana", x ="Media", y = "Mediana")
dev.off()

#Representacion de las medias de cada gen en todos los tejidos (histograma)
png("Media_histograma.png", width = 7, height = 7, units = "in", res = 110)
ggplot(res3, aes(x=Mean)) +
  geom_histogram(bins = 40,colour="white")+
  theme(axis.title.x =element_text(size=16, face="bold"), axis.title.y = element_text(size=24, face="bold"),
        axis.text.x = element_text(face="bold", size=14),
        axis.text.y=element_text(face="bold", size=14), title =element_text(size=18, face="bold"))+
  labs(title="Distribución de la media de los genes", x ="Valores de expresión", y = "")
dev.off()


######Formas de clusterizar####
###   PCA  #####
library(tidyverse)
#El siguiente paso solo es para representar los gáficos en español.
a1 <- as.data.frame(lapply(MatrizFiltrada, gsub, pattern = "OVARY", replacement = "OVARIO", fixed = TRUE))
a2 <- as.data.frame(lapply(a1, gsub, pattern = "SKIN", replacement = "PIEL", fixed = TRUE))
a2 <- as.data.frame(lapply(a2, gsub, pattern = "KIDNEY", replacement = "RIÑON", fixed = TRUE))
a2 <- as.data.frame(lapply(a2, gsub, pattern = "CENTRAL_NERVOUS_SYSTEM", replacement = "SISTEMA_NERVIOSO_CENTRAL", fixed = TRUE))
a2 <- as.data.frame(lapply(a2, gsub, pattern = "FIBROBLAST", replacement = "FIBROBLASTO", fixed = TRUE))
a2 <- as.data.frame(lapply(a2, gsub, pattern = "AUTONOMIC_GANGLIA", replacement = "GANGLIOS_AUTONOMICOS", fixed = TRUE))
a2 <- as.data.frame(lapply(a2, gsub, pattern = "BILIARY_TRACT", replacement = "TRACTO_BILIAR", fixed = TRUE))
a2 <- as.data.frame(lapply(a2, gsub, pattern = "OESOPHAGUS", replacement = "ESOFAGO", fixed = TRUE))
a2 <- as.data.frame(lapply(a2, gsub, pattern = "HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", replacement = "TEJIDO_HEMATOPOYETICO_Y_LINFOIDE", fixed = TRUE))
a2 <- as.data.frame(lapply(a2, gsub, pattern = "LUNG", replacement = "PULMON", fixed = TRUE))
a2 <- as.data.frame(lapply(a2, gsub, pattern = "BREAST", replacement = "MAMA", fixed = TRUE))
a2 <- as.data.frame(lapply(a2, gsub, pattern = "THYROID", replacement = "TIROIDE", fixed = TRUE))
a2 <- as.data.frame(lapply(a2, gsub, pattern = "STOMACH", replacement = "ESTOMAGO", fixed = TRUE))
a2 <- as.data.frame(lapply(a2, gsub, pattern = "ENDOMETRIUM", replacement = "ENDOMETRIO", fixed = TRUE))
a2 <- as.data.frame(lapply(a2, gsub, pattern = "LIVER", replacement = "HIGADO", fixed = TRUE))
a2 <- as.data.frame(lapply(a2, gsub, pattern = "LARGE_INTESTINE", replacement = "INTESTINO_DELGADO", fixed = TRUE))
a2 <- as.data.frame(lapply(a2, gsub, pattern = "URINARY_TRACT", replacement = "TRACTO_URINARIO", fixed = TRUE))
a2 <- as.data.frame(lapply(a2, gsub, pattern = "BONE", replacement = "HUESO", fixed = TRUE))
a2 <- as.data.frame(lapply(a2, gsub, pattern = "SOFT_TISSUE", replacement = "TEJIDO_BLANDO", fixed = TRUE))
a2 <- as.data.frame(lapply(a2, gsub, pattern = "PROSTATE", replacement = "PROSTATA", fixed = TRUE))
a2 <- as.data.frame(lapply(a2, gsub, pattern = "UPPER_AERODIGESTIVE_TRACT", replacement = "TRACTO_AERODIGESTIVO_SUPERIOR", fixed = TRUE))

write.table (a2, file = "MatrizFiltradaEspanol.csv", col.names = TRUE, sep="\t")
MatrizFiltradaEspanol <- read.table ('./MatrizFiltradaEspanol.csv', header = TRUE, sep= '\t')
#p1 <- MatrizFiltrada[,apply(MatrizFiltrada, 2, function(x) all(x > 0))] Esta es la línea si no se traducen los tejidos al castellano
p1 <- MatrizFiltradaEspanol[,apply(MatrizFiltradaEspanol, 2, function(x) all(x > 0))] 
pca.result <- prcomp(p1[,-c(1:3)],center = TRUE, scale. = TRUE)
head(summary(pca.result))
dim(pca.result$x)
dim(pca.out$x)

#################### -----PCA, sin hemato y pulmon-----   #######################
library(dplyr)
MatrizFiltradaSinH <- MatrizFiltrada%>%filter(!(TISSUE=="HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"))
MFsinHyL<-MatrizFiltradaSinH%>%filter(!(TISSUE=="LUNG"))
p2 <- MFsinHyL[,apply(MFsinHyL, 2, function(x) all(x > 0))] 
pca.result2 <- prcomp(p2[,-c(1:3)],center = TRUE, scale. = TRUE)

#################################################################################################################################################
########################################## PCA representado por ggplot  ##########################################################################
##################################################################################################################################################
pca.result3 <- as.data.frame(pca.result$x)
#pca.result3 <- cbind(MatrizFiltrada$TISSUE,pca.result3) si los tejidos no estan en casterllano
pca.result3 <- cbind(MatrizFiltradaEspanol$TISSUE,pca.result3)
colnames(pca.result3)[1]<-"TISSUE"
pca.result3$TISSUE<-as.character(pca.result3$TISSUE)
#pca.result3 <- cbind(MatrizFiltrada$CCLE, pca.result3) si los tejidos no están en castellano
pca.result3 <- cbind(MatrizFiltradaEspanol$CCLE, pca.result3)

#Tejidos <- MatrizFiltrada$TISSUE si los tejidos no estan en castellano
Tejidos <- MatrizFiltradaEspanol$TISSUE
rainbowcols <- rainbow(5, s = 0.9)
colores<-colorRampPalette(rainbowcols)(24)
names(colores)<-unique(Tejidos)

library(ggplot2)
library(ggrepel)

percentage <- round(pca.result$sdev / sum(pca.result$sdev) * 100, 2) #Calula los porcentajes de cada CPs.
TEJIDOS <- pca.result3$TISSUE
png("PCA2.png", width = 25, height = 15, units = "in", res = 110)
ggplot() +
  geom_point(data = pca.result3, aes(x = pca.result3$PC1, y = pca.result3$PC2, colour = TEJIDOS), size=3) +
  theme(legend.position="top", legend.key.size = unit(1.0, "cm"),#legend.title = element_text(size=6, face="bold"),
        axis.title.x = element_text(size=38, face="bold"), legend.title=element_blank(), axis.title.y = element_text(size=35, face="bold"),
        legend.background = element_rect(fill="lightgrey", size=4.5, linetype="solid"),
        panel.background = element_rect(fill='darkgrey', colour='darkgrey'), legend.text=element_text(size=22, face = 'bold')) +
  guides(colour = guide_legend(override.aes = list(size=16)))+
  #guides(fill = guide_legend(nrow = 7))+
  #geom_text_repel(aes(x = pca.result3[(pca.result3$PC1>175),]$PC1, y = pca.result3[(pca.result3$PC1>175),]$PC2,
                      #colour = pca.result3[(pca.result3$PC1>175),]$TISSUE),
                  #label=pca.result3[(pca.result3$PC1>175),]$`MatrizFiltrada$CCLE`, size=7.5, face="bold", color='black') +
  scale_colour_manual(values=colores) +
  xlab(paste("PC1 (", paste(as.character(percentage[1]), "%)", sep="") )) +
  ylab(paste("PC2 (", paste(as.character(percentage[2]), "%)", sep="") ))
dev.off()
#save.image("./20190801_filtrado.RData")

######## PCA sin hemato y Lung  #####
pca.result4 <- as.data.frame(pca.result2$x)
pca.result4 <- cbind(MFsinHyL$TISSUE, pca.result4)
colnames(pca.result4)[1]<-"TISSUE"
pca.result4$TISSUE<-as.character(pca.result4$TISSUE)
pca.result4 <- cbind(MFsinHyL$CCLE, pca.result4)

library(RColorBrewer)
Tejidos3 <- MFsinHyL$TISSUE
rainbowcols <- rainbow(5, s = 0.9)
colores2<-colorRampPalette(rainbowcols)(22)
names(colores2)<-unique(Tejidos3)

percentage2 <- round(pca.result2$sdev / sum(pca.result2$sdev) * 100, 2) #Calula los porcentajes de cada CPs.
TEJIDOS <- pca.result4$TISSUE
png("PCA_sinHematoYPumon2.png", width = 15, height = 15, units = "in", res = 110)
ggplot() +
  geom_point(data = pca.result4, aes(x = pca.result4$PC1, y = pca.result4$PC2, colour = TEJIDOS), size=3.5) +
  theme(legend.position="top", legend.key.size = unit(0.8, "cm"),#legend.title = element_text(size=6, face="bold"),
        axis.title.x = element_text(size=22, face="bold"), legend.title=element_blank(), axis.title.y = element_text(size=22, face="bold"),
        legend.background = element_rect(fill="lightgrey", size=0.5, linetype="solid"),
        panel.background = element_rect(fill='darkgrey', colour='darkgrey'), legend.text=element_text(size=14, face = 'bold')) +
  geom_text_repel(aes(x = pca.result4[(pca.result4$PC1>200),]$PC1, y = pca.result4[(pca.result4$PC1>200),]$PC2, 
                      colour = pca.result4[(pca.result4$PC1>200),]$TISSUE),
                  label=pca.result4[(pca.result4$PC1>200),]$`MFsinHyL$CCLE`, size=5.5, face="bold", color='black') +
  guides(colour = guide_legend(override.aes = list(size=10)))+
  scale_colour_manual(values=colores2) +
  xlab(paste("PC1 (", paste(as.character(percentage2[1]), "%)", sep="") )) +
  ylab(paste("PC2 (", paste(as.character(percentage2[2]), "%)", sep="") )) 
dev.off()

#########################################################################################################################################
############################################################### CLUSTER #################################################################
#########################################################################################################################################


clust.med <- hclust ( dist (MatrizFiltrada, method="euclidean" ) , method = "average") #Calcula en un objeto la distancia media entre observaciones para luego poder representarlo en un dendrograma.
dev.new(width=300, height=100)
plot( clust.med, main="Análisis Cluster jerárquico", labels =MatrizFiltrada$TISSUE, ylab="Distance", xlab="", sub="", cex=0.6)
png("cluster_jerarquico.png", width = 15, height = 15, units = "in", res = 110)
plot( clust.med, main="Análisis Cluster jerárquico", labels=FALSE, ylab="Distance", xlab="", sub="", cex=0.6)
dev.off()

library(dendextend)

dend2 <- clust.med %>% as.dendrogram
dend2 %>% plot #Te saca solo la parte del dendograma del clustering.

rainbowcols <- rainbow(5, s = 0.9)
colores<-colorRampPalette(rainbowcols)(24)
names(colores)<-unique(Tejidos)

#Le esta asociando a cada tejido un color diferente
col2<-as.data.frame(colores)
col2$Tissue<-rownames(col2)
col3<-as.data.frame(as.character(MatrizFiltrada[,-c(1,3:ncol(MatrizFiltrada))]))
colnames(col3)<-"Tissue"
col4<-merge(col3,col2, by="Tissue")

png("clustering_feo2.png", width = 7, height = 7, units = "in", res = 110)
dend2%>%
  set("labels_colors", "white") %>% 
  set("labels_cex", 0.01) %>%
  plot
colored_bars(colors = col4$colores, dend = dend2%>%set("labels_colors", "white"), sort_by_labels_order = T, size=30)
dev.off()

###########################################################################################################################################
################################### K-means (k=23, que son los tejidos)####################################################################
###########################################################################################################################################

#km.out <- kmeans(MatrizFiltrada[,-c(1,2,3)] , cent=23 )
km.out <- kmeans(MatrizFiltradaEspanol[,-c(1,2,3)] , cent=23 )
dim(km.out)
names (km.out)
names(MatrizFiltrada_clustered)
#MatrizFiltrada_clustered <- data.frame(MatrizFiltrada, cluster=factor(km.out$cluster))
MatrizFiltrada_clustered <- data.frame(MatrizFiltradaEspanol, cluster=factor(km.out$cluster))
save.image("./20190801_filtrado.RData")

#Contamos el número que da por cluster para compararlo con el número de lineas por tejido y ver cuanto aciertos hay.
library(dplyr)
listado_tej <- unique(MatrizFiltrada_clustered$TISSUE)
MatrizCuentas <- data.frame(MatrizFiltrada_clustered$TISSUE, MatrizFiltrada_clustered$cluster)
numero_tejidos<- c(1:length(listado_tej))
colnames(MatrizCuentas)<-c('Tissue','cluster')

m<-data.frame(matrix(ncol = length(numero_tejidos), nrow = length(listado_tej)))
colnames(m)<-as.character((1:length(listado_tej)))
rownames(m)<-sort(listado_tej, decreasing = FALSE)

for (tejido in listado_tej){
  print(tejido)
  z1 <-MatrizCuentas %>% filter(Tissue==tejido) #Agrupo por tejido
  for (numero in numero_tejidos){
    numeroClasificacion<-colSums(z1==numero) %>% .[-c(1)]
    m[tejido, numero]<-numeroClasificacion
    }
}

tablaFreqTej <- data.frame(table(MatrizFiltrada_clustered$TISSUE))
v <- c(subset(tablaFreqTej, Freq > 5)) #Listado de los tejidos que con mas de 5 lineas por tejido.
v <- data.frame(subset(tablaFreqTej, Freq > 5))

ListaTejidoFreq<-v
colnames(ListaTejidoFreq) <- c("TISSUE", "Nº_lineas")
nº_lineas<-ListaTejidoFreq[,'Nº_lineas']
Mkmeans <-cbind(nº_lineas, m)

#########  Representacion numero tejidos por linea ###############
library(ggplot2)
library(RColorBrewer)
library(ggsci)

Tejidos2 <- rownames(Mkmeans)
Tejidos2[20]<-"TEJIDO_HEMATOPOYETICO\nY_LINFOIDE"
Tejidos2[6]<-"GANGLIOS\nAUTONOMICOS"
Tejidos2[18]<-"SISTEMA_NERVIOSO\nCENTRAL"
Tejidos2[22]<-"TRACTO_AERODIGESTIVO\nSUPERIOR"
rainbowcols <- rainbow(5, s = 0.9)
colores2<-colorRampPalette(rainbowcols)(24)
names(colores2)<-unique(Tejidos2)


#Ejecutar si los nombres de los tejidos estan en ingles ###
Tejidos2 <- rownames(Mkmeans)
Tejidos2[9]<-"HAEMATOPOIETIC\nAND_LYMPHOID_TISSUE"
Tejidos2[1]<-"AUTONOMIC\nGANGLIA"
Tejidos2[5]<-"CENTRAL_NERVOUS\nSYSTEM"
Tejidos2[23]<-"UPPER\nAERODIGESTIVE_TRACT"
colores2<-colorRampPalette(rainbowcols)(24)
names(colores2)<-unique(Tejidos2)


png("Distribucion_tejidos.png", width = 20, height = 10, units = "in", res = 110)
ggplot(Mkmeans, aes(Tejidos2, (nº_lineas), fill = Tejidos2)) +
  geom_bar(stat = "identity") + 
  ylab("Número de lineas en cada tejido")+
  scale_fill_manual (values=colores2)+
  theme(legend.position="none", #legend.key.size = unit(0.8, "cm"), legend.title = element_text(size=6, face="bold"),
        axis.title.x = element_blank(), legend.title=element_blank(), axis.title.y = element_text(size=25, face="bold"),
        legend.background = element_rect(fill="lightgrey", size=4.5, linetype="solid"),
        panel.background = element_rect(fill='darkgrey', colour='darkgrey'), axis.text.x = element_text(angle = 65, hjust = 1, face="bold", size=22),
        axis.line = element_line(color = "darkblue", size = 0.5, linetype = "solid"), 
        axis.text.y = element_text(face="bold", size=20))
dev.off()

#Representacion del k-means- En un eje, el dendograma de como se ajustan los cluster segun k-means, y el otro eje los tejidos. Lo ideal si cada clustr fuera una tejido es la diagonal de un color y el resto de otro.
#heatmap
install.packages("gplots")
library(RColorBrewer)
display.brewer.all()
col=colorRampPalette(brewer.pal(9,"RdYlGn"))

# # probando colores
png("Heatmap_cluster2.png", width = 30, height = 15, units = "in", res = 110)
heatmap.2(as.matrix(Mkmeans[,-c(1)]),
          trace="none",
          #breaks = seq(0,100,2),
          # col=col2,
          col=col2,
          cexRow=2,
          cexCol = 2,
          margins = c(10, 40),
          keysize = 1.2,
          srtCol=0
)
dev.off()
library(gplots)
col<- colorRampPalette(brewer.pal(9,"Blues"))(200)
col2<-c(col3[1] , col[2:200])
col3<-c("#fbfbfb")
m <- as.data.frame(Mkmeans/(Mkmeans[,1]))
getwd()
setwd("/home/estrella/TFM/Datos/depmap_19Q1_relabeled/FiltradoDatosYAnalisisDatosCrudos_(Agrupado por tejido)/")

png("Heatmap_cluster_normalizado2_2.png", width = 30, height = 15, units = "in", res = 110)
heatmap.2(as.matrix(m[,-c(1)]),
          trace="none", 
          key=TRUE,
          Colv=FALSE, #Quita el dendograma de las columnas
          #Rowv=FALSE,
          dendrogram = "row",
          #breaks = seq(0,200,0.1),
          col=col2,
          cexRow=2.3,
          cexCol = 2.8,
          margins = c(5, 42),
          keysize =0.85,
          srtCol=5,
          key.xlab="Valor", key.ylab="", key.title="", 
          key.par=list(mar=c(4,1,1,2), cex=1.0, cex.lab=2.0, cex.axis=2.0)      
)
dev.off()

#Representacion del k-means con barplot (cada barra, representa el total de lineas por tejido. Cada color son las lineas asignada a u cluster en concreto.)

rainbowcols <- rainbow(5, s = 0.9)
colores3<-colorRampPalette(rainbowcols)(23)

MatrizCuentas2<-MatrizCuentas
MatrizCuentas2$total<-"1"
MatrizCuentas2$cluster
x <- MatrizCuentas2$Tissue
Tissue2 <- as.data.frame(recode(x, "HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"="HAEMATOPOIETIC\nAND_LYMPHOID_TISSUE", "AUTONOMIC_GANGLIA"="AUTONOMIC\nGANGLIA", "CENTRAL_NERVOUS_SYSTEM"="CENTRAL_NERVOUS\nSYSTEM", "UPPER_AERODIGESTIVE_TRACT"="UPPER\nAERODIGESTIVE_TRACT"))
colnames(Tissue2)<-"Tissue2"
MatrizCuentas2 <- cbind(MatrizCuentas2, Tissue2)

png("Distribucion_tejidos_cluster_barplot.png", width = 22, height = 12, units = "in", res = 110)
ggplot(MatrizCuentas2, aes(Tissue2, total, fill = cluster))+
  geom_bar(stat = "identity") + 
  theme_minimal()  +
  ylab("NUMERO DE LINEAS ASIGNADAS A CADA CLUSTER")+
  labs (fill = "CLUSTERS")+
  scale_fill_manual (values=c(colores3))+
  theme(legend.position="right", legend.title = element_text(size=16, face="bold"),
        axis.title.x = element_blank(), axis.title.y = element_text(size=18, face="bold"),
        legend.background = element_rect(fill="white", size=0.8, linetype="solid"),
        panel.background = element_rect(fill='white', colour='darkgrey'), axis.text.x = element_text(angle = 65, hjust = 1, face="bold", size=16),
        axis.text.y=element_blank(), legend.text=element_text(size=14, face = 'bold'),
        axis.line = element_line(color = "darkblue", size = 0.5, linetype = "solid"))
dev.off()
getwd()
save.image("./20190801_filtrado.RData")
