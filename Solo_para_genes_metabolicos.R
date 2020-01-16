##############################################################################################################
##########################REPETICION DEL ANALISIS SOLO PARA GENES METABOLICOS ################################
##############################################################################################################

getwd()
setwd("/home/estrella/TFM/Datos/depmap_19Q1_relabeled/GenesMetabolicos_(Agrupado por tejido)/")
#matrizDatosCrudos<- read.table('/home/masterbioinfo/M19Q2_metabol.csv', header = TRUE, sep= '\t')
matrizDatosCrudos<- read.table ('../M19Q2_metabol.csv', header = TRUE, sep= '\t')
library(dplyr)
library(reshape2)
div<-colsplit(matrizDatosCrudos$CCLE, "_", c("ID","TISSUE")) #Divides la columna cell_line(que tiene ID_tej) en dos columnas.
MatrizDatos <- cbind(div,matrizDatosCrudos) #Le unes lo que acabas de separar.

#------ AGRUPAR POR TEJIDO-----------------
unique(MatrizDatos$TISSUE)
tablaFreqTej <- data.frame(table(MatrizDatos$TISSUE))

v <- c(subset(tablaFreqTej, Freq > 5)) #Listado de los tejidos que con mas de 5 lineas por tejido.
v <- data.frame(subset(tablaFreqTej, Freq > 5))
MatrizMetabolismo <-MatrizDatos[MatrizDatos$TISSUE%in%as.character(v$Var1),]
unique(MatrizFiltrada$TISSUE)
#Esto solo es para cambiar los nombres de los tejidos a castellano


#Grafico de numero de lineas por tejido(es igual que la matriz total).
par(mfrow=c(2,1.9))
barplot (prop.table(table(MatrizMetabolismo$TISSUE)), col=c("blue"), main="Lineas celulares agrupadas por tejidos", ylim=c(0,0.2), ylab = "Frecuencias", las=2, cex.lab=1)
barplot (table(MatrizMetabolismo$TISSUE), col=c("blue"), main="Lineas celulares agrupadas por tejidos", ylab = "Nº lineas", las=2, cex.lab=1)

#--------MEDIA Y MEDIANA DE LOS GENES SEGUN TEJIDO-------------
########  MEDIA  ######for (tejido in listado_tej){ #DESDE AQUI!!!!!!!!
listado_tej<- as.character(v$Var1) #Listado de los tejidos que me valen
MetabolGenesMedia <-data.frame()
for (tejido in listado_tej) {
  print(tejido)
  mtejido<-MatrizMetabolismo[which(MatrizMetabolismo$TISSUE==tejido),][,4:length(MatrizMetabolismo)]
  media_genes<-c()
  media_genes<-apply(mtejido,2,mean)
  media_genes<-t(as.data.frame(media_genes))
  rownames(media_genes)<-tejido
  MetabolGenesMedia<-rbind(MetabolGenesMedia, media_genes)
}

########  MEDIANA  ######

MetabolGenesMediana <-data.frame()
for (tejido in listado_tej) {
  print(tejido)
  mtejido<-MatrizMetabolismo[which(MatrizMetabolismo$TISSUE==tejido),][,4:length(MatrizMetabolismo)]
  mediana_genes<-c()
  mediana_genes<-apply(mtejido,2,median)
  mediana_genes<-t(as.data.frame(mediana_genes))
  rownames(mediana_genes)<-tejido
  MetabolGenesMediana<-rbind(MetabolGenesMediana, mediana_genes)
}

#-------ELIMINANDO GENES NO EXPRESADOS (GENES QUE SE EXPRESEN EN MENOS DE UN 10% DE LAS LINEAS)-----------------------------
library(dplyr)
listado_genesMetabol <- colnames(MatrizMetabolismo[-c(1)])

MGMetabolNoExpresados<-data.frame() #Cuenta el numero de ceros. Numero de veces que no se expresa un gen.
for (tejido in listado_tej){ 
  print(tejido)
  t1 <-MatrizMetabolismo%>%filter(TISSUE==tejido) #Agrupo por tejido
  vacio_tejido<-cbind(tejido, colSums(t1==0)  %>% .[-(1:3)] %>% as.data.frame() %>% t())
  MGMetabolNoExpresados <-rbind(MGMetabolNoExpresados,vacio_tejido)
}

colnames(tablaFreqTej) <- c("tejido", "Freq")
MGMetabolNoExpresadosDef <- merge(tablaFreqTej, MGMetabolNoExpresados, by="tejido") #Nos uniría la columna de frecuencias con el numero total de lineas por tejido.
rownames(MGMetabolNoExpresadosDef)<-MGMetabolNoExpresadosDef$tejido #Cambiamos el nombre a rownames la columna tejido. Para que solo se queden las filas y columnas con las que queremos operar.
MGMetabolNoExpresadosDef<-MGMetabolNoExpresadosDef[-c(1)]
str(MGMetabolNoExpresadosDef) #Tenemos como factores las lineas. Por eso le quitamos la primera columna (tej), y el resto que lo transforme a numerico. y lo guarde como data frame.
p1<-apply(MGMetabolNoExpresadosDef,2,as.numeric) %>% as.data.frame() #Al hacer esto te quita los row names por eso hacemos lo siguiente.
any(is.na(p1))
rownames(p1)<-rownames(MGMetabolNoExpresadosDef)
frecuenciasNoExpr <- p1[-c(1)]/p1$`Freq` 
frecuenciasNoExpr[frecuenciasNoExpr > 0.1] <- NA #Cambio los valores de frecuencia que sean mayores a 0.1 por NA. Si obtengo una proporcion de ceros mayor del 0.9, quiere decir que solo se expresa en el 10% de las lineas y lo considero como no expresado.
MatrizGenesMetabolNoExpr <-frecuenciasNoExpr

#------------IDENTIFICAMOS GENES ALTAMENTE EXPRESADOS EN TODOS LOS TEJIDOS--------------
# Definimos los genes altamente expresados como aquellos que estan por encima del cuartil 75%.
#De la matriz que contiene las medias de cada gen por tejido (media_genes2). Hacer de cada columna (gen) el cuartil. Categorizamos cada variable segun el valor de los cuartiles.

library(dplyr)

p75 <-as.data.frame(t(apply(MetabolGenesMedia, 2, quantile, prob = 0.75)))
MetabolAltaExpr <- MetabolGenesMedia
MetabolAltaExprDef <- data.frame(row.names = rownames(MetabolGenesMedia))
for (genes in colnames(MetabolGenesMedia)){
  print (genes)
  MetabolAltaExpr2<-sweep(MetabolAltaExpr[genes],1, p75[[genes]], ">") #Aparece como "FALSE", no altamente expresado, "TRUE", altamente expresado.
  MetabolAltaExprDef<-cbind(MetabolAltaExprDef, MetabolAltaExpr2)
}


#------------GENES ALTAMENTE EXPRESADOS EN UN TEJIDO Y EN OTRO NO EXPRESADO--------------------
# Tienes que seleccionar mismo gen, en un tejido altamente expresado y en otro no expresado. Idea es fusionar las dos matrices anteriores y seleccionar.

MetabolExprDef<-MatrizGenesMetabolNoExpr
for (tejido in listado_tej){
  print(tejido)
  for (genes in colnames(MetabolExprDef)){
    # print (genes)
    # genes= "CFTR"
    ifelse (
      is.na(MetabolExprDef[tejido,genes]),
      next,
      ifelse(
        MetabolAltaExprDef[tejido,genes] == "TRUE",
        MetabolExprDef[tejido,genes] <- "SI",
        MetabolExprDef[tejido,genes] <- "NO")
      )
  }
  }

library(reshape2)
MEDefParaFiltarGenes<-melt(as.matrix(MetabolExprDef))
colnames(MEDefParaFiltarGenes)<-c("Tejido", "Gen", "ValorExpr")
str(MEDefParaFiltarGenes)

library(plyr)
a1<-count(MEDefParaFiltarGenes, c('Gen','ValorExpr')) #Contamos por gen, cuanto SI, NO y NA hay.

#library(reshape2) 
a2<-dcast(a1, Gen ~ ValorExpr) #Volvemos a cambiar la forma de los datos. Le dices que te coja lo anterior y que gen este modulado por vlor expresion.
rownames(a2) <- a2[,1]
a2 <- a2[-c(1)]
a2[is.na(a2)] <- 0

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
ClasificacionGENESMetabol<-filtradoGenes(a2, ntejidos)

TablaClasificacionGenes <- data.frame(length(ClasificacionGENESMetabol$GenesNoExpr),length(ClasificacionGENESMetabol$GenesAlto),
                                      length(ClasificacionGENESMetabol$GenesPocoExpr), length(ClasificacionGENESMetabol$GenePocoYnoExpr),
                                      length(ClasificacionGENESMetabol$GenesAltoYBajaExpr), length(ClasificacionGENESMetabol$GenesAltosYnoExpr), 
                                      length(ClasificacionGENESMetabol$OtrosGenes))

colnames(TablaClasificacionGenes) <-names(ClasificacionGENESMetabol[-c(8)])
rownames(TablaClasificacionGenes)<- "NumeroGenes"
TablaClasificacionGenes <-t(TablaClasificacionGenes)
TablaClasificacionGenes <- as.data.frame((TablaClasificacionGenes))

# porciones <- c(length(ClasificacionGENESMetabol$GenesNoExpr),length(ClasificacionGENESMetabol$GenesAltosYnoExpr),
#                length(ClasificacionGENESMetabol$GenesAlto), length(ClasificacionGENESMetabol$OtrosGenes),
#                length(ClasificacionGENESMetabol$GenesPocoExpr), length(ClasificacionGENESMetabol$GenePocoYnoExpr),
#                length(ClasificacionGENESMetabol$GenesAltoYBajaExpr))
ClasificacionGenes <- c("Genes no expresados en ningun tejido", "Genes con alta expresion en todos los tejidos", 
              "Genes poco expresados en todos los tejidos", "Genes poco expresados en algun tejido y no expresados en otroa", 
             "Genes con alta expresion en algunos tejidos y baja expresion en otros", "Genes con alta expresion en algunos tejidos y no expresados en otros",
             "Otros genes")

barplot(TablaClasificacionGenes$NumeroGenes)
library(ggplot2)
GENES <- rownames(TablaClasificacionGenes)
png("ClasificacionGenes2.png", width = 10, height = 10, units = "in", res = 110)
# ggplot(TablaClasificacionGenes, aes(rownames(TablaClasificacionGenes), log10(NumeroGenes), fill = GENES))+#, colour = "white"))+
colores<-viridisLite::viridis(7)
png("ClasificacionGenes2.png", width = 10, height = 10, units = "in", res = 110)
# ggplot(TablaClasificacionGenes, aes(rownames(TablaClasificacionGenes), log10(NumeroGenes), fill = GENES))+#, colour = "white"))+
ggplot(TablaClasificacionGenes, aes(rownames(TablaClasificacionGenes), (NumeroGenes), fill = GENES))+#, colour = "white"))+
  geom_bar(stat = "identity") + 
  scale_fill_manual (values=c(colores))+
  # theme_minimal()  +
  ylab("Número de genes en cada grupo")+
  theme(legend.position="top", legend.key.size = unit(0.8, "cm"), #legend.title = element_text(size=10, face="bold"),
        axis.title.x = element_blank(), legend.title=element_blank(), axis.title.y = element_text(size=24, face="bold"),
        legend.background = element_rect(fill="lightgrey", size=0.5, linetype="solid"),
        panel.background = element_rect(fill='lightgrey', colour='darkgrey'), axis.text.x = element_text(angle = 50, hjust = 1, face="bold", size=16),
        axis.text.y=element_text(face="bold", size=14), legend.text=element_text(size=14),
        axis.line = element_line(color = "darkblue", size = 0.5, linetype = "solid"))
dev.off()

save.image("./20190710_Metabolismo.RData")

############################################################################################################################################
#################################### Summary gen a gen sobre su media y su mediana  ########################################################
############################################################################################################################################

summary(MetabolGenesMedia$CYP51A1)
resumen <- data.frame (apply(MetabolGenesMedia, 2, summary))

library(magrittr)
# library(tidyverse)
library(dplyr)
library(ggplot2)

res2<-resumen %>% t() %>% as.data.frame() #traspnes la matriz resumen y la pasas a data.frame
hist(res2$Mean)
hist(res2$Mean,breaks = 40)

res3a<-res2[order(res2$Mean),]
#Representacion distribucion de medias frente a medianas, de cada gen en todos los tejidos.
png("MediavsMediana.png", width = 7, height = 7, units = "in", res = 110)
ggplot(res2%>% arrange(Median), aes(x=seq(1:length(res2$Median)),y=Median)) +
  geom_point()
ggplot(res2, aes(x=Mean,y=Median)) +
  geom_point()
dev.off()
#Representacion de las medias de cada gen en todos los tejidos (histograma)
png("Media_histograma.png", width = 7, height = 7, units = "in", res = 110)
ggplot(res3, aes(x=Mean)) +
  geom_histogram(bins = 40,colour="white")
dev.off()

######Formas de clusterizar####
###   PCA  #####
library(tidyverse)
#El siguiente paso solo es para representar los gáficos en español.
a1 <- as.data.frame(lapply(MatrizMetabolismo, gsub, pattern = "OVARY", replacement = "OVARIO", fixed = TRUE))
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

write.table (a2, file = "MatrizMetalbolFiltradaEspanol.csv", col.names = TRUE, sep="\t")
MatrizMetalbolFiltradaEspanol <- read.table ('./MatrizMetalbolFiltradaEspanol.csv', header = TRUE, sep= '\t')
#p1 <- MatrizMetabolismo[,apply(MatrizMetabolismo, 2, function(x) all(x > 0))] Esta es la línea si no se traducen los tejidos al castellano
p1 <- MatrizMetalbolFiltradaEspanol[,apply(MatrizMetalbolFiltradaEspanol, 2, function(x) all(x > 0))] 
pca.result <- prcomp(p1[,-c(1:3)],center = TRUE, scale. = TRUE)
head(summary(pca.result))
dim(pca.result$x)
dim(pca.out$x)

#################### -----PCA, sin hemato y pulmon-----   #######################
library(dplyr)
MatrizMetabolSinH <- MatrizMetabolismo%>%filter(!(TISSUE=="HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"))
MMsinHyL<-MatrizMetabolSinH%>%filter(!(TISSUE=="LUNG"))
pca.out3 <- prcomp( MMsinHyL[,-c(1:2)], scale=T )

################################################################################################################################################
######################################### PCA representado por plot, según el profe de bioestadistica ##########################################
################################################################################################################################################
Cols = function (vec)
   { cols = rainbow(length(unique((vec))))
   cols [ as.numeric ( as.factor (vec)) ]
}
tejidos<-as.character(MatrizMetabolismo$TISSUE)
plot ( pca.out$x[ , 1], pca.out$x[ , 2], xlab="CP1", ylab="CP2", cex=1.2, col=Cols(MatrizMetabolismo$TISSUE) , pch=16 ) #Representación gráfica de los CPAs. No parece que se vea ningún grupo separado en especial.
plot ( pca.out$x[ , 1], pca.out$x[ , 2], xlab="CP1", ylab="CP2", cex=1.2, col=Cols(tejidos) , pch=16 ) #Representación gráfica de los CPAs. No parece que se vea ningún grupo separado en especial.

#################################################################################################################################################
########################################## PCA representado por ggplot  ##########################################################################
##################################################################################################################################################
pca.result3 <- as.data.frame(pca.result$x)
#pca.result3 <- cbind(MatrizFiltrada$TISSUE,pca.result3) si los tejidos no estan en casterllano
pca.result3 <- cbind(MatrizMetalbolFiltradaEspanol$TISSUE,pca.result3)
colnames(pca.result3)[1]<-"TISSUE"
pca.result3$TISSUE<-as.character(pca.result3$TISSUE)
#pca.result3 <- cbind(MatrizFiltrada$CCLE, pca.result3) si los tejidos no están en castellano
pca.result3 <- cbind(MatrizMetalbolFiltradaEspanol$CCLE, pca.result3)

#Tejidos <- MatrizFiltrada$TISSUE si los tejidos no estan en castellano
Tejidos <- MatrizMetalbolFiltradaEspanol$TISSUE
rainbowcols <- rainbow(5, s = 0.9)
colores<-colorRampPalette(rainbowcols)(24)
names(colores)<-unique(Tejidos)

library(ggplot2)
library(ggrepel)

percentage <- round(pca.result$sdev / sum(pca.result$sdev) * 100, 2) #Calula los porcentajes de cada CPs.
TEJIDOS <- pca.result3$TISSUE
png("PCA2_metabol.png", width = 25, height = 15, units = "in", res = 110)
ggplot()+
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


######## PCA sin hemato y Lung######
pca.out4 <- as.data.frame(pca.out3$x)
pca.out4 <- cbind(MMsinHyL$TISSUE,pca.out4)
colnames(pca.out4)[1]<-"TISSUE"
pca.out4$TISSUE<-as.character(pca.out4$TISSUE)

MatrizFiltradaSinH <- MatrizFiltrada%>%filter(!(Tissue=="HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"))
MFsinHyL<-MatrizFiltradaSinH%>%filter(!(Tissue=="LUNG"))

library(RColorBrewer)
Tejidos3 <- MFsinHyL$Tissue
rainbowcols <- rainbow(5, s = 0.9)
colores2<-colorRampPalette(rainbowcols)(21)
names(colores2)<-unique(Tejidos3)

percentage2 <- round(pca.out3$sdev / sum(pca.out3$sdev) * 100, 2) #Calula los porcentajes de cada CPs.
TEJIDOS <- pca.out4$TISSUE
png("PCA_sinHematoYPumon.png", width = 15, height = 15, units = "in", res = 110)
ggplot() +
  geom_point(data = pca.out4, aes(x = pca.out4$PC1, y = pca.out4$PC2, colour = TEJIDOS), size=3.5) +
  theme(legend.position="top", legend.key.size = unit(0.8, "cm"), #legend.title = element_text(size=6, face="bold"),
        axis.title.x = element_text(size=12, face="bold"), legend.title=element_blank(), axis.title.y = element_text(size=22, face="bold"),
        legend.background = element_rect(fill="lightgrey", size=0.5, linetype="solid"),
        panel.background = element_rect(fill='darkgrey', colour='darkgrey'), legend.text=element_text(size=14, face = 'bold')) +
  guides(colour = guide_legend(override.aes = list(size=10)))+
  scale_colour_manual(values=colores2) +
  xlab(paste("PC1 (", paste(as.character(percentage2[1]), "%)", sep="") )) +
  ylab(paste("PC2 (", paste(as.character(percentage2[2]), "%)", sep="") )) 
dev.off()

## CLUSTER ## ----
clust.med <- hclust ( dist (MatrizMetabolismo, method="euclidean" ) , method = "average") #Calcula en un objeto la distancia media entre observaciones para luego poder representarlo en un dendrograma.
dev.new(width=300, height=100)
plot( clust.med, main="Análisis Cluster jerárquico", labels =MatrizMetabolismo$TISSUE, ylab="Distance", xlab="", sub="", cex=0.6)
plot( clust.med, main="Análisis Cluster jerárquico", labels=FALSE, ylab="Distance", xlab="", sub="", cex=0.6)
plot(clust.med)
colored_bars(colores, clust.med)

library(dendextend)

dend2 <- clust.med %>% as.dendrogram
dend2 %>% plot #Te saca solo la parte del dendograma del clustering.

library(RColorBrewer)
RColorBrewer::display.brewer.all()
rainbowcols <- rainbow(5, s = 0.9)
colores<-colorRampPalette(rainbowcols)(23)
names(colores)<-unique(MatrizMetabolismo$TISSUE)

#Le esta asociando a cada tejido un color diferente
col2<-as.data.frame(colores)
col2$TISSUE<-rownames(col2)
col3<-as.data.frame(as.character(MatrizMetabolismo[,-c(2:ncol(MatrizMetabolismo))]))
colnames(col3)<-"TISSUE"
col4<-merge(col3,col2, by="TISSUE")

png("clustering_feo.png", width = 7, height = 7, units = "in", res = 110)
dend2%>%
  set("labels_colors", "white") %>% 
  set("labels_cex", 0.01) %>%
  plot
colored_bars(colors = col4$colores, dend = dend2%>%set("labels_colors", "white"), sort_by_labels_order = T)
dev.off()

###########################################################################################################################################
################################### K-means (k=23, que son los tejidos)####################################################################
###########################################################################################################################################

#km.out <- kmeans(MatrizMetabolismo[,-c(1)] , cent=23 ) Si los tejidos estan en ingles
km.out <- kmeans(MatrizMetalbolFiltradaEspanol[,-c(1,2,3)] , cent=23 )
dim(km.out)
names (km.out)
#MatrizMetabolismo_clustered <- data.frame(MatrizMetabolismo, cluster=factor(km.out$cluster))
MatrizMetabolismo_clustered <- data.frame(MatrizMetalbolFiltradaEspanol, cluster=factor(km.out$cluster))

library(dplyr)
MatrizCuentas <- data.frame(MatrizMetabolismo_clustered$TISSUE, MatrizMetabolismo_clustered$cluster)
listado_tej <- unique(MatrizMetabolismo_clustered$TISSUE)
numero_tejidos<- c(1:length(listado_tej))
colnames(MatrizCuentas)<-c('Tissue','cluster')

m<-data.frame(matrix(ncol = length(numero_tejidos), nrow = length(listado_tej)))
colnames(m)<-as.character((1:length(listado_tej)))
rownames(m)<-sort(listado_tej, decreasing = FALSE)
for (tejido in listado_tej){
  print(tejido)
  z1 <-MatrizCuentas %>% filter (Tissue==tejido) #Agrupo por tejido
  for (numero in numero_tejidos){
    numeroClasificacion<-colSums(z1==numero) %>% .[-c(1)]
    m[tejido, numero]<-numeroClasificacion
  }
}
tablaFreqTej <- data.frame(table(MatrizMetabolismo_clustered$TISSUE))
v <- c(subset(tablaFreqTej, Freq > 5)) #Listado de los tejidos que con mas de 5 lineas por tejido.
v <- data.frame(subset(tablaFreqTej, Freq > 5))

ListaTejidoFreq<-as.data.frame(v)
colnames(ListaTejidoFreq) <- c("TISSUE", "Nº_lineas")
nº_lineas<-ListaTejidoFreq[,'Nº_lineas']
Mkmeans <-cbind(nº_lineas, m)

library(ggplot2)

Tejidos2 <- rownames(Mkmeans)
Tejidos2[20]<-"TEJIDO_HEMATOPOYETICO\nY_LINFOIDE"
Tejidos2[6]<-"GANGLIOS\nAUTONOMICOS"
Tejidos2[18]<-"SISTEMA_NERVIOSO\nCENTRAL"
Tejidos2[22]<-"TRACTO_AERODIGESTIVO\nSUPERIOR"
#### Ejecutar si los tejidos estan en ingles #######
#Tejidos2 <- rownames(Mkmeans)
#Tejidos2[8]<-"HAEMATOPOIETIC\nAND_LYMPHOID_TISSUE"
#Tejidos2[1]<-"AUTONOMIC\nGANGLIA"
#Tejidos2[5]<-"CENTRAL_NERVOUS\nSYSTEM"
#Tejidos2[22]<-"UPPER\nAERODIGESTIVE_TRACT"
rainbowcols <- rainbow(5, s = 0.9)
colores2<-colorRampPalette(rainbowcols)(24)
names(colores2)<-unique(Tejidos2)




png("Distribucion_tejidos_Metabolismo.png", width = 20, height = 10, units = "in", res = 110)
ggplot(Mkmeans, aes(Tejidos2, (nº_lineas), fill = Tejidos2))+
  geom_bar(stat = "identity") + 
  ylab("Número de lineas en cada tejido")+
  scale_fill_manual (values=c(colores))+
  theme(legend.position="none",
        axis.title.x = element_blank(), legend.title=element_blank(), axis.title.y = element_text(size=22, face="bold"),
        legend.background = element_rect(fill="lightgrey", size=0.5, linetype="solid"),
        panel.background = element_rect(fill='darkgrey', colour='darkgrey'), axis.text.x = element_text(angle = 65, hjust = 1, face="bold", size=16),
        axis.line = element_line(color = "darkblue", size = 0.5, linetype = "solid"),
        axis.text.y = element_text(face="bold", size=14))

dev.off()

#Representacion del k-means- En un eje, el dendograma de como se ajustan los cluster segun k-means, y el otro eje los tejidos. Lo ideal si cada clustr fuera una tejido es la diagonal de un color y el resto de otro.
#heatmap
install.packages("gplots")
library(gplots)
library(RColorBrewer)
display.brewer.all()
col<- colorRampPalette(brewer.pal(9,"Blues"))(50)
col2<-c(col3[1] , col[2:50])
col3<-c("#fbfbfb")

# probando colores
png("Heatmap_cluster_Metabolicos.png", width = 30, height = 15, units = "in", res = 110)
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

col<- colorRampPalette(brewer.pal(9,"Blues"))(200)
col3<-c("#fbfbfb")
col2<-c(col3[1] , col[2:200])

m <- as.data.frame(Mkmeans/(Mkmeans[,1]))
setwd("/home/estrella/TFM/Datos/depmap_19Q1_relabeled/GenesMetabolicos_(Agrupado por tejido)/")

png("Heatmap_cluster_normalizado_metabol.png", width = 30, height = 15, units = "in", res = 110)
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
  scale_fill_manual (values=c(colores3))+
  labs (fill = "CLUSTERS")+
  theme(legend.position="right", legend.title = element_text(size=16, face="bold"),
        axis.title.x = element_blank(), axis.title.y = element_text(size=18, face="bold"),
        legend.background = element_rect(fill="lightgrey", size=0.8, linetype="solid"),
        panel.background = element_rect(fill='darkgrey', colour='darkgrey'), axis.text.x = element_text(angle = 65, hjust = 1, face="bold", size=16),
        axis.text.y=element_blank(), legend.text=element_text(size=14, face = 'bold'),
        axis.line = element_line(color = "darkblue", size = 0.5, linetype = "solid"))
dev.off()

save.image("./20190722_metabol.RData")


# haciendo boxplot de todos los genes en todos los tejidos ----
library(ggplot2)
ggplot(matrizDatosCrudos) + geom_boxplot(aes(y=matrizDatosCrudos$TSPAN6))

library(reshape2)
library(dplyr)
meltData <- (matrizDatosCrudos[,-c(1:2)]) %>% .[,1:10] %>% melt()
p <- ggplot(meltData, aes(factor(variable), value)) 
p <- p + geom_boxplot()

means <- matrizDatosCrudos[,-c(1:2)] %>% .[,1:10]  %>% colMeans() %>% as.data.frame() 
means$name <- rownames(means)
myorder <- means$name[order(-((means$.)))]

medians <- matrizDatosCrudos[,-c(1:2)] %>% .[,1:10]  %>% apply(., 2, FUN = median) %>% as.data.frame() 
medians$name <- rownames(medians)
myorder <- medians$name[order(-((medians$.)))]
ggplot(meltData, aes(factor(variable), value))+ geom_boxplot()  + scale_x_discrete(limits=myorder)

meltData <- (matrizDatosCrudos[,-c(1:2)]) %>% melt()
medians <- matrizDatosCrudos[,-c(1:2)] %>% apply(., 2, FUN = median) %>% as.data.frame() 
medians$name <- rownames(medians)
myorder <- medians$name[order(-((medians$.)))]
p<-ggplot(meltData, aes(factor(variable), value))+ geom_boxplot()  + scale_x_discrete(limits=myorder)
ggsave("boxplot_todosgenes.png",p)
ggsave("boxplot_todosgenes_largo.png",p,width = 100,height = 9,units = "cm", dpi = 320)
ggsave("boxplot_todosgenes_largo.pdf",p,width = 100,height = 9,units = "cm", device = cairo_pdf)


