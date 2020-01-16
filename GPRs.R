setwd("/home/estrella/TFM/Datos/depmap/GPR/")
library(data.table)
GP_crudo <- read.table (file='/home/estrella/TFM/Datos/Merge_GP_expression.csv', header = TRUE, sep='\t')
filas <- as.character(GP_crudo$reaction_id)
GP_crudo <- as.data.frame(lapply(GP_crudo[,-c(1)], as.numeric))
rownames(GP_crudo) <- filas

GP_crudo[GP_crudo==-1000] <- NA

######## MEDIA, MEDIANA, VARIANZA, SD Y CV de la matriz con los valores de la GPRs de expresion ############################################
#Calculo media de cada reaccion.
GP_crudo_t <-as.data.frame(t(GP_crudo))
GP_crudo_t<-GP_crudo_t[colSums(!is.na(GP_crudo_t)) > 0]

media_GPR<-as.data.frame(apply(GP_crudo_t,2,mean))
colnames(media_GPR)<- 'media_Rx'
 
#Calculo de la varianza de cada reaccion.
varianza_GPR <-as.data.frame(apply(GP_crudo_t,2,var))
colnames(varianza_GPR)<-'varianza_Rx'

#Calculo de la varianza de cada reaccion.
sd_GPR <-as.data.frame(apply(GP_crudo_t,2,sd))
colnames(sd_GPR)<-'sd_Rx'

#Uno media y varianzas y sd.
media_var_sd <- cbind(media_GPR, varianza_GPR, sd_GPR)
summary(media_var_sd)

#Calculo el coeficiente de variacion de cada reaccion.
CV_GPR<- as.data.frame(media_var_sd$sd_Rx/media_var_sd$media_Rx)
colnames(CV_GPR)<-"cv"
media_var_sd <- cbind(media_var_sd, CV_GPR)

#################################################################################################################################################################
########  PARA LA MATRIZ DE CERES ###############################################################################################################################
library(data.table)
GP_crudo_ceres <- read.table (file='/home/estrella/TFM/Datos/Merge_GP_ceres.csv', header = TRUE, sep='\t', stringsAsFactors = F)
filas <- as.character(GP_crudo_ceres$reaction_id)
GP_crudo_ceres <- as.data.frame(lapply(GP_crudo_ceres[,-c(1)], as.numeric))
rownames(GP_crudo_ceres) <- filas
GP_crudo_ceres[GP_crudo_ceres==-1000] <- NA


######## MEDIA, MEDIANA, VARIANZA, SD Y CV de la matriz con los valores de la GPRs de CERES ############################################
#Calculo media de cada reaccion.
GP_crudo_ceres_t <-as.data.frame(t(GP_crudo_ceres))
GP_crudo_ceres_t<- GP_crudo_ceres_t[colSums(!is.na(GP_crudo_ceres_t)) > 0]
media_GPR_ceres<-as.data.frame(apply(GP_crudo_ceres_t,2,mean))
colnames(media_GPR_ceres)<- 'media_Rx_CERES'

#Calculo de la varianza de cada reaccion.
varianza_GPR_ceres <-as.data.frame(apply(GP_crudo_ceres_t,2,var))
colnames(varianza_GPR_ceres)<-'varianza_Rx_CERES'

#Calculo de la varianza de cada reaccion.
sd_GPR_ceres <-as.data.frame(apply(GP_crudo_ceres_t,2,sd))
colnames(sd_GPR_ceres)<-'sd_Rx_CERES'

#Uno media y varianzas y sd y cv.
media_var_sd_CERES <- cbind(media_GPR_ceres, varianza_GPR_ceres, sd_GPR_ceres)
summary(media_var_sd)

#Calculo el coeficiente de variacion de cada reaccion.
CV_GPR_ceres <- as.data.frame(media_var_sd_CERES$sd_Rx_CERES/media_var_sd_CERES$media_Rx_CERES)
colnames(CV_GPR_ceres)<-"cv_CERES"
media_var_sd_CERES <- cbind(media_var_sd_CERES, CV_GPR_ceres)


t<- rownames(media_var_sd_CERES)
media_var_sd_CERES <-cbind(t,media_var_sd_CERES)
Matriz_medias_ExpYCeres<- merge(media_var_sd, media_var_sd_CERES, by="t", all.x = T)

############################################################################################################
###########################  GRAFICAS MEDIA  ###############################################################
############################################################################################################

#~~~~~~~~~~~~~~~~ Plot  varianzaVsMedia (expresion) de todas las Rx ~~~~~~~~~~~~~~~~~~~~~
library(ggplot2)
Varianza <-media_var_sd$varianza_Rx
Media<-media_var_sd$media_Rx
png("MediavsVarianzaGP.png", width = 13, height = 13, units = "in", res = 110)
ggplot(media_var_sd, aes(x=Media,y=Varianza)) +
  geom_point(colour = "chartreuse4")+
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.text=element_text(size=16), axis.title=element_text(size=28,face="bold"), plot.title=element_text(size=32, face="bold"))+
  labs(x='Media', y='Varianza')+
  ggtitle("Representación media vs varianza de los valores de las GPR de expresión")
# ylim(c(0,18))+
# xlim(c(0,18))
dev.off()

#~~~~~~~~~~~~~~~~ Plot histograma de RX según su media (expresion) ~~~~~~~~~~~~~~~~~~~~~~~~~~

library(ggplot2)
RX <- media_var_sd$media_Rx
h<-hist(media_var_sd$media_Rx, breaks = 50, plot = FALSE) #De esta forma agrupa los Rx de 200.
png("HistogramaPGR.png", width = 13, height = 13, units = "in", res = 110)
plot(h, xaxt = "n", xlab = "Reacciones", ylab = "Frecuencia", cex.lab=2, cex.main=2.5,
     main = "Histograma de las medias de las GPR sobre datos de expresión", col = "chartreuse4")
dev.off()


cuartiles<-quantile(media_var_sd$media_Rx)
png("HistogramaPGR_2.png", width = 13, height = 13, units = "in", res = 110)
ggplot(Matriz_medias_ExpYCeres,aes(x=media_Rx))+
  geom_histogram(binwidth=0.1, fill="chartreuse4",color="black")+
  geom_density(aes(y=0.1 * ..count..), colour="red")+
 theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.text=element_text(size=16), axis.title=element_text(size=28,face="bold"), plot.title=element_text(size=32, face="bold"))+
  geom_vline (aes(xintercept=cuartiles[2]), color="blue", linetype="dashed", size=0.8)+
  geom_vline (aes(xintercept=cuartiles[4]), color="blue", linetype="dashed", size=0.8)+
  geom_vline (aes(xintercept=cuartiles[3]), color="blue", linetype="dashed", size=0.8)+
  labs(x='Expresión', y='Densidad')
dev.off()


# #library("qplot")
# Matriz_medias_ExpYCeres_2<- Matriz_medias_ExpYCeres[order(Matriz_medias_ExpYCeres$media_Rx),]
# #ggplot(Matriz_medias_ExpYCeres_2, aes(t, media_Rx)) + geom_bar(stat = "identity")
# ggplot(Matriz_medias_ExpYCeres_2, aes(x=media_Rx)) + geom_density()

#~~~~~~~~~~~~~~~~ Plot media vs cv ( expresion) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
png("MediavsCoeficienteVariacion.png", width = 13, height = 13, units = "in", res = 110)
ggplot(media_var_sd, aes(x=Media,y=cv)) +
  geom_point(colour = "goldenrod")+
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.text=element_text(size=16), axis.title=element_text(size=28,face="bold"), plot.title=element_text(size=32, face="bold"))+
  labs(x='Media', y='Coeficiente Variacion')+
  ggtitle("Representación media vs Coeficiente de variación")
# geom_text(aes(x = media_var_metabol[(media_var_metabol$cv>10),]$media_gen, y = media_var_metabol[(media_var_metabol$cv>10),]$cv, 
#               colour = media_var_metabol[(media_var_metabol$cv>10),]$cv),label=media_var_metabol[(media_var_metabol$cv>10),]$`media_var_metabol$cv`, size=3) +
# geom_text_repel(aes(x = media_var_metabol[(media_var_metabol$cv>10),]$media_gen, y = media_var_metabol[(media_var_metabol$cv>10),]$cv, 
#                     colour = media_var_metabol[(media_var_metabol$media_gen>10),]$cv) ,
#                 label=media_var_metabol[(media_var_metabol$media_gen>10),]$`media_var_metabol$cv`, size=5.5, face="bold", color='black')

dev.off()
save.image("./20190801_GPR.RData")

#~~~~~~~~~~~~~~~~ Plot histograma de RX según su media (esenciabilidad) ~~~~~~~~~~~~~~~~~~~~~~~~~~

RX_ceres <- media_var_sd_CERES$media_Rx_CERES
h<-hist(media_var_sd_CERES$media_Rx_CERES, breaks = 50, plot = FALSE) #De esta forma agrupa los Rx de 200.
png("HistogramaPGR_ceres.png", width = 13, height = 13, units = "in", res = 110)
plot(h, xaxt = "n", xlab = "Reacciones", ylab = "Frecuencia", cex.lab=2, cex.main=2.5,
     main = "Histograma de las medias de las GPR sobre datos de esenciabilidad", col = "chartreuse4")
dev.off()

#~~~~~~~~~~~~~~~~ Plot mediaVSVarianza (esenciabilidad) ~~~~~~~~~~~~~~~~~~~~~~~~~~

Varianza <-media_var_sd_CERES$varianza_Rx_CERES
Media<-media_var_sd_CERES$media_Rx_CERES
png("MediavsVarianzaCERES.png", width = 13, height = 13, units = "in", res = 110)
ggplot(media_var_sd_CERES, aes(x=Media,y=Varianza)) +
  geom_point(colour = "chartreuse4")+
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.text=element_text(size=16), axis.title=element_text(size=28,face="bold"), plot.title=element_text(size=32, face="bold"))+
  labs(x='Media', y='Varianza')+
  ggtitle("Representación media vs varianza de los valores de las GPR de CERES")
# ylim(c(0,18))+
# xlim(c(0,18))
dev.off()

#~~~~~~~~~~~~~~~~ Plot ExpresionVsCERES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Oredeno la matriz de expresion
GP_ceres <-GP_crudo_ceres_t[,-c(1:3)]
GP_expr<- GP_crudo_t[,-c(1:3)]
#Ordeno las columna por orden alfabetico.
GP_ceres<-GP_ceres[ , order(names(GP_ceres))]
GP_expr<-GP_expr[ , order(names(GP_expr))]
#Compruebo que el numero de columnas son las mismas.
# z<-data.frame(colnames(GP_expr))# Genes ordenados que salen de la matriz expresion con genes metabolicos.
# z2<-data.frame(colnames(GP_ceres)) #Genes ordenados que salen de la matriz de ceres con genes metabolicos.
#Quito las lineas en expresion que no hay en CERES.
GP_expr_igual<-GP_expr[rownames(GP_expr) %in% as.character(rownames(GP_ceres)),]
GP_ceres_igual<-GP_ceres[rownames(GP_ceres) %in% as.character(rownames(GP_expr_igual)),]
#Para ordenar de forma alfabetica las filas de las dos matrices.
GP_expr_igual<-GP_expr_igual[order(rownames(GP_expr_igual)) ,]
GP_ceres_igual<-GP_ceres_igual[order(rownames(GP_ceres_igual)) ,]

library(reshape2)
Expresion<-melt(as.matrix(GP_expr_igual))
colnames(Expresion)<-c("linea", "RX", "expresion")
Esencial<-melt(as.matrix(GP_ceres_igual))
colnames(Esencial)<-c("linea", "RX", "ceres")
Expr_ceres<-data.frame()
Expr_ceres<-cbind(Expresion, Esencial[,3])
colnames(Expr_ceres) <-c("linea", "RX", "expresion", "ceres")


library(ggplot2)
library(ggrepel)
png("ExpresionVsCeresGP.png", width = 13, height = 13, units = "in", res = 110)
ggplot(Expr_ceres, aes(x=expresion , y=ceres)) +
  geom_point(colour = "darkblue")+
  #ggtitle("Distribución expresion vs esenciabilidad")+
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.text=element_text(size=18), axis.title=element_text(size=30,face="bold"), plot.title=element_text(size=34, face="bold"), 
        axis.text.x = element_text(size=20, face="bold"),
        axis.text.y = element_text(size=20, face="bold"))+ 
  xlab("log2 TPM")+
  ylab("CERES score")+
  geom_hline(yintercept = -1, color="red", linetype="dashed", size=1.5)+
  scale_y_continuous(breaks = c(round(min(Expr_ceres$ceres),digits = 0), -1,0,1, round(max(Expr_ceres$ceres),digits = 0)))
dev.off()

library(dplyr)
Expr_ceres %>% filter(ceres< -1 & expresion == 0) %>% arrange(expresion)

##########################################################################################################################################
#################     PCA   ##############################################################################################################
##########################################################################################################################################
#~~~~~~~~~~~ EXPRESION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(dplyr)
# DA UN PROBLEMA: cannot rescale a constant/zero column to unit variance. Quito por eso los NAs.
#Se supone que es pq alguna columna es una constante o es cero y por tanto su varianza es cero.
#La inteligencia colectiva propone, quitar esos ceros (que aqui son los NA) y luego hacer un PCA.
pca.out <- prcomp(GP_crudo_t, scale=T )
PoV <- (pca.out$sdev^2/sum(pca.out$sdev^2))*100 #Proporcion de varianza explicada por cada CP1
head(PoV)

library(dplyr)
library(reshape2)
CCLE<-as.data.frame(rownames(GP_crudo_t)) #Le añadimos las columna CCLE(linea_tejido), para crear dos: ID y TISSUE. 
div<-colsplit(CCLE$`rownames(GP_crudo_t)`, "_", c("ID","TISSUE"))
colnames(CCLE) <-"CCLE"
GP_crudo_t <- cbind(CCLE, div, GP_crudo_t)

pca.out2 <- as.data.frame(pca.out$x)
pca.out2 <- cbind(GP_crudo_t$TISSUE, pca.out2)
colnames(pca.out2)[1]<-"TISSUE"
pca.out2$TISSUE<-as.character(pca.out2$TISSUE)
pca.out2 <- cbind(GP_crudo_t$CCLE, pca.out2)
colnames(pca.out2)[1]<-"CCLE"

#~~~~~~~~  CERES  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(dplyr)
pca.out.ceres <- prcomp(GP_crudo_ceres_t, scale=T )
PoV.ceres <- (pca.out.ceres$sdev^2/sum(pca.out.ceres$sdev^2))*100 #Proporcion de varianza explicada por cada CP1
head(PoV.ceres)

library(dplyr)
library(reshape2)
CCLE<-as.data.frame(rownames(GP_crudo_ceres_t)) #Le añadimos las columna CCLE(linea_tejido), para crear dos: ID y TISSUE. 
div<-colsplit(CCLE$`rownames(GP_crudo_ceres_t)`, "_", c("ID","TISSUE"))
colnames(CCLE) <-"CCLE"
GP_crudo_ceres_t <- cbind(CCLE, div, GP_crudo_ceres_t)

pca.out2.ceres <- as.data.frame(pca.out.ceres$x)
pca.out2.ceres <- cbind(GP_crudo_ceres_t$TISSUE, pca.out2.ceres)
colnames(pca.out2.ceres)[1]<-"TISSUE"
pca.out2.ceres$TISSUE<-as.character(pca.out2.ceres$TISSUE)
pca.out2.ceres <- cbind(GP_crudo_ceres_t$CCLE, pca.out2.ceres)
colnames(pca.out2.ceres)[1]<-"CCLE"

#½½½½½½½½½½½½½½½½½½½½½½½½½½½½ ESTA REPRESENTACIÓN ES CON TODOS LOS TEJIDOS, SIN QUITAR AQUELLOS CON MENOS DE CINCO LINEAS CELULARES  ½½½½½½½½½½½½½½½½½½½
#~~~~~~~~~~~~~~~~~~~~~~~ EXPRESION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(ggplot2)
library(ggrepel)
TEJIDOS <- pca.out2$TISSUE
listaTej <- data.frame(table(GP_crudo_t$TISSUE))
#rainbowcols <- rainbow(5, s = 0.9)
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
colores<-colorRampPalette(color)(31)
names(colores)<-unique(listaTej$Var1)
#percentage <- round(pca.out$sdev / sum(pca.out$sdev) * 100, 2) #Calula los porcentajes de cada CPs.
percentage <- round(PoV, 2)
png("PCA_PG_expresion_TODOStej.png", width = 17, height = 15, units = "in", res = 110)
ggplot() +
  geom_point(data = pca.out2, aes(x = pca.out2$PC1, y = pca.out2$PC2, colour = TEJIDOS), size=3) +
  theme(legend.position="top", legend.key.size = unit(1, "cm"), #legend.title = element_text(size=6, face="bold"),
        axis.title.x = element_text(size=24, face="bold"), legend.title=element_blank(), axis.title.y = element_text(size=24, face="bold"),
        legend.background = element_rect(fill="gray90", size=4.5, linetype="solid"),
        panel.background = element_rect(fill='gray22', colour='darkgrey'), legend.text=element_text(size=14, face = 'bold')) +
  guides(colour = guide_legend(override.aes = list(size=5)))+
  # geom_text(aes(x = pca.out2[(pca.out2$PC1>100),]$PC1, y = pca.out2[(pca.out2$PC1>100),]$PC2, colour = pca.out2[(pca.out2$PC1>100),]$TISSUE),label=pca.out2[(pca.out2$PC1>100),]$`GP_crudo_t$ID`, size=3) +
  # geom_text_repel(aes(x = pca.out2[(pca.out2$PC1>100),]$PC1, y = pca.out2[(pca.out2$PC1>100),]$PC2,
  #     colour = pca.out2[(pca.out2$PC1>100),]$TISSUE),
  # label=pca.out2[(pca.out2$PC1>100),]$`GP_crudo_t$CCLE`, size=5.5, face="bold", color='black') +
  # theme_minimal() +
  scale_colour_manual(values=colores) +
  xlab(paste("PC1 (", paste(as.character(percentage[1]), "%)", sep="") )) +
  ylab(paste("PC2 (", paste(as.character(percentage[2]), "%)", sep="") ))
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~ CERES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(ggplot2)
library(ggrepel)
TEJIDOS.ceres <- pca.out2.ceres$TISSUE
listaTej.ceres <- data.frame(table(GP_crudo_ceres_t$TISSUE))
#rainbowcols <- rainbow(5, s = 0.9)
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
colores<-colorRampPalette(color)(25)
names(colores)<-unique(listaTej.ceres$Var1)
#percentage <- round(pca.out$sdev / sum(pca.out$sdev) * 100, 2) #Calula los porcentajes de cada CPs.
percentage.ceres <- round(PoV.ceres, 2)
png("PCA_PG_CERES_TODOStej.png", width = 17, height = 15, units = "in", res = 110)
ggplot() +
  geom_point(data = pca.out2.ceres, aes(x = pca.out2.ceres$PC1, y = pca.out2.ceres$PC2, colour = TEJIDOS.ceres), size=5) +
  theme(legend.position="top", legend.key.size = unit(1, "cm"), #legend.title = element_text(size=6, face="bold"),
        axis.title.x = element_text(size=24, face="bold"), legend.title=element_blank(), axis.title.y = element_text(size=24, face="bold"),
        legend.background = element_rect(fill="gray90", size=5.5, linetype="solid"),
        panel.background = element_rect(fill='gray22', colour='darkgrey'), legend.text=element_text(size=13, face = 'bold')) +
  guides(colour = guide_legend(override.aes = list(size=5)))+
  # geom_text(aes(x = pca.out2.ceres[(pca.out2.ceres$PC1>100),]$PC1, y = pca.out2.ceres[(pca.out2.ceres$PC1>100),]$PC2, colour = pca.out2.ceres[(pca.out2.ceres$PC1>100),]$TISSUE),label=pca.out2.ceres[(pca.out2.ceres$PC1>100),]$`GP_crudo_ceres_t$ID`, size=3) +
  # geom_text_repel(aes(x = pca.out2.ceres[(pca.out2.ceres$PC1>100),]$PC1, y = pca.out2.ceres[(pca.out2.ceres$PC1>100),]$PC2,
  #     colour = pca.out2.ceres[(pca.out2.ceres$PC1>100),]$TISSUE),
  # label=pca.out2.ceres[(pca.out2.ceres$PC1>100),]$`GP_crudo_ceres_t$CCLE`, size=5.5, face="bold", color='black') +
  # theme_minimal() +
  scale_colour_manual(values=colores) +
  xlab(paste("PC1 (", paste(as.character(percentage.ceres[1]), "%)", sep="") )) +
  ylab(paste("PC2 (", paste(as.character(percentage.ceres[2]), "%)", sep="") ))
dev.off()

save.image("./20190802_GPR.RData")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%% QUITO LOS TEJIDOS CON MENOS DE CINCO LINEAS CELULAREs   %%%%%%%%%%%%%%%%%%%%%%%%%%%
#################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~Expresion  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
v <- data.frame(subset(listaTej, Freq > 5))#Listado de los tejidos que con mas de 5 lineas por tejido.
GP_crudo_t_filtrada <-GP_crudo_t[GP_crudo_t$TISSUE%in%as.character(v$Var1),]
pca.out.filt <- prcomp(GP_crudo_t_filtrada[, -c(1:3)], scale=T ) #Quito los NA de la matriz filtrada
PoV.filt <- (pca.out.filt$sdev^2/sum(pca.out.filt$sdev^2))*100 #Proporcion de varianza explicada por cada CP1
head(PoV)

pca.out2.filt <- as.data.frame(pca.out.filt$x)
pca.out2.filt <- cbind(GP_crudo_t_filtrada$TISSUE, pca.out2.filt)
colnames(pca.out2.filt)[1]<-"TISSUE"
pca.out2.filt$TISSUE<-as.character(pca.out2.filt$TISSUE)
pca.out2.filt <- cbind(GP_crudo_t_filtrada$CCLE, pca.out2.filt)
colnames(pca.out2.filt)[1]<-"CCLE"

TEJIDOS <- pca.out2.filt$TISSUE
listaTej.filt <- data.frame(table(GP_crudo_t_filtrada$TISSUE))
# rainbowcols <- rainbow(5, s = 0.9)
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
colores2<-colorRampPalette(color)(24)
names(colores2)<-unique(listaTej.filt$Var1)
#percentage <- round(pca.out$sdev / sum(pca.out$sdev) * 100, 2) #Calula los porcentajes de cada CPs.
percentage2 <- round(PoV.filt, 2)
png("PCA_PG_expresion_filtradatej.png", width = 17, height = 15, units = "in", res = 110)
ggplot() +
  geom_point(data = pca.out2.filt, aes(x = pca.out2.filt$PC1, y= pca.out2.filt$PC2, colour = TEJIDOS), size=3) +
  theme(legend.position="top", legend.key.size = unit(1, "cm"), #legend.title = element_text(size=6, face="bold"),
        axis.title.x = element_text(size=24, face="bold"), legend.title=element_blank(), axis.title.y = element_text(size=24, face="bold"),
        legend.background = element_rect(fill="gray90", size=4.5, linetype="solid"),
        panel.background = element_rect(fill='gray22', colour='darkgrey'), legend.text=element_text(size=14, face = 'bold')) +
  guides(colour = guide_legend(override.aes = list(size=5)))+
  # geom_text(aes(x = pca.out2.filt[(pca.out2.filt$PC1>100),]$PC1, y = pca.out2.filt[(pca.out2$PC1>100),]$PC2, colour = pca.out2.filt[(pca.out2.filt$PC1>100),]$TISSUE),label=pca.out2.filt[(pca.out2.filt$PC1>100),]$`GP_crudo_t_filtrada$ID`, size=3) +
  # geom_text_repel(aes(x = pca.out2.filt[(pca.out2.filt$PC1>100),]$PC1, y = pca.out2.filt[(pca.out2.filt$PC1>100),]$PC2,
  #     colour = pca.out2.filt[(pca.out2.filt$PC1>100),]$TISSUE),
  # label=pca.out2.filt[(pca.out2.filt$PC1>100),]$`GP_crudo_t_filtrada$CCLE`, size=5.5, face="bold", color='black') +
  # theme_minimal() +
  scale_colour_manual(values=colores2) +
  xlab(paste("PC1 (", paste(as.character(percentage2[1]), "%)", sep="") )) +
  ylab(paste("PC2 (", paste(as.character(percentage2[2]), "%)", sep="") ))
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~  CERES  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
v2 <- data.frame(subset(listaTej.ceres, Freq > 5)) #Listado de los tejidos que con mas de 5 lineas por tejido.
GP_crudo_ceres_t_filtrada <-GP_crudo_ceres_t[GP_crudo_ceres_t$TISSUE%in%as.character(v2$Var1),]
pca.out.ceres.filt <- prcomp(GP_crudo_ceres_t_filtrada[, -c(1:3)], scale=T ) #Quito los NA de la matriz filtrada
PoV.ceres.filt <- (pca.out.ceres.filt$sdev^2/sum(pca.out.ceres.filt$sdev^2))*100 #Proporcion de varianza explicada por cada CP1
head(PoV.ceres.filt)

pca.out2.ceres.filt <- as.data.frame(pca.out.ceres.filt$x)
pca.out2.ceres.filt <- cbind(GP_crudo_ceres_t_filtrada$TISSUE, pca.out2.ceres.filt)
colnames(pca.out2.ceres.filt)[1]<-"TISSUE"
pca.out2.ceres.filt$TISSUE<-as.character(pca.out2.ceres.filt$TISSUE)
pca.out2.ceres.filt <- cbind(GP_crudo_ceres_t_filtrada$CCLE, pca.out2.ceres.filt)
colnames(pca.out2.ceres.filt)[1]<-"CCLE"


TEJIDOS.ceres <- pca.out2.ceres.filt$TISSUE
listaTej.filt.ceres <- data.frame(table(GP_crudo_ceres_t_filtrada$TISSUE))
#rainbowcols <- rainbow(5, s = 0.9)
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
colores2<-colorRampPalette(color)(19)
names(colores2)<-unique(listaTej.filt.ceres$Var1)
#percentage <- round(pca.out$sdev / sum(pca.out$sdev) * 100, 2) #Calula los porcentajes de cada CPs.
percentage.ceres2 <- round(PoV.ceres.filt, 2)
png("PCA_PG_CERES_filtradatej.png", width = 17, height = 15, units = "in", res = 110)
ggplot() +
  geom_point(data = pca.out2.ceres.filt, aes(x = pca.out2.ceres.filt$PC1, y = pca.out2.ceres.filt$PC2, colour = TEJIDOS.ceres), size=5) +
  theme(legend.position="top", legend.key.size = unit(1, "cm"), #legend.title = element_text(size=6, face="bold"),
        axis.title.x = element_text(size=28, face="bold"), legend.title=element_blank(), axis.title.y = element_text(size=28, face="bold"),
        legend.background = element_rect(fill="gray90", size=5, linetype="solid"),
        panel.background = element_rect(fill='gray22', colour='darkgrey'), legend.text=element_text(size=14, face = 'bold')) +
  guides(colour = guide_legend(override.aes = list(size=5)))+
  # geom_text(aes(x = pca.out2.ceres.filt[(pca.out2.ceres.filt$PC1>100),]$PC1, y = pca.out2.ceres.filt[(pca.out2.ceres.filt$PC1>100),]$PC2, colour = pca.out2.ceres.filt[(pca.out2.ceres.filt$PC1>100),]$TISSUE),label=pca.out2.ceres.filt[(pca.out2.ceres.filt$PC1>100),]$`pca.out2.ceres.filt$ID`, size=3) +
  # geom_text_repel(aes(x = pca.out2.ceres.filt[(pca.out2.ceres.filt$PC1>100),]$PC1, y = pca.out2.ceres.filt[(pca.out2.ceres.filt$PC1>100),]$PC2,
  #     colour = pca.out2.ceres.filt[(pca.out2.ceres.filt$PC1>100),]$TISSUE),
  # label=pca.out2.ceres.filt[(pca.out2.ceres.filt$PC1>100),]$`GP_crudo_ceres_t_filtrada$CCLE`, size=5.5, face="bold", color='black') +
  # theme_minimal() +
  scale_colour_manual(values=colores2) +
  xlab(paste("PC1 (", paste(as.character(percentage2[1]), "%)", sep="") )) +
  ylab(paste("PC2 (", paste(as.character(percentage2[2]), "%)", sep="") ))
dev.off()

########################################################################################################################
#~~~~~~~~~~~~~~~~~~ Quitando hemato ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~ Expresion  (con ceres no lo voy a hacer, pq no creo que se clarifique nada)~~~~~~~~~~~~~~~~~~~~~~~
GP_crudo_t_filtradaNH <- GP_crudo_t_filtrada%>%filter(!(TISSUE == "HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"))
pca.out.NH <- prcomp(GP_crudo_t_filtradaNH[, -c(1:3)], scale=T ) #Quito los NA de la matriz filtrada
PoV.NH <- (pca.out.NH$sdev^2/sum(pca.out.NH$sdev^2))*100 #Proporcion de varianza explicada por cada CP1
head(PoV.NH)

pca.out2.NH <- as.data.frame(pca.out.NH$x)
pca.out2.NH <- cbind(GP_crudo_t_filtradaNH$TISSUE, pca.out2.NH)
colnames(pca.out2.NH)[1]<-"TISSUE"
pca.out2.NH$TISSUE<-as.character(pca.out2.NH$TISSUE)
pca.out2.NH <- cbind(GP_crudo_t_filtradaNH$CCLE, pca.out2.NH)
colnames(pca.out2.NH)[1]<-"CCLE"

TEJIDOS.NH <- pca.out2.NH$TISSUE
listaTej.NH <- data.frame(table(GP_crudo_t_filtradaNH$TISSUE))
# rainbowcols <- rainbow(5, s = 0.9)
# colores2<-colorRampPalette(rainbowcols)(23)
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
colores2<-colorRampPalette(color)(23)
names(colores2)<-unique(listaTej.NH$Var1)
#percentage <- round(pca.out$sdev / sum(pca.out$sdev) * 100, 2) #Calula los porcentajes de cada CPs.
percentage.NH <- round(PoV.NH, 2)
png("PCA_PG_sinHemato.png", width = 17, height = 15, units = "in", res = 110)
ggplot() +
  geom_point(data = pca.out2.NH, aes(x = pca.out2.NH$PC1, y = pca.out2.NH$PC2, colour = TEJIDOS.NH), size=5) +
  theme(legend.position="top", legend.key.size = unit(1, "cm"), #legend.title = element_text(size=6, face="bold"),
        axis.title.x = element_text(size=28, face="bold"), legend.title=element_blank(), axis.title.y = element_text(size=28, face="bold"),
        legend.background = element_rect(fill="gray90", size=12, linetype="solid"),
        panel.background = element_rect(fill='gray22', colour='darkgrey'), legend.text=element_text(size=14, face = 'bold', color = "black")) +
  guides(colour = guide_legend(override.aes = list(size=5)))+
  # geom_text(aes(x = pca.out2.NH[(pca.out2.NH$PC1>100),]$PC1, y = pca.out2.NH[(pca.out2.NH$PC1>100),]$PC2, colour = pca.out2.NH[(pca.out2.NH$PC1>100),]$TISSUE),label=pca.out2.NH[(pca.out2.NH$PC1>100),]$`pca.out2.NH$ID`, size=3) +
  # geom_text_repel(aes(x = pca.out2.NH[(pca.out2.NH$PC1>100),]$PC1, y = pca.out2.NH[(pca.out2.NH$PC1>100),]$PC2,
  #     colour = pca.out2.NH[(pca.out2.NH$PC1>100),]$TISSUE),
  # label=pca.out2.NH[(pca.out2.NH$PC1>100),]$`GP_crudo_t_filtradaNH$CCLE`, size=5.5, face="bold", color='black') +
  # theme_minimal() +
  scale_colour_manual(values=colores2) +
  xlab(paste("PC1 (", paste(as.character(percentage.NH[1]), "%)", sep="") )) +
  ylab(paste("PC2 (", paste(as.character(percentage.NH[2]), "%)", sep="") ))
dev.off()

save.image("./20190802_GPR.RData")
############################################################################################################################################
################################   k-means (K=24, nº tejidos) ##############################################################################
############################################################################################################################################
km.out <- kmeans(GP_crudo_t_filtrada[,-c(1,2,3)] , cent=24)
names (km.out)

GPFiltrada_clustered <- data.frame(GP_crudo_t_filtrada, cluster=factor(km.out$cluster))
colnames(listaTej.filt) <- c("TISSUE", "Nº_lineas")
listado_tej <-listaTej.filt$TISSUE
library(dplyr)
MatrizCuentasGP <- data.frame(GPFiltrada_clustered$TISSUE, GPFiltrada_clustered$cluster)
numero_tejidos<- c(1:length(listado_tej))
colnames(MatrizCuentasGP)<-c('Tissue','cluster')

m<-data.frame(matrix(ncol = length(numero_tejidos), nrow = length(listado_tej)))
colnames(m)<-as.character((1:length(listado_tej)))
rownames(m)<-listado_tej
for (tejido in listado_tej){
  print(tejido)
  z1 <-MatrizCuentasGP %>% filter(Tissue==tejido) #Agrupo por tejido
  for (numero in numero_tejidos){
    numeroClasificacion<-colSums(z1==numero) %>% .[-c(1)]
    m[tejido, numero]<-numeroClasificacion
  }
}


nº_lineas<-listaTej.filt[,'Nº_lineas']
Mkmeans <-cbind(nº_lineas, m)

library(gplots)
col<- colorRampPalette(brewer.pal(9,"Blues"))(200)
col2<-c(col3[1] , col[2:200])
col3<-c("#fbfbfb")
png("Heatmap_cluster_GP.png", width = 30, height = 15, units = "in", res = 110)
heatmap.2(as.matrix(Mkmeans[,-c(1)]),
          trace="none",
          #breaks = seq(0,100,2),
          col=col2,
          cexRow=2,
          cexCol = 2,
          margins = c(10, 40),
          keysize = 1.2,
          srtCol=0
)
dev.off()
m <- as.data.frame(Mkmeans/(Mkmeans[,1])) #Para normalizar por el numero de lineas.
png("Heatmap_cluster_GP_normalizado.png", width = 30, height = 15, units = "in", res = 110)
heatmap.2(as.matrix(m[,-c(1)]),
          trace="none", 
          #breaks = seq(0,150,0.1),
          # col=col2,
          col=col2,
          cexRow=2,
          cexCol = 2,
          margins = c(10, 40),
          keysize = 1.2,
          srtCol=0
          )
dev.off()
save.image("./20190802_GPR.RData")

########################################################################################################################################
########################### GRUPOS SEGUN LA REPRESENTACION DE MEDIA VS VARIANZA #########################################################
########################################################################################################################################
library(ggplot2)
library(ggExtra)
#~~~~~~~~~~~~~~~~~~~~~ Grafica con los boxplot  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
IQR <- IQR(media_var_sd$varianza_Rx)
upper_whisker = (quantile(Matriz_medias_ExpYCeres$varianza_Rx)[4]+ 1.5 * IQR)
Varianza <-media_var_sd$varianza_Rx
Media<-media_var_sd$media_Rx
png("MediavsVarianzaGP_boxplot.png", width = 13, height = 13, units = "in", res = 110)
scatter<- ggplot(media_var_sd, aes(x=Media,y=Varianza)) +
  geom_point(colour = "chartreuse4")+
  geom_vline (aes(xintercept=quantile(media_var_sd$media_Rx)[2]), color="blue", linetype="dashed", size=0.8)+
  geom_vline (aes(xintercept=quantile(media_var_sd$media_Rx)[3]), color="blue", linetype="dashed", size=0.8)+
  geom_vline (aes(xintercept=quantile(media_var_sd$media_Rx)[4]), color="blue", linetype="dashed", size=0.8)+
  geom_hline (aes(yintercept=quantile(media_var_sd$varianza_Rx)[4]), color="blue", linetype="dashed", size=0.8)+
  geom_hline (aes(yintercept=upper_whisker), color="blue", linetype="dashed", size=0.8)+
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.text=element_text(size=16), axis.title=element_text(size=28,face="bold"), plot.title=element_text(size=32, face="bold"))+
  labs(x='Media', y='Varianza')
ggMarginal(scatter, type="boxplot", size=10, colour = 'black', fill = 'darkgrey', xparams = list(size=0.5), yparams = list(size=0.5))
dev.off()

#~~~~~~~~~~~~~~~~ Representacion conjunta mediaVsvarianza genes metabol y GPRs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
media_var_metabol<-read.table("media_var_metabol.csv", header = T, sep = "\t")
Media2<- media_var_metabol$media_gen
Varianza2<- media_var_metabol$varianza_gen
png("MediavsVarianza_GenesMetabol_GPR.png", width = 13, height = 13, units = "in", res = 110)
ggplot(media_var_sd, aes(x=Media,y=Varianza)) +
  geom_point(data=media_var_metabol, aes(x=Media2,y=Varianza2), colour = "darkred")+
  geom_point(colour = "chartreuse4")+
  geom_vline (aes(xintercept=quantile(media_var_sd$media_Rx)[2]), color="blue", linetype="dashed", size=0.8)+
  geom_vline (aes(xintercept=quantile(media_var_sd$media_Rx)[3]), color="blue", linetype="dashed", size=0.8)+
  geom_vline (aes(xintercept=quantile(media_var_sd$media_Rx)[4]), color="blue", linetype="dashed", size=0.8)+
  geom_hline (aes(yintercept=quantile(media_var_sd$varianza_Rx)[4]), color="blue", linetype="dashed", size=0.8)+
  geom_hline (aes(yintercept=upper_whisker), color="blue", linetype="dashed", size=0.8)+
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.text=element_text(size=16), axis.title=element_text(size=28,face="bold"), plot.title=element_text(size=32, face="bold"))+
  labs(x='Media', y='Varianza')
dev.off()

#save.image("./20190802_GPR.RData")
#~~~~~~~~~~~~~~~~~~~~~ Grupos  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(dplyr)
t<- rownames(media_var_sd)
media_var_sd <-cbind(t,media_var_sd)
#GRUPO 1:
grupo1 <-data.frame(Matriz_medias_ExpYCeres%>%filter(media_Rx<quantile(Matriz_medias_ExpYCeres$media_Rx)[2] & 
                                                       varianza_Rx<=quantile(Matriz_medias_ExpYCeres$varianza_Rx)[4]))
write.table(grupo1, "grupo1.txt", quote = F, sep = "\t", row.names = F, col.names = T)

#GRUPO 2:
grupo2<-data.frame(Matriz_medias_ExpYCeres%>%filter(media_Rx>quantile(Matriz_medias_ExpYCeres$media_Rx)[2]&media_Rx<quantile(Matriz_medias_ExpYCeres$media_Rx)[3] &
                   varianza_Rx<=quantile(Matriz_medias_ExpYCeres$varianza_Rx)[4]))
write.table(grupo2, "grupo2.txt", quote = F, sep = "\t", row.names = F, col.names = T)

#GRUPO 3:
grupo3<-data.frame(Matriz_medias_ExpYCeres%>%filter(media_Rx>quantile(Matriz_medias_ExpYCeres$media_Rx)[3]&media_Rx<quantile(Matriz_medias_ExpYCeres$media_Rx)[4] &
                                           varianza_Rx<=quantile(Matriz_medias_ExpYCeres$varianza_Rx)[4]))
write.table(grupo3, "grupo3.txt", quote = F, sep = "\t", row.names = F, col.names = T)

#GRUPO 4:
grupo4<-data.frame(Matriz_medias_ExpYCeres%>%filter(media_Rx>quantile(Matriz_medias_ExpYCeres$media_Rx)[4] &
                                           varianza_Rx<=quantile(Matriz_medias_ExpYCeres$varianza_Rx)[4]))
write.table(grupo4, "grupo4.txt", quote = F, sep = "\t", row.names = F, col.names = T)

#GRUPO 5:
grupo5 <-data.frame(Matriz_medias_ExpYCeres%>%filter(media_Rx<quantile(Matriz_medias_ExpYCeres$media_Rx)[2] & 
                                            varianza_Rx>quantile(Matriz_medias_ExpYCeres$varianza_Rx)[4]&varianza_Rx<upper_whisker))
write.table(grupo5, "grupo5.txt", quote = F, sep = "\t", row.names = F, col.names = T)

#GRUPO 6:
grupo6<-data.frame(Matriz_medias_ExpYCeres%>%filter(media_Rx>quantile(Matriz_medias_ExpYCeres$media_Rx)[2]&media_Rx<quantile(Matriz_medias_ExpYCeres$media_Rx)[3] &
                                           varianza_Rx>quantile(Matriz_medias_ExpYCeres$varianza_Rx)[4]&varianza_Rx<upper_whisker))
write.table(grupo6, "grupo6.txt", quote = F, sep = "\t", row.names = F, col.names = T)

#GRUPO 7:
grupo7<-data.frame(Matriz_medias_ExpYCeres%>%filter(media_Rx>quantile(Matriz_medias_ExpYCeres$media_Rx)[3]&media_Rx<quantile(Matriz_medias_ExpYCeres$media_Rx)[4] &
                                           varianza_Rx>quantile(Matriz_medias_ExpYCeres$varianza_Rx)[4]&varianza_Rx<upper_whisker))
write.table(grupo7, "grupo7.txt", quote = F, sep = "\t", row.names = F, col.names = T)

#GRUPO 8:
grupo8<-data.frame(Matriz_medias_ExpYCeres%>%filter(media_Rx>quantile(Matriz_medias_ExpYCeres$media_Rx)[4] &
                                           varianza_Rx>quantile(Matriz_medias_ExpYCeres$varianza_Rx)[4]&varianza_Rx<upper_whisker))
write.table(grupo8, "grupo8.txt", quote = F, sep = "\t", row.names = F, col.names = T)

#GRUPO 9:
grupo9 <-data.frame(Matriz_medias_ExpYCeres%>%filter(media_Rx<quantile(Matriz_medias_ExpYCeres$media_Rx)[2] & 
                                            varianza_Rx>upper_whisker))
write.table(grupo9, "grupo9.txt", quote = F, sep = "\t", row.names = F, col.names = T)

#GRUPO 10:
grupo10<-data.frame(Matriz_medias_ExpYCeres%>%filter(media_Rx>quantile(Matriz_medias_ExpYCeres$media_Rx)[2]&media_Rx<quantile(Matriz_medias_ExpYCeres$media_Rx)[3] &
                                           varianza_Rx>upper_whisker))
write.table(grupo10, "grupo10.txt", quote = F, sep = "\t", row.names = F, col.names = T)

#GRUPO 11:
grupo11<-data.frame(Matriz_medias_ExpYCeres%>%filter(media_Rx>quantile(Matriz_medias_ExpYCeres$media_Rx)[3]&media_Rx<quantile(Matriz_medias_ExpYCeres$media_Rx)[4] &
                                           varianza_Rx>upper_whisker))
write.table(grupo11, "grupo11.txt", quote = F, sep = "\t", row.names = F, col.names = T)

#GRUPO 12:
grupo12<-data.frame(Matriz_medias_ExpYCeres%>%filter(media_Rx>quantile(Matriz_medias_ExpYCeres$media_Rx)[4] &
                                           varianza_Rx>upper_whisker))
write.table(grupo12, "grupo12.txt", quote = F, sep = "\t", row.names = F, col.names = T)

#save.image("./20190805_GPR.RData")
#~~~~~~~~~~~~~~~~~~~~~~~~~~ Representación de los grupos ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GRUPO 1:
png("EjemGrupo1.png", width = 13, height = 13, units = "in", res = 110)
ggplot(GP_crudo_t, aes(x=GP_crudo_t$ALADGLNexR)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#0d1c32") + 
  geom_density(fill = "#d05046", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(GP_crudo_t$ALADGLNexR)), color="black", linetype="dashed", size=1)+
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores Expresión", y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS EXPRESION TRAS GPR EN LA REACCION ALADGLNexR")
dev.off()


png("EjemGrupo1Ceres.png", width = 13, height = 13, units = "in", res = 110)
ggplot(data=GP_crudo_ceres_t, aes(x=GP_crudo_ceres_t$ALADGLNexR)) +
  geom_density(color="darkblue", fill="lightblue", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=1) +
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS ESENCIABILIDAD TRAS GPR EN LA REACCION ALADGLNexR")
dev.off()

png("EjemGrupo1_1.png", width = 13, height = 13, units = "in", res = 110)
ggplot(GP_crudo_t, aes(x=GP_crudo_t$CHOLATEt2)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#0d1c32") + 
  geom_density(fill = "#d05046", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(GP_crudo_t$CHOLATEt2)), color="black", linetype="dashed", size=1)+
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores Expresión", y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS EXPRESION TRAS GPR EN LA REACCION CHOLATEt2")
dev.off()


png("EjemGrupo1_1Ceres.png", width = 13, height = 13, units = "in", res = 110)
ggplot(data=GP_crudo_ceres_t, aes(x=GP_crudo_ceres_t$CHOLATEt2)) +
  geom_density(color="darkblue", fill="lightblue", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=1) +
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS ESENCIABILIDAD TRAS GPR EN LA REACCION CHOLATEt2")
dev.off()

# GRUPO 2:
png("EjemGrupo2.png", width = 13, height = 13, units = "in", res = 110)
ggplot(GP_crudo_t, aes(x=GP_crudo_t$r0470)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#0d1c32") + 
  geom_density(fill = "#d05046", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(GP_crudo_t$r0470)), color="black", linetype="dashed", size=1)+
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores Expresión", y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS EXPRESION TRAS GPR EN LA REACCION r0470")
dev.off()


png("EjemGrupo2Ceres.png", width = 13, height = 13, units = "in", res = 110)
ggplot(data=GP_crudo_ceres_t, aes(x=GP_crudo_ceres_t$r0470)) +
  geom_density(color="darkblue", fill="lightblue", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=1) +
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS ESENCIABILIDAD TRAS GPR EN LA REACCION r0470")
dev.off()

png("EjemGrupo2_1.png", width = 13, height = 13, units = "in", res = 110)
ggplot(GP_crudo_t, aes(x=GP_crudo_t$GLYATm)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#0d1c32") + 
  geom_density(fill = "#d05046", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(GP_crudo_t$GLYATm)), color="black", linetype="dashed", size=1)+
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores Expresión", y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS EXPRESION TRAS GPR EN LA REACCION GLYATm")
dev.off()


png("EjemGrupo2_1Ceres.png", width = 13, height = 13, units = "in", res = 110)
ggplot(data=GP_crudo_ceres_t, aes(x=GP_crudo_ceres_t$GLYATm)) +
  geom_density(color="darkblue", fill="lightblue", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=1) +
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS ESENCIABILIDAD TRAS GPR EN LA REACCION GLYATm")
dev.off()

# GRUPO 3:
png("EjemGrupo3.png", width = 13, height = 13, units = "in", res = 110)
ggplot(GP_crudo_t, aes(x=GP_crudo_t$DMATT)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#0d1c32") + 
  geom_density(fill = "#d05046", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(GP_crudo_t$DMATT)), color="black", linetype="dashed", size=1)+
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores Expresión", y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS EXPRESION TRAS GPR EN LA REACCION DMATT")
dev.off()


png("EjemGrupo3Ceres.png", width = 13, height = 13, units = "in", res = 110)
ggplot(data=GP_crudo_ceres_t, aes(x=GP_crudo_ceres_t$DMATT)) +
  geom_density(color="darkblue", fill="lightblue", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=1) +
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS ESENCIABILIDAD TRAS GPR EN LA REACCION DMATT")
dev.off()

png("EjemGrupo3_1.png", width = 13, height = 13, units = "in", res = 110)
ggplot(GP_crudo_t, aes(x=GP_crudo_t$GLYOX)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#0d1c32") + 
  geom_density(fill = "#d05046", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(GP_crudo_t$GLYOX)), color="black", linetype="dashed", size=1)+
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores Expresión", y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS EXPRESION TRAS GPR EN LA REACCION GLYOX")
dev.off()


png("EjemGrupo3_1Ceres.png", width = 13, height = 13, units = "in", res = 110)
ggplot(data=GP_crudo_ceres_t, aes(x=GP_crudo_ceres_t$GLYOX)) +
  geom_density(color="darkblue", fill="lightblue", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=1) +
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS ESENCIABILIDAD TRAS GPR EN LA REACCION GLYOX")
dev.off()

# GRUPO 4:
png("EjemGrupo4.png", width = 13, height = 13, units = "in", res = 110)
ggplot(GP_crudo_t, aes(x=GP_crudo_t$DTMPK)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#0d1c32") + 
  geom_density(fill = "#d05046", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(GP_crudo_t$DTMPK)), color="black", linetype="dashed", size=1)+
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores Expresión", y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS EXPRESION TRAS GPR EN LA REACCION DTMPK")
dev.off()


png("EjemGrupo4Ceres.png", width = 13, height = 13, units = "in", res = 110)
ggplot(data=GP_crudo_ceres_t, aes(x=GP_crudo_ceres_t$DTMPK)) +
  geom_density(color="darkblue", fill="lightblue", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=1) +
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS ESENCIABILIDAD TRAS GPR EN LA REACCION DTMPK")
dev.off()

png("EjemGrupo4_1.png", width = 13, height = 13, units = "in", res = 110)
ggplot(GP_crudo_t, aes(x=GP_crudo_t$GLXO1)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#0d1c32") + 
  geom_density(fill = "#d05046", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(GP_crudo_t$GLXO1)), color="black", linetype="dashed", size=1)+
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores Expresión", y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS EXPRESION TRAS GPR EN LA REACCION GLXO1")
dev.off()


png("EjemGrupo4_1Ceres.png", width = 13, height = 13, units = "in", res = 110)
ggplot(data=GP_crudo_ceres_t, aes(x=GP_crudo_ceres_t$GLXO1)) +
  geom_density(color="darkblue", fill="lightblue", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=1) +
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS ESENCIABILIDAD TRAS GPR EN LA REACCION GLXO1")
dev.off()

# GRUPO 5:
png("EjemGrupo5.png", width = 13, height = 13, units = "in", res = 110)
ggplot(GP_crudo_t, aes(x=GP_crudo_t$SERPT)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#0d1c32") + 
  geom_density(fill = "#d05046", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(GP_crudo_t$SERPT)), color="black", linetype="dashed", size=1)+
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores Expresión", y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS EXPRESION TRAS GPR EN LA REACCION SERPT")
dev.off()


png("EjemGrupo5Ceres.png", width = 13, height = 13, units = "in", res = 110)
ggplot(data=GP_crudo_ceres_t, aes(x=GP_crudo_ceres_t$SERPT)) +
  geom_density(color="darkblue", fill="lightblue", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=1) +
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS ESENCIABILIDAD TRAS GPR EN LA REACCION SERPT")
dev.off()

png("EjemGrupo5_1.png", width = 13, height = 13, units = "in", res = 110)
ggplot(GP_crudo_t, aes(x=GP_crudo_t$FUT12g)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#0d1c32") + 
  geom_density(fill = "#d05046", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(GP_crudo_t$FUT12g)), color="black", linetype="dashed", size=1)+
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores Expresión", y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS EXPRESION TRAS GPR EN LA REACCION FUT12g")
dev.off()


png("EjemGrupo5_1Ceres.png", width = 13, height = 13, units = "in", res = 110)
ggplot(data=GP_crudo_ceres_t, aes(x=GP_crudo_ceres_t$FUT12g)) +
  geom_density(color="darkblue", fill="lightblue", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=1) +
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS ESENCIABILIDAD TRAS GPR EN LA REACCION FUT12g")
dev.off()

# GRUPO 6:
png("EjemGrupo6.png", width = 13, height = 13, units = "in", res = 110)
ggplot(GP_crudo_t, aes(x=GP_crudo_t$C14STRr)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#0d1c32") + 
  geom_density(fill = "#d05046", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(GP_crudo_t$C14STRr)), color="black", linetype="dashed", size=1)+
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores Expresión", y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS EXPRESION TRAS GPR EN LA REACCION C14STRr")
dev.off()


png("EjemGrupo6Ceres.png", width = 13, height = 13, units = "in", res = 110)
ggplot(data=GP_crudo_ceres_t, aes(x=GP_crudo_ceres_t$C14STRr)) +
  geom_density(color="darkblue", fill="lightblue", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=1) +
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS ESENCIABILIDAD TRAS GPR EN LA REACCION C14STRr")
dev.off()

png("EjemGrupo6_1.png", width = 13, height = 13, units = "in", res = 110)
ggplot(GP_crudo_t, aes(x=GP_crudo_t$r0120)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#0d1c32") + 
  geom_density(fill = "#d05046", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(GP_crudo_t$r0120)), color="black", linetype="dashed", size=1)+
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores Expresión", y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS EXPRESION TRAS GPR EN LA REACCION r0120")
dev.off()


png("EjemGrupo6_1Ceres.png", width = 13, height = 13, units = "in", res = 110)
ggplot(data=GP_crudo_ceres_t, aes(x=GP_crudo_ceres_t$r0120)) +
  geom_density(color="darkblue", fill="lightblue", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=1) +
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS ESENCIABILIDAD TRAS GPR EN LA REACCION r0120")
dev.off()

# GRUPO 7:
png("EjemGrupo7.png", width = 13, height = 13, units = "in", res = 110)
ggplot(GP_crudo_t, aes(x=GP_crudo_t$UGCG)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#0d1c32") + 
  geom_density(fill = "#d05046", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(GP_crudo_t$UGCG)), color="black", linetype="dashed", size=1)+
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores Expresión", y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS EXPRESION TRAS GPR EN LA REACCION UGCG")
dev.off()


png("EjemGrupo7Ceres.png", width = 13, height = 13, units = "in", res = 110)
ggplot(data=GP_crudo_ceres_t, aes(x=GP_crudo_ceres_t$UGCG)) +
  geom_density(color="darkblue", fill="lightblue", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=1) +
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS ESENCIABILIDAD TRAS GPR EN LA REACCION UGCG")
dev.off()

png("EjemGrupo7_1.png", width = 13, height = 13, units = "in", res = 110)
ggplot(GP_crudo_t, aes(x=GP_crudo_t$UGCG)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#0d1c32") + 
  geom_density(fill = "#d05046", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(GP_crudo_t$UGCG)), color="black", linetype="dashed", size=1)+
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores Expresión", y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS EXPRESION TRAS GPR EN LA REACCION UGCG")
dev.off()


png("EjemGrupo7_1Ceres.png", width = 13, height = 13, units = "in", res = 110)
ggplot(data=GP_crudo_ceres_t, aes(x=GP_crudo_ceres_t$FAOXC101x)) +
  geom_density(color="darkblue", fill="lightblue", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=1) +
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS ESENCIABILIDAD TRAS GPR EN LA REACCION FAOXC101x")
dev.off()

# GRUPO 8:
png("EjemGrupo8.png", width = 13, height = 13, units = "in", res = 110)
ggplot(GP_crudo_t, aes(x=GP_crudo_t$HMGCOASi)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#0d1c32") + 
  geom_density(fill = "#d05046", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(GP_crudo_t$HMGCOASi)), color="black", linetype="dashed", size=1)+
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores Expresión", y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS EXPRESION TRAS GPR EN LA REACCION HMGCOASi")
dev.off()

png("EjemGrupo8Ceres.png", width = 13, height = 13, units = "in", res = 110)
ggplot(data=GP_crudo_ceres_t, aes(x=GP_crudo_ceres_t$HMGCOASi)) +
  geom_density(color="darkblue", fill="lightblue", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=1) +
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS ESENCIABILIDAD TRAS GPR EN LA REACCION HMGCOASi")
dev.off()

png("EjemGrupo8_1.png", width = 13, height = 13, units = "in", res = 110)
ggplot(GP_crudo_t, aes(x=GP_crudo_t$GTHDH)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#0d1c32") + 
  geom_density(fill = "#d05046", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(GP_crudo_t$GTHDH)), color="black", linetype="dashed", size=1)+
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores Expresión", y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS EXPRESION TRAS GPR EN LA REACCION GTHDH")
dev.off()

png("EjemGrupo8_1Ceres.png", width = 13, height = 13, units = "in", res = 110)
ggplot(data=GP_crudo_ceres_t, aes(x=GP_crudo_ceres_t$GTHDH)) +
  geom_density(color="darkblue", fill="lightblue", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=1) +
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS ESENCIABILIDAD TRAS GPR EN LA REACCION GTHDH")
dev.off()

# GRUPO 9:
png("EjemGrupo9.png", width = 13, height = 13, units = "in", res = 110)
ggplot(GP_crudo_t, aes(x=GP_crudo_t$PROD2m)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#0d1c32") + 
  geom_density(fill = "#d05046", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(GP_crudo_t$PROD2m)), color="black", linetype="dashed", size=1)+
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores Expresión", y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS EXPRESION TRAS GPR EN LA REACCION PROD2m")
dev.off()

png("EjemGrupo9Ceres.png", width = 13, height = 13, units = "in", res = 110)
ggplot(data=GP_crudo_ceres_t, aes(x=GP_crudo_ceres_t$PROD2m)) +
  geom_density(color="darkblue", fill="lightblue", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=1) +
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS ESENCIABILIDAD TRAS GPR EN LA REACCION PROD2m")
dev.off()

png("EjemGrupo9_1.png", width = 13, height = 13, units = "in", res = 110)
ggplot(GP_crudo_t, aes(x=GP_crudo_t$HISDC)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#0d1c32") + 
  geom_density(fill = "#d05046", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(GP_crudo_t$HISDC)), color="black", linetype="dashed", size=1)+
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores Expresión", y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS EXPRESION TRAS GPR EN LA REACCION HISDC")
dev.off()

png("EjemGrupo9_1Ceres.png", width = 13, height = 13, units = "in", res = 110)
ggplot(data=GP_crudo_ceres_t, aes(x=GP_crudo_ceres_t$HISDC)) +
  geom_density(color="darkblue", fill="lightblue", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=1) +
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS ESENCIABILIDAD TRAS GPR EN LA REACCION HISDC")
dev.off()

# GRUPO 10:
png("EjemGrupo10.png", width = 13, height = 13, units = "in", res = 110)
ggplot(GP_crudo_t, aes(x=GP_crudo_t$CYSGLYexR)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#0d1c32") + 
  geom_density(fill = "#d05046", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(GP_crudo_t$CYSGLYexR)), color="black", linetype="dashed", size=1)+
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores Expresión", y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS EXPRESION TRAS GPR EN LA REACCION CYSGLYexR")
dev.off()

png("EjemGrupo10Ceres.png", width = 13, height = 13, units = "in", res = 110)
ggplot(data=GP_crudo_ceres_t, aes(x=GP_crudo_ceres_t$CYSGLYexR)) +
  geom_density(color="darkblue", fill="lightblue", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=1) +
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS ESENCIABILIDAD TRAS GPR EN LA REACCION CYSGLYexR")
dev.off()

png("EjemGrupo10_1.png", width = 13, height = 13, units = "in", res = 110)
ggplot(GP_crudo_t, aes(x=GP_crudo_t$ESTRONEGLCt)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#0d1c32") + 
  geom_density(fill = "#d05046", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(GP_crudo_t$ESTRONEGLCt)), color="black", linetype="dashed", size=1)+
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores Expresión", y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS EXPRESION TRAS GPR EN LA REACCION ESTRONEGLCt")
dev.off()

png("EjemGrupo10_1Ceres.png", width = 13, height = 13, units = "in", res = 110)
ggplot(data=GP_crudo_ceres_t, aes(x=GP_crudo_ceres_t$ESTRONEGLCt)) +
  geom_density(color="darkblue", fill="lightblue", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=1) +
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS ESENCIABILIDAD TRAS GPR EN LA REACCION ESTRONEGLCt")
dev.off()

# GRUPO 11:
png("EjemGrupo11.png", width = 13, height = 13, units = "in", res = 110)
ggplot(GP_crudo_t, aes(x=GP_crudo_t$ARGSS)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#0d1c32") + 
  geom_density(fill = "#d05046", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(GP_crudo_t$ARGSS)), color="black", linetype="dashed", size=1)+
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores Expresión", y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS EXPRESION TRAS GPR EN LA REACCION ARGSS")
dev.off()

png("EjemGrupo11Ceres.png", width = 13, height = 13, units = "in", res = 110)
ggplot(data=GP_crudo_ceres_t, aes(x=GP_crudo_ceres_t$ARGSS)) +
  geom_density(color="darkblue", fill="lightblue", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=1) +
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS ESENCIABILIDAD TRAS GPR EN LA REACCION ARGSS")
dev.off()

png("EjemGrupo11_1.png", width = 13, height = 13, units = "in", res = 110)
ggplot(GP_crudo_t, aes(x=GP_crudo_t$MACOXO)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#0d1c32") + 
  geom_density(fill = "#d05046", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(GP_crudo_t$MACOXO)), color="black", linetype="dashed", size=1)+
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores Expresión", y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS EXPRESION TRAS GPR EN LA REACCION MACOXO")
dev.off()

png("EjemGrupo11_1Ceres.png", width = 13, height = 13, units = "in", res = 110)
ggplot(data=GP_crudo_ceres_t, aes(x=GP_crudo_ceres_t$MACOXO)) +
  geom_density(color="darkblue", fill="lightblue", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=1) +
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS ESENCIABILIDAD TRAS GPR EN LA REACCION MACOXO")
dev.off()

# GRUPO 12:
png("EjemGrupo12.png", width = 13, height = 13, units = "in", res = 110)
ggplot(GP_crudo_t, aes(x=GP_crudo_t$PSERT)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#0d1c32") + 
  geom_density(fill = "#d05046", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(GP_crudo_t$PSERT)), color="black", linetype="dashed", size=1)+
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores Expresión", y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS EXPRESION TRAS GPR EN LA REACCION PSERT")
dev.off()

png("EjemGrupo12Ceres.png", width = 13, height = 13, units = "in", res = 110)
ggplot(data=GP_crudo_ceres_t, aes(x=GP_crudo_ceres_t$PSERT)) +
  geom_density(color="darkblue", fill="lightblue", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=1) +
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS ESENCIABILIDAD TRAS GPR EN LA REACCION PSERT")
dev.off()

png("EjemGrupo12_1.png", width = 13, height = 13, units = "in", res = 110)
ggplot(GP_crudo_t, aes(x=GP_crudo_t$LALDO2x)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#0d1c32") + 
  geom_density(fill = "#d05046", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(GP_crudo_t$LALDO2x)), color="black", linetype="dashed", size=1)+
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores Expresión", y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS EXPRESION TRAS GPR EN LA REACCION LALDO2x")
dev.off()

png("EjemGrupo12_1Ceres.png", width = 13, height = 13, units = "in", res = 110)
ggplot(data=GP_crudo_ceres_t, aes(x=GP_crudo_ceres_t$LALDO2x)) +
  geom_density(color="darkblue", fill="lightblue", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=1) +
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=14), plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN DATOS ESENCIABILIDAD TRAS GPR EN LA REACCION LALDO2x")
dev.off()
save.image("./20190805_GPR.RData")

#*******************************************************************************************************************************************
#GP con expresion muy baja y alta esenciabilidad.

library(dplyr)
library(reshape2)
cuartiles<-quantile(media_var_sd$media_Rx)
GP_mediaBaja_CeresAlto <-data.frame(Expr_ceres%>%filter(ceres<=(-0.9) & expresion<=cuartiles[2]))
write.table(GP_mediaBaja_CeresAlto, "GP_mediaBaja_CeresAlto.txt", quote = F, sep = "\t", row.names = F, col.names = T)
GP_mediaBaja_CeresAlto2 <-data.frame(Expr_ceres%>%filter(ceres<=(-0.9) & expresion>=1 & expresion<=cuartiles[2]))
div<-colsplit(GP_mediaBaja_CeresAlto2$linea, "_", c("ID","TISSUE"))
GP_mediaBaja_CeresAlto2<-cbind(div, GP_mediaBaja_CeresAlto2)
write.table(GP_mediaBaja_CeresAlto2, "GP_mediaBaja_CeresAlto2.txt", quote = F, sep = "\t", row.names = F, col.names = T)
setwd("/home/estrella/TFM/Datos/depmap/GPR/")

#Core de Rx que son esenciales en todas las cell_line (no me sale seleccionar columnas de GP_crudo_ceres_t), asi que voy a seleccionar de la matriz de las medias.

Core_ceres<-as.data.frame(Matriz_medias_ExpYCeres%>%filter(media_Rx_CERES<=(-0.5)))
write.table(Core_ceres, "Core_ceres.txt", quote = F, sep = "\t", row.names = F, col.names = T)

IQR <- IQR(Matriz_medias_ExpYCeres$varianza_Rx)
cuartiles<-quantile(Matriz_medias_ExpYCeres$media_Rx)


##########################################################################################################################
##########################################################################################################################
## ............................................. IDEAS THRESHOLDING ......................................................
##########################################################################################################################
#########################################################################################################################

#Dos thresholds globales (cuartiles de la distribucion de las medias): Uno por debajo de Q2 y otro por encima de Q3+(Q2-Q1)
cuantiles <-quantile(Matriz_medias_ExpYCeres$media_Rx, c(.05, .10, .90, .80, .70))
P05 <-cuantiles[1]
P10 <-cuantiles[2]
P90 <-cuantiles[3]
P80 <-cuantiles[4]
P70 <-cuantiles[5]
#TH_l <-cuartiles[2]
#TH_u <-(cuartiles[4]+(cuartiles[3]-cuartiles[2]))
# cuantiles_RXs<-as.data.frame(apply(GP_crudo_t,2,quantile)) 

cuartiles<-quantile(media_var_sd$media_Rx)
png("Histograma_con_TH.png", width = 20, height = 13, units = "in", res = 110)
ggplot(Matriz_medias_ExpYCeres,aes(x=media_Rx))+
  geom_histogram(binwidth=0.1, fill="blue",color="black")+
  geom_density(aes(y=0.47* ..count..), fill="darksalmon", colour="darksalmon", alpha=0.4)+
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.text=element_text(size=16), axis.title=element_text(size=28,face="bold"), plot.title=element_text(size=32, face="bold"))+
  geom_vline (aes(xintercept=P05), color="blue", linetype="dashed", size=0.8)+
  geom_text(aes(x=(P05+0.2), label="P05", y=305), colour="blue", size=7)+
  geom_vline (aes(xintercept=P10), color="red", linetype="dashed", size=0.8)+
  geom_text(aes(x=(P10+0.2), label="P10", y=290), colour="red", size=7)+
  geom_vline (aes(xintercept=P90), color="yellow", linetype="dashed", size=0.8)+
  geom_text(aes(x=(P90+0.2), label="P90", y=305), colour="yellow", size=7)+
  geom_vline (aes(xintercept=P80), color="green", linetype="dashed", size=0.8)+
  geom_text(aes(x=(P80+0.2), label="P80", y=290), colour="green", size=7)+
  geom_vline (aes(xintercept=P70), color="orange", linetype="dashed", size=0.8)+
  geom_text(aes(x=(P70+0.2), label="P70", y=305), colour="orange", size=7)+
  #geom_vline (aes(xintercept=cuartiles[3]), color="blue", size=0.8)
  labs(x='Media de expresión', y='Densidad')
dev.off()


library(dplyr)
RX_valen<-c(rownames(na.omit(GP_crudo_ceres)[c(1:10),c(1:10)]))
mprueba <- GP_crudo[rownames(GP_crudo)%in% as.character(RX_valen), c(1:10)]
Rx="13DAMPPOX"
CL="X143B_BONE" 
# x=mprueba
# lista_CL<-colnames(mprueba)
# lista_Rx<-rownames(mprueba)

  # for (CL in lista_CL){
  #   print(CL)
  #  for (Rx in lista_Rx){
  #    print(Rx)
  #     if (GP_crudo_esti[Rx,CL] <(as.numeric(TH_l))) {
  #       GP_crudo_esti[Rx, CL] <- 0
  #     }else{
  #         if (GP_crudo_esti[Rx,CL]>(as.numeric(TH_u))) {
  #           GP_crudo_esti[Rx, CL] <- 0
  #         }else{
  #             # if (GP_crudo_esti[Rx,CL]<(as.numeric(cuantiles_RXs[[Rx]][1])) 
  #             #     & GP_crudo_esti[Rx, CL]>(as.numeric(cuantiles_RXs[[Rx]][4]))) {
  #             if (GP_crudo_esti[Rx,CL]>media_var_sd[Rx,2]) {
  #               GP_crudo_esti[Rx,CL] <- 0
  #             }else{
  #               GP_crudo_esti[Rx,CL] <- 1
  #          }
  #         }
  #       }
  #     }
  #   }

# GP_crudo_esti_2<- GP_crudo_esti
# GP_crudo_esti_3<- GP_crudo_esti #Es la misma que arriba pero con 1, en lugar de menos 1.
# GP_crudo_esti_4<- GP_crudo_esti #Poniendo el Thlocal como lo que está por encima de la media de cada Rx.

#De momento seran TRES ZONAS: -1 (nada confidente, debajo del P_20), 0 (no sé confidencia), 1 (muy confidente, arriba P_90).
Percentiles_20<-as.numeric(quantile(media_var_sd$media_Rx, c(0.20)))
Percentiles_90<-as.numeric(quantile(media_var_sd$media_Rx, c(0.90)))

GP_crudo_esti <- na.omit(GP_crudo)
lista_CL<-colnames(GP_crudo_esti)
lista_Rx<-rownames(GP_crudo_esti)

for (CL in lista_CL){
  print(CL)
 for (Rx in lista_Rx){
   print(Rx)
    if (GP_crudo_esti[Rx,CL] <Percentiles_20) {
      GP_crudo_esti[Rx, CL] <- -1
    }else{
        if (GP_crudo_esti[Rx,CL]>Percentiles_90) {
          GP_crudo_esti[Rx, CL] <- 1
        }else{
              GP_crudo_esti[Rx,CL] <- 0
         }
        }
      }
    }


Estimacion_1<-GP_crudo_esti
write.table (Estimacion_1, file = "Estimacion_1.csv", col.names = TRUE, sep="\t")
setwd("/home/estrella/TFM/Datos/depmap/GPR/")
save.image("./20190809_GPR.RData")

library(reshape2)
Estimacion_1_vector<-melt(as.matrix(Estimacion_1))
colnames(Estimacion_1_vector)<-c("Rx", "CL", "valorTH")

GP_crudo_ceres<-na.omit(GP_crudo_ceres)
library(ggplot2)
#Distribución de los datos de esenciabilidad de todas las Rx para una cell line específica (X143B_BONE).
#Entiendo que sobre esta representación, si mi nuemro de -1 es similar a lo que es esencial entonces estoy estimando bien.


png("Distribucion_esenciabilidad_X42MGBA_CENTRAL_NERVOUS_SYSTEM.png", width = 13, height = 13, units = "in", res = 110)
ggplot(GP_crudo_ceres, aes(GP_crudo_ceres$X42MGBA_CENTRAL_NERVOUS_SYSTEM))+
  # geom_histogram(binwidth=0.2, fill="lightblue",color="black")+
  geom_density(aes(y=0.2* ..count..), fill="blue", colour="blue", alpha=0.4)+
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.text=element_text(size=16), axis.title=element_text(size=28,face="bold"), plot.title=element_text(size=32, face="bold"))+
  labs(x='Esenciabilidad', y='Densidad')+
  ggtitle("Distribución de la esenciabilidad de la todas las Rx para la linea X42MGBA_CENTRAL_NERVOUS_SYSTEM")
dev.off()


png("Estimacion_4_X42MGBA_CENTRAL_NERVOUS_SYSTEM.png", width = 13, height = 13, units = "in", res = 110)
ggplot(Estimacion_1, aes(Estimacion_1$X253J_URINARY_TRACT))+
  # geom_histogram(binwidth=0.2, fill="lightblue",color="black")+
  geom_density(aes(y=0.2* ..count..), fill="blue", colour="blue", alpha=0.4)+
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.text=element_text(size=16), axis.title=element_text(size=28,face="bold"), plot.title=element_text(size=32, face="bold"))+
  labs(x='Esenciabilidad', y='Densidad')+
ggtitle("Estimación para todas las Rx para la linea X42MGBA_CENTRAL_NERVOUS_SYSTEM")
dev.off()


# wilcox.test ( GP_crudo_ceres$X42MGBA_CENTRAL_NERVOUS_SYSTEM ~ #Con que la compara, lo que me sale de 0 y -1 no es continuo.)

# #~~~~~~~~~~~~~~~~~ Plot zonas confidencias ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ggplot(GP_crudo, (aes(x=GP_crudo$X42MGBA_CENTRAL_NERVOUS_SYSTEM)))+
#   geom_histogram(binwidth=0.9, fill="lightblue",color="black")+
#   geom_line(aes(y = ..density.., colour = 'Blue'), stat = 'density')
# 
# GP_crudo<-na.omit(GP_crudo)
# ggplot(GP_crudo) +
#   geom_histogram(aes(x=X42MGBA_CENTRAL_NERVOUS_SYSTEM, y=..density..), bins=30)+
#   geom_density(color="blue")
# dev.off()


ggplot()+
  geom_histogram(GP_crudo, aes(x=colnames(GP_crudo)), fill="lightblue",color="black")
