getwd()
setwd("/home/estrella/TFM/Datos/depmap/AnalisisGenes/")
matrizDatosCrudos<- read.table (file='/home/estrella/TFM/Datos/depmap/M19Q2_etiquetada.csv', header = TRUE, sep='\t')

#####################################################################################################################################################
######################  ------------------- CLASIFICACION DE GENES SEGUN MEDIA Y VARIANZA-----------------------#####################################
#####################################################################################################################################################

#OBJETIVO: Olvidarnos de la clasificación por tejidos. Mirar las cell line en conjunto y calcular la media y varianza, para cada gen. Para despúes 
#          standarizar las varianzas y poder comparar los genes entre si. Así podremos clasificarlos en: High, moderate y low confidence, para meterlos 
#          en el modelo metabolico.

library(dplyr)
#Calculo media de cada gen.
Datos2 <- as.data.frame(lapply(matrizDatosCrudos[,-c(1)], as.numeric)) #Quito primera con el nombre de las cell line.
media_genes<-as.data.frame(apply(Datos2,2,mean))
colnames(media_genes)<- 'media_gen'

#Calculo de la varianza de cada gen.
varianza_genes <-as.data.frame(apply(Datos2,2,var))
sd_genes <-as.data.frame(apply(Datos2,2,sd))
colnames(varianza_genes)<-'varianza_gen'
colnames(sd_genes)<-'sd_gen'
#Uno media y varianzas.
media_var <- cbind(varianza_genes, media_genes,sd_genes)
media_var<-cbind(rownames(media_var), media_var)
colnames(media_var[1]) <-'genes'

summary (media_var)
which.max(media_var$varianza_gen)
media_var[which.max(media_var$varianza_gen),1] #Gen con la varianza mas alta
media_var[which.max(media_var$media_gen),1] #Gen con la media mas alta

#Plot distribucion de genes según su media
library(ggplot2)
GENES <- media_var[,1]
h<-hist(media_var$media_gen,breaks = 500, plot = FALSE) #De esta forma agrupa los genes de 500 en 500 por media y mira la frecuencia de cada grupo.
plot(h, xaxt = "n", xlab = "Genes", ylab = "Frecuencia",
     main = "", col = "lightblue")
save.image("./20190710_genes.RData")

#Plot media vs varianza
Media<- media_var$media_gen
Varianza<- media_var$varianza_gen
SD<- media_var$sd_gen
png("MediavsVarianza.png", width = 13, height = 13, units = "in", res = 110)
ggplot(media_var, aes(x=Media,y=Varianza)) +
  geom_point(colour = "darkblue")+
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.text=element_text(size=16), axis.title=element_text(size=28,face="bold"), plot.title=element_text(size=32, face="bold"))+
  labs(x='Media', y='Varianza')+
  ggtitle("Representación media vs mediana")+
  ylim(c(0,18))+
  xlim(c(0,18))
dev.off()

library(ggplot2)
library(ggExtra)

#Esto es para integrar al diagrama de disersión los boxplot para hacer grupos.
scatter<-ggplot(media_var, aes(x=Media,y=Varianza)) +
  geom_point(colour = "darkblue")+
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.text=element_text(size=16), axis.title=element_text(size=28,face="bold"), plot.title=element_text(size=32, face="bold"))+
  labs(x='Media', y='Varianza')+
  ylim(c(0,18))+
  xlim(c(0,18))

ggMarginal(scatter, type="boxplot", size=10, colour = 'black', fill = 'darkgrey', xparams = list(size=0.5), yparams = list(size=0.5))
   
#Plot media vs SD
ggplot(media_var, aes(x=Media,y=SD)) +
  geom_point(colour = "darkblue")+
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))+
  labs(x='Media', y='SD')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################################### Superposicion de las dos graficas MediaVsVarianza     ##################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MatrizMetabolismo2<- read.table('/home/estrella/TFM/Datos/depmap/M19Q2_metabol.csv', header = TRUE, sep= '\t')

media_genes<-as.data.frame(apply(MatrizMetabolismo2[-c(1)],2,mean))
colnames(media_genes)<- 'media_gen'

#Calculo de la varianza de cada gen. METABOLICOS
varianza_genes <-as.data.frame(apply(MatrizMetabolismo2[-c(1)],2,var))
sd_genes <-as.data.frame(apply(MatrizMetabolismo2[-c(1)],2,sd))
colnames(varianza_genes)<-'varianza_gen'
colnames(sd_genes)<-'sd_gen'
#Uno media y varianzas.
media_var_metabol <- cbind(varianza_genes, media_genes,sd_genes)
media_var_metabol<-cbind(rownames(media_var_metabol), media_var_metabol)
colnames(media_var[1]) <-'genes'
CV<- as.data.frame(media_var_metabol$sd_gen/media_var_metabol$media_gen)
colnames(CV)<-"cv"
media_var_metabol<-cbind(media_var_metabol, CV)
write.table(media_var_metabol,"/home/estrella/TFM/Datos/depmap/GPR/media_var_metabol.csv", quote = T, sep='\t')

################## Plot mediavsVarianza de los genes metabolicos ################################################
library(ggplot2)
Media2<- media_var_metabol$media_gen
Varianza2<- media_var_metabol$varianza_gen
SD2<- media_var_metabol$sd_gen
png("MediavsVarianza_metabol.png", width = 13, height = 13, units = "in", res = 110)
ggplot(media_var_metabol, aes(x=Media2,y=Varianza2)) +
  geom_point(colour = "darkred")+
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.text=element_text(size=16), axis.title=element_text(size=28,face="bold"), plot.title=element_text(size=32, face="bold"))+
  labs(x='Media', y='Varianza')+
  ggtitle("Representación media vs mediana")+
  ylim(c(0,18))+
  xlim(c(0,18))
dev.off()

################ Plot coeficiente de variación #####################################################################
gen<- media_var_metabol$`rownames(media_var_metabol)`
cv <-media_var_metabol$cv
png("Coeficiente de variacion.png", width = 13, height = 13, units = "in", res = 110)
ggplot(media_var_metabol, aes(x=cv)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#0d1c32") + 
  geom_density(fill = "#d05046", alpha = 0.3) + 
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=20), plot.title = element_text(size=22, face="bold"))+
  labs(x="Valores coeficiente de variación", y="Densidad")+
  ggtitle("Distribución del coeficiente de variación")
dev.off()

library(ggrepel)
png("MediavsCoeficienteVariacion_metabol.png", width = 13, height = 13, units = "in", res = 110)
ggplot(media_var_metabol, aes(x=Media2,y=cv)) +
  geom_point(colour = "darkred")+
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.text=element_text(size=16), axis.title=element_text(size=28,face="bold"), plot.title=element_text(size=32, face="bold"))+
  labs(x='Media', y='Coeficiente Variacion')+
  ggtitle("Representación media vs Coeficiente de variación")
dev.off()

muyvariantes<-media_var_metabol[media_var_metabol$cv>10,]
write.table(muyvariantes,"genes_Con_CV_muyAlto.txt")

################ Plot mediaVsVarianza METABOLICOS con boxplot  #####################################################################
IQR <- IQR(media_var_metabol$varianza_gen)
upper_whisker = (quantile(media_var_metabol$varianza_gen)[4]+ 1.5 * IQR)
png("MediavsVarianza_metabol_boxplot.png", width = 13, height = 13, units = "in", res = 110)
scatter2 <- ggplot(media_var_metabol, aes(x=Media2,y=Varianza2)) +
  geom_point(colour = "darkred")+
  geom_vline (aes(xintercept=quantile(media_var_metabol$media_gen)[2]), color="blue", linetype="dashed", size=1.5)+
  geom_vline (aes(xintercept=quantile(media_var_metabol$media_gen)[3]), color="blue", linetype="dashed", size=1.5)+
  geom_vline (aes(xintercept=quantile(media_var_metabol$media_gen)[4]), color="blue", linetype="dashed", size=1.5)+
  geom_hline (aes(yintercept=quantile(media_var_metabol$varianza_gen)[4]), color="blue", linetype="dashed", size=1.5)+
  geom_hline (aes(yintercept=upper_whisker), color="blue", linetype="dashed", size=1.5)+
  theme (panel.background = element_rect(fill='lightgrey', colour='darkgrey'),axis.text.x = element_text(size=20, face="bold"),
         axis.text.y = element_text(size=20, face="bold"), axis.title=element_text(size=28,face="bold"))+
  labs(x='Mean', y='Variance')
ggMarginal(scatter2, type="boxplot", size=10, colour = 'black', fill = 'darkgrey', xparams = list(size=0.5), yparams = list(size=0.5))
dev.off()

#######################################################################################################################################################
#######################################################################################################################################################
############################################    BARPLOT DE ESPRESION FRENTE ESENCIABILIDAD   ##########################################################
#######################################################################################################################################################
#######################################################################################################################################################

MMetabolismo<- MatrizMetabolismo2[,-c(1)]
MMetabolismo<-MMetabolismo[ , order(names(MMetabolismo))]
rownames(MMetabolismo)<-MatrizMetabolismo2$Cell_line
t<-data.frame(colnames(MMetabolismo)) # Genes ordenados que salen de la matriz expresion con genes metabolicos.
t2<-data.frame(colnames(Ceres_metabol)) #Genes ordenados que salen de la matriz de ceres con genes metabolicos. 
#Como la matriz de ceres tiene más genes metabólicos, quitamos aquellos que no están en los daros de expresión(MMetabolismo).
Ceres_metabol2<-data.frame(Ceres_metabol[,colnames(Ceres_metabol) %in% as.character(t$colnames.MMetabolismo.)])
t3<-data.frame(colnames(Ceres_metabol2))
MMetabolismo2<-data.frame(MMetabolismo[,colnames(MMetabolismo) %in% as.character(t3$colnames.Ceres_metabol2.)])#Cuando compruebas el numero de genes con t3, hay menos genes en ceres que la de expresion. Por eso le quitas los genes que no están 
MMetabolismo2<-MMetabolismo2[rownames(MMetabolismo2) %in% as.character(rownames(Ceres_metabol2)),]#El número de lineas en Ceres tb es menor, por eso te quedas solo con las que coinciden. 
#Salen 4 lineas menos en la de expresión por eso, lo corriges.
Ceres_metabol2<-Ceres_metabol2[rownames(Ceres_metabol2) %in% as.character(rownames(MMetabolismo2)),]

MMetabolismo3<-MMetabolismo2[order(rownames(MMetabolismo2)) ,]#Para ordenar de forma alfabetica las filas de las dos matrices.
Ceres_metabol3<-Ceres_metabol2[order(rownames(Ceres_metabol2)) ,]

library(reshape2)
Expresion<-melt(as.matrix(MMetabolismo3))
colnames(Expresion)<-c("linea", "gen", "expresion")
Esencial<-melt(as.matrix(Ceres_metabol3))
colnames(Esencial)<-c("linea", "gen", "ceres")
Expr_ceres<-data.frame()
Expr_ceres<-cbind(Expresion, Esencial[,3])
colnames(Expr_ceres) <-c("linea", "gen", "expresion", "ceres")


library(ggplot2)
library(ggrepel)
png("ExpresionVsCeres.png", width = 13, height = 13, units = "in", res = 110)
ggplot(Expr_ceres, aes(x=expresion , y=ceres)) +
  geom_point(colour = "darkblue")+
  ggtitle("Distribución expresion vs esenciabilidad")+
  theme(panel.background = element_rect(fill='lightgrey', colour='darkgrey'),
        axis.text=element_text(size=18), axis.title=element_text(size=30,face="bold"), plot.title=element_text(size=34, face="bold"))+ 
  xlab("Expresion")+
  ylab("Esenciabilidad")+
  geom_hline(yintercept = -1, color="red", linetype="dashed", size=0.8)+
  scale_y_continuous(breaks = c(round(min(Expr_ceres$ceres),digits = 0), -1,0,1, round(max(Expr_ceres$ceres),digits = 0)))
dev.off()