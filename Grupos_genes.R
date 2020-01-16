#####################################################################################################################################################
######################  ------------------- CLASIFICACION DE GENES SEGUN MEDIA Y VARIANZA-----------------------#####################################
#####################################################################################################################################################
########## GRUPO1 ######## media baja, varianza baja
library(dplyr)
grupo1_metabol_new <-data.frame(media_var_metabol%>%filter(media_gen<quantile(media_var_metabol$media_gen)[2] & 
                                                             varianza_gen<=quantile(media_var_metabol$varianza_gen)[4]))
colnames(grupo1_metabol_new)<- c('genes', 'varianza_gen', 'media_gen', 'sd_gen', 'CV_gen')
write.table(grupo1_metabol_new, "grupo1_metabol_new.txt", quote = F, sep = "\t", row.names = F, col.names = T)

#Seleccionamos los genes metabolicos del grupo1 en ceres
grupo1_ceres_new<-data.frame(media_var_ceres_metabol%>%filter(media_var_ceres_metabol$genes%in%grupo1_metabol_new$genes))

#Unimos la seleccion de genes con los datos de ceres para ver la esencialidad de los genes
grupo1_ceres_metabol_new <-merge(grupo1_metabol_new, grupo1_ceres_new, by='genes', all=T)
write.table(grupo1_ceres_metabol_new, "grupo1_ceres_metabol_new.txt", quote = F, sep = "\t", row.names = F, col.names = T)

########## GRUPO2 ######## Media alta, varianza baja
grupo2_metabol_new <-data.frame(media_var_metabol%>%filter(media_gen>=quantile(media_var_metabol$media_gen)[2]&media_gen<=quantile(media_var_metabol$media_gen)[3] 
                                                           & varianza_gen<=quantile(media_var_metabol$varianza_gen)[4]))
colnames(grupo2_metabol_new)<- c('genes', 'varianza_gen', 'media_gen', 'sd_gen', 'CV_gen')
write.table(grupo2_metabol_new, "grupo2_metabol_new.txt", quote = F, sep = "\t", row.names = F, col.names = T)

#Seleccionamos los genes metabolicos del grupo1 en ceres
grupo2_ceres_new<-data.frame(media_var_ceres_metabol%>%filter(media_var_ceres_metabol$genes%in%grupo2_metabol_new$genes))

#Unimos la seleccion de genes con los datos de ceres para ver la esencialidad de los genes
grupo2_ceres_metabol_new <-merge(grupo2_metabol_new, grupo2_ceres_new, by='genes', all=T)
write.table(grupo2_ceres_metabol_new, "grupo2_ceres_metabol_new.txt", quote = F, sep = "\t", row.names = F, col.names = T)

########## GRUPO3 ######## Media muy alta y varianzas baja
grupo3_metabol_new <-data.frame(media_var_metabol%>%filter(media_gen>=quantile(media_var_metabol$media_gen)[3]& media_gen<=quantile(media_var_metabol$media_gen)[4]&
                                                             varianza_gen<=quantile(media_var_metabol$varianza_gen)[3]))
colnames(grupo3_metabol_new)<- c('genes', 'varianza_gen', 'media_gen', 'sd_gen', 'CV_gen')
write.table(grupo3_metabol_new, "grupo3_metabol_new.txt", quote = F, sep = "\t", row.names = F, col.names = T)

#Seleccionamos los genes metabolicos del grupo1 en ceres
grupo3_ceres_new<-data.frame(media_var_ceres_metabol%>%filter(media_var_ceres_metabol$genes%in%grupo3_metabol_new$genes))

#Unimos la seleccion de genes con los datos de ceres para ver la esencialidad de los genes
grupo3_ceres_metabol_new <-merge(grupo3_metabol_new, grupo3_ceres_new, by='genes', all=T)
write.table(grupo3_ceres_metabol_new, "grupo3_ceres_metabol_new.txt", quote = F, sep = "\t", row.names = F, col.names = T)

########## GRUPO4 ######## Media baja y varianzas alta
grupo4_metabol_new <-data.frame(media_var_metabol%>%filter(media_gen>=quantile(media_var_metabol$media_gen)[4]
                                                           & varianza_gen<=quantile(media_var_metabol$varianza_gen)[3]))
colnames(grupo4_metabol_new)<- c('genes', 'varianza_gen', 'media_gen', 'sd_gen', 'CV_gen')
write.table(grupo4_metabol_new, "grupo4_metabol_new.txt", quote = F, sep = "\t", row.names = F, col.names = T)

#Seleccionamos los genes metabolicos del grupo1 en ceres
grupo4_ceres_new<-data.frame(media_var_ceres_metabol%>%filter(media_var_ceres_metabol$genes%in%grupo4_metabol_new$genes))

#Unimos la seleccion de genes con los datos de ceres para ver la esencialidad de los genes
grupo4_ceres_metabol_new <-merge(grupo4_metabol_new, grupo4_ceres_new, by='genes', all=T)
write.table(grupo4_ceres_metabol_new, "grupo4_ceres_metabol_new.txt", quote = F, sep = "\t", row.names = F, col.names = T)

########## GRUPO5 ######## Media muy alta y varianzas alta
grupo5_metabol_new <-data.frame(media_var_metabol%>%filter(media_gen<quantile(media_var_metabol$media_gen)[2] & 
                                                             varianza_gen>quantile(media_var_metabol$varianza_gen)[4]& varianza_gen<=upper_whisker))
colnames(grupo5_metabol_new)<- c('genes', 'varianza_gen', 'media_gen', 'sd_gen', 'CV_gen')
write.table(grupo5_metabol_new, "grupo5_metabol_new.txt", quote = F, sep = "\t", row.names = F, col.names = T)

#Seleccionamos los genes metabolicos del grupo1 en ceres
grupo5_ceres_new<-data.frame(media_var_ceres_metabol%>%filter(media_var_ceres_metabol$genes%in%grupo5_metabol_new$genes))

#Unimos la seleccion de genes con los datos de ceres para ver la esencialidad de los genes
grupo5_ceres_metabol_new <-merge(grupo5_metabol_new, grupo5_ceres_new, by='genes', all=T)
write.table(grupo5_ceres_metabol_new, "grupo5_ceres_metabol_new.txt", quote = F, sep = "\t", row.names = F, col.names = T)

########## GRUPO6 ######## Media baja y varianzas muy alta
grupo6_metabol_new <-data.frame(media_var_metabol%>%filter(media_gen>=quantile(media_var_metabol$media_gen)[2]&media_gen<=quantile(media_var_metabol$media_gen)[3]
                                                           & varianza_gen>quantile(media_var_metabol$varianza_gen)[4]& varianza_gen<=upper_whisker))
colnames(grupo6_metabol_new)<- c('genes', 'varianza_gen', 'media_gen', 'sd_gen', 'CV_gen')
write.table(grupo6_metabol_new, "grupo6_metabol_new.txt", quote = F, sep = "\t", row.names = F, col.names = T)

#Seleccionamos los genes metabolicos del grupo1 en ceres
grupo6_ceres_new<-data.frame(media_var_ceres_metabol%>%filter(media_var_ceres_metabol$genes%in%grupo6_metabol_new$genes))

#Unimos la seleccion de genes con los datos de ceres para ver la esencialidad de los genes
grupo6_ceres_metabol_new <-merge(grupo6_metabol_new, grupo6_ceres_new, by='genes', all=T)
write.table(grupo6_ceres_metabol_new, "grupo6_ceres_metabol_new.txt", quote = F, sep = "\t", row.names = F, col.names = T)

########## GRUPO7 ######## Media alta y varianzas muy alta
grupo7_metabol_new <-data.frame(media_var_metabol%>%filter(media_gen>=quantile(media_var_metabol$media_gen)[3]& media_gen<=quantile(media_var_metabol$media_gen)[4]
                                                           & varianza_gen>quantile(media_var_metabol$varianza_gen)[4]& varianza_gen<=upper_whisker))
colnames(grupo7_metabol_new)<- c('genes', 'varianza_gen', 'media_gen', 'sd_gen', 'CV_gen')
write.table(grupo7_metabol_new, "grupo7_metabol_new.txt", quote = F, sep = "\t", row.names = F, col.names = T)

#Seleccionamos los genes metabolicos del grupo1 en ceres
grupo7_ceres_new<-data.frame(media_var_ceres_metabol%>%filter(media_var_ceres_metabol$genes%in%grupo7_metabol_new$genes))

#Unimos la seleccion de genes con los datos de ceres para ver la esencialidad de los genes
grupo7_ceres_metabol_new <-merge(grupo7_metabol_new, grupo7_ceres_new, by='genes', all=T)
write.table(grupo7_ceres_metabol_new, "grupo7_ceres_metabol_new.txt", quote = F, sep = "\t", row.names = F, col.names = T)

########## GRUPO8 ######## Media baja y varianzas alta
grupo8_metabol_new <-data.frame(media_var_metabol%>%filter(media_gen>=quantile(media_var_metabol$media_gen)[4]
                                                           & varianza_gen>quantile(media_var_metabol$varianza_gen)[4]& varianza_gen<=upper_whisker))
colnames(grupo8_metabol_new)<- c('genes', 'varianza_gen', 'media_gen', 'sd_gen', 'CV_gen')
write.table(grupo8_metabol_new, "grupo8_metabol_new.txt", quote = F, sep = "\t", row.names = F, col.names = T)

#Seleccionamos los genes metabolicos del grupo1 en ceres
grupo8_ceres_new<-data.frame(media_var_ceres_metabol%>%filter(media_var_ceres_metabol$genes%in%grupo8_metabol_new$genes))

#Unimos la seleccion de genes con los datos de ceres para ver la esencialidad de los genes
grupo8_ceres_metabol_new <-merge(grupo8_metabol_new, grupo8_ceres_new, by='genes', all=T)
write.table(grupo8_ceres_metabol_new, "grupo8_ceres_metabol_new.txt", quote = F, sep = "\t", row.names = F, col.names = T)

########## GRUPO9 ######## Media muy alta y varianzas alta
grupo9_metabol_new <-data.frame(media_var_metabol%>%filter(media_gen<quantile(media_var_metabol$media_gen)[2] & 
                                                             varianza_gen>upper_whisker))
colnames(grupo9_metabol_new)<- c('genes', 'varianza_gen', 'media_gen', 'sd_gen', 'CV_gen')
write.table(grupo9_metabol_new, "grupo9_metabol_new.txt", quote = F, sep = "\t", row.names = F, col.names = T)

#Seleccionamos los genes metabolicos del grupo1 en ceres
grupo9_ceres_new<-data.frame(media_var_ceres_metabol%>%filter(media_var_ceres_metabol$genes%in%grupo9_metabol_new$genes))

#Unimos la seleccion de genes con los datos de ceres para ver la esencialidad de los genes
grupo9_ceres_metabol_new <-merge(grupo9_metabol_new, grupo9_ceres_new, by='genes', all=T)
write.table(grupo9_ceres_metabol_new, "grupo9_ceres_metabol_new.txt", quote = F, sep = "\t", row.names = F, col.names = T)

########## GRUPO10 ######## Media baja y varianzas muy alta
grupo10_metabol_new <-data.frame(media_var_metabol%>%filter(media_gen>=quantile(media_var_metabol$media_gen)[2]&media_gen<=quantile(media_var_metabol$media_gen)[3]
                                                           & varianza_gen>upper_whisker))
colnames(grupo10_metabol_new)<- c('genes', 'varianza_gen', 'media_gen', 'sd_gen', 'CV_gen')
write.table(grupo10_metabol_new, "grupo10_metabol_new.txt", quote = F, sep = "\t", row.names = F, col.names = T)

#Seleccionamos los genes metabolicos del grupo1 en ceres
grupo10_ceres_new<-data.frame(media_var_ceres_metabol%>%filter(media_var_ceres_metabol$genes%in%grupo10_metabol_new$genes))

#Unimos la seleccion de genes con los datos de ceres para ver la esencialidad de los genes
grupo10_ceres_metabol_new <-merge(grupo10_metabol_new, grupo10_ceres_new, by='genes', all=T)
write.table(grupo10_ceres_metabol_new, "grupo10_ceres_metabol_new.txt", quote = F, sep = "\t", row.names = F, col.names = T)

########## GRUPO11 ######## 
grupo11_metabol_new <-data.frame(media_var_metabol%>%filter(media_gen>=quantile(media_var_metabol$media_gen)[3]& media_gen<=quantile(media_var_metabol$media_gen)[4]
                                                          & varianza_gen>upper_whisker))
colnames(grupo11_metabol_new)<- c('genes', 'varianza_gen', 'media_gen', 'sd_gen', 'CV_gen')
write.table(grupo11_metabol_new, "grupo11_metabol_new.txt", quote = F, sep = "\t", row.names = F, col.names = T)

#Seleccionamos los genes metabolicos del grupo1 en ceres
grupo11_ceres_new<-data.frame(media_var_ceres_metabol%>%filter(media_var_ceres_metabol$genes%in%grupo11_metabol_new$genes))

#Unimos la seleccion de genes con los datos de ceres para ver la esencialidad de los genes
grupo11_ceres_metabol_new <-merge(grupo11_metabol_new, grupo11_ceres_new, by='genes', all=T)
write.table(grupo11_ceres_metabol_new, "grupo11_ceres_metabol_new.txt", quote = F, sep = "\t", row.names = F, col.names = T)

########## GRUPO12 ######## Media baja y varianzas alta
grupo12_metabol_new <-data.frame(media_var_metabol%>%filter(media_gen>=quantile(media_var_metabol$media_gen)[4]
                                                          & varianza_gen>upper_whisker))
colnames(grupo12_metabol_new)<- c('genes', 'varianza_gen', 'media_gen', 'sd_gen', 'CV_gen')
write.table(grupo12_metabol_new, "grupo12_metabol_new.txt", quote = F, sep = "\t", row.names = F, col.names = T)

#Seleccionamos los genes metabolicos del grupo1 en ceres
grupo12_ceres_new<-data.frame(media_var_ceres_metabol%>%filter(media_var_ceres_metabol$genes%in%grupo12_metabol_new$genes))

#Unimos la seleccion de genes con los datos de ceres para ver la esencialidad de los genes
grupo12_ceres_metabol_new <-merge(grupo12_metabol_new, grupo12_ceres_new, by='genes', all=T)
write.table(grupo12_ceres_metabol_new, "grupo12_ceres_metabol_new.txt", quote = F, sep = "\t", row.names = F, col.names = T)

save.image("./20190717_Grupos.RData")

#####################################################################################################################################################################
################################################# REPRESENTACION GRAFICA DE ALGUNOS GENES DE CADA GRUPO #############################################################
#####################################################################################################################################################################
#~~~~~~~~~~~~~ grupo 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
png("EjemGrupo1_new.png", width = 15, height = 15, units = "in", res = 110)
ggplot(MatrizMetabolismo2, aes(x=MatrizMetabolismo2$PLA2G10)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#666666") + 
  geom_density(fill = "#669900", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(MatrizMetabolismo2$PLA2G10)), color="black", linetype="dashed", size=2)+
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,face="bold"), 
        axis.text=element_text(size=35, color="black", face = "bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="log2 TPM", y="Densidad")+
  ggtitle("DISTRIBUCIÓN EXPRESIÓN GEN PLA2G10")
dev.off()

png("EjemGrupo1_Ceres_new.png", width = 15, height = 15, units = "in", res = 110) 
ggplot(data=Ceres_metabol, aes(x=Ceres_metabol$PLA2G10)) +
  geom_density(color="darkblue", fill="#6699FF", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=2) +
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,color="black", face="bold"), 
        axis.text=element_text(size=35, color="black", face="bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN ESENCIABILIDAD GEN PLA2G10")
dev.off()

png("EjemGrupo1_1_new.png", width = 24, height = 18, units = "in", res = 110)
ggplot(MatrizMetabolismo2, aes(x=MatrizMetabolismo2$SLC4A9)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#666666") + 
  geom_density(fill = "#669900", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(MatrizMetabolismo2$PLA2G10)), color="black", linetype="dashed", size=2)+
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,face="bold"), 
        axis.text=element_text(size=35, color="black", face = "bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="log2 TPM", y="Densidad")+
  ggtitle("DISTRIBUCIÓN EXPRESIÓN GEN SLC4A9")
dev.off()

png("EjemGrupo1_1Ceres_new.png", width = 24, height = 18, units = "in", res = 110)
ggplot(data=Ceres_metabol, aes(x=Ceres_metabol$SLC4A9)) +
  geom_density(color="darkblue", fill="#6699FF", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=2) +
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,color="black", face="bold"), 
        axis.text=element_text(size=35, color="black", face="bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN ESENCIABILIDAD GEN SLC4A9")
dev.off()

#~~~~~~~~~~~~~ grupo 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
png("EjemGrupo2_new.png", width = 15, height = 15, units = "in", res = 110)
ggplot(MatrizMetabolismo2, aes(x=MatrizMetabolismo2$SLC4A5)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#666666") + 
  geom_density(fill = "#669900", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(MatrizMetabolismo2$PLA2G10)), color="black", linetype="dashed", size=2)+
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,face="bold"), 
        axis.text=element_text(size=35, color="black", face = "bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="log2 TPM", y="Densidad")+
  ggtitle("DISTRIBUCIÓN EXPRESIÓN GEN SLC4A5")
dev.off()

png("EjemGrupo2_Ceres_new.png", width = 15, height = 15, units = "in", res = 110)
ggplot(data=Ceres_metabol, aes(x=Ceres_metabol$SLC4A5)) +
  geom_density(color="darkblue", fill="#6699FF", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=2) +
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,color="black", face="bold"), 
        axis.text=element_text(size=35, color="black", face="bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN ESENCIABILIDAD GEN SLC4A5")
dev.off()

png("EjemGrupo2_1_new.png", width = 15, height = 15, units = "in", res = 110)
ggplot(MatrizMetabolismo2, aes(x=MatrizMetabolismo2$HYAL3)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#666666") + 
  geom_density(fill = "#669900", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(MatrizMetabolismo2$PLA2G10)), color="black", linetype="dashed", size=2)+
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,face="bold"), 
        axis.text=element_text(size=35, color="black", face = "bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="log2 TPM", y="Densidad")+
  ggtitle("DISTRIBUCIÓN EXPRESIÓN GEN HYAL3")
dev.off()

png("EjemGrupo2_1Ceres_new.png", width = 15, height = 15, units = "in", res = 110)
ggplot(data=Ceres_metabol, aes(x=Ceres_metabol$HYAL3)) +
  geom_density(color="darkblue", fill="#6699FF", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=2) +
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,color="black", face="bold"), 
        axis.text=element_text(size=35, color="black", face="bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN ESENCIABILIDAD GEN HYAL3")
dev.off()

#~~~~~~~~~~~~~ grupo 3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
png("EjemGrupo3_new.png", width = 15, height = 15, units = "in", res = 110)
ggplot(MatrizMetabolismo2, aes(x=MatrizMetabolismo2$ATP6V1A)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#666666") + 
  geom_density(fill = "#669900", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(MatrizMetabolismo2$ATP6V1A)), color="black", linetype="dashed", size=2)+
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,face="bold"), 
        axis.text=element_text(size=35, color="black", face = "bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="log2 TPM", y="Densidad")+
  ggtitle("DISTRIBUCIÓN EXPRESIÓN GEN ATP6V1A")
dev.off()

png("EjemGrupo3_Ceres_new.png", width = 15, height = 15, units = "in", res = 110)
ggplot(data=Ceres_metabol, aes(x=Ceres_metabol$ATP6V1A)) +
  geom_density(color="darkblue", fill="#6699FF", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=2) +
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,color="black", face="bold"), 
        axis.text=element_text(size=35, color="black", face="bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN ESENCIABILIDAD GEN ATP6V1A")
dev.off()

png("EjemGrupo3_1_new.png", width = 15, height = 15, units = "in", res = 110)
ggplot(MatrizMetabolismo2, aes(x=MatrizMetabolismo2$ACAD8)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#666666") + 
  geom_density(fill = "#669900", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(MatrizMetabolismo2$PLA2G10)), color="black", linetype="dashed", size=2)+
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,face="bold"), 
        axis.text=element_text(size=35, color="black", face = "bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="log2 TPM", y="Densidad")+
  ggtitle("DISTRIBUCIÓN EXPRESIÓN GEN ACAD8")
dev.off()

png("EjemGrupo3_1Ceres_new.png", width = 15, height = 15, units = "in", res = 110)
ggplot(data=Ceres_metabol, aes(x=Ceres_metabol$ACAD8)) +
  geom_density(color="darkblue", fill="#6699FF", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=2) +
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,color="black", face="bold"), 
        axis.text=element_text(size=35, color="black", face="bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN ESENCIABILIDAD GEN ACAD8")
dev.off()

#~~~~~~~~~~~~~ grupo 4 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
png("EjemGrupo4_new.png", width = 15, height = 15, units = "in", res = 110)
ggplot(MatrizMetabolismo2, aes(x=MatrizMetabolismo2$ALDOA)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#666666") + 
  geom_density(fill = "#669900", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(MatrizMetabolismo2$ALDOA)), color="black", linetype="dashed", size=2)+
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,face="bold"), 
        axis.text=element_text(size=35, color="black", face = "bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="log2 TPM", y="Densidad")+
  ggtitle("DISTRIBUCIÓN EXPRESIÓN GEN ALDOA")
dev.off()

png("EjemGrupo4_Ceres_new.png", width = 15, height = 15, units = "in", res = 110)
ggplot(data=Ceres_metabol, aes(x=Ceres_metabol$ALDOA)) +
  geom_density(color="darkblue", fill="#6699FF", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=2) +
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,color="black", face="bold"), 
        axis.text=element_text(size=35, color="black", face="bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN ESENCIABILIDAD GEN ALDOA")
dev.off()

png("EjemGrupo4_1_new.png", width = 15, height = 15, units = "in", res = 110)
ggplot(MatrizMetabolismo2, aes(x=MatrizMetabolismo2$DGUOK)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#666666") + 
  geom_density(fill = "#669900", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(MatrizMetabolismo2$PLA2G10)), color="black", linetype="dashed", size=2)+
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,face="bold"), 
        axis.text=element_text(size=35, color="black", face = "bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="log2 TPM", y="Densidad")+
  ggtitle("DISTRIBUCIÓN EXPRESIÓN GEN DGUOK")
dev.off()

png("EjemGrupo4_1Ceres_new.png", width = 15, height = 15, units = "in", res = 110)
ggplot(data=Ceres_metabol, aes(x=Ceres_metabol$DGUOK)) +
  geom_density(color="darkblue", fill="#6699FF", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=2) +
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,color="black", face="bold"), 
        axis.text=element_text(size=35, color="black", face="bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN ESENCIABILIDAD GEN DGUOK")
dev.off()

#~~~~~~~~~~~~~ grupo 5 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
png("EjemGrupo5_new.png", width = 15, height = 15, units = "in", res = 110)
ggplot(MatrizMetabolismo2, aes(x=MatrizMetabolismo2$MPO)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#666666") + 
  geom_density(fill = "#669900", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(MatrizMetabolismo2$MPO)), color="black", linetype="dashed", size=2)+
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,face="bold"), 
        axis.text=element_text(size=35, color="black", face = "bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="log2 TPM", y="Densidad")+
  ggtitle("DISTRIBUCIÓN EXPRESIÓN GEN MPO")
dev.off()

png("EjemGrupo5_Ceres_new.png", width = 15, height = 15, units = "in", res = 110)
ggplot(data=Ceres_metabol, aes(x=Ceres_metabol$MPO)) +
  geom_density(color="darkblue", fill="#6699FF", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=2) +
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,color="black", face="bold"), 
        axis.text=element_text(size=35, color="black", face="bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN ESENCIABILIDAD GEN MPO")
dev.off()

png("EjemGrupo5_1_new.png", width = 15, height = 15, units = "in", res = 110)
ggplot(MatrizMetabolismo2, aes(x=MatrizMetabolismo2$FABP1)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#666666") + 
  geom_density(fill = "#669900", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(MatrizMetabolismo2$PLA2G10)), color="black", linetype="dashed", size=2)+
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,face="bold"), 
        axis.text=element_text(size=35, color="black", face = "bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="log2 TPM", y="Densidad")+
  ggtitle("DISTRIBUCIÓN EXPRESIÓN GEN FABP1")
dev.off()

png("EjemGrupo5_1Ceres_new.png", width = 15, height = 15, units = "in", res = 110)
ggplot(data=Ceres_metabol, aes(x=Ceres_metabol$FABP1)) +
  geom_density(color="darkblue", fill="#6699FF", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=2) +
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,color="black", face="bold"), 
        axis.text=element_text(size=35, color="black", face="bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN ESENCIABILIDAD GEN FABP1")
dev.off()

#~~~~~~~~~~~~~ grupo 6 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
png("EjemGrupo6_new.png", width = 15, height = 15, units = "in", res = 110)
ggplot(MatrizMetabolismo2, aes(x=MatrizMetabolismo2$HSD11B2)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#666666") + 
  geom_density(fill = "#669900", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(MatrizMetabolismo2$PLA2G10)), color="black", linetype="dashed", size=2)+
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,face="bold"), 
        axis.text=element_text(size=35, color="black", face = "bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="log2 TPM", y="Densidad")+
  ggtitle("DISTRIBUCIÓN EXPRESIÓN GEN HSD11B2")
dev.off()

png("EjemGrupo6_Ceres_new.png", width = 15, height = 15, units = "in", res = 110)
ggplot(data=Ceres_metabol, aes(x=Ceres_metabol$HSD11B2)) +
  geom_density(color="darkblue", fill="#6699FF", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=2) +
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,color="black", face="bold"), 
        axis.text=element_text(size=35, color="black", face="bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN ESENCIABILIDAD GEN HSD11B2")
dev.off()

png("EjemGrupo6_1_new.png", width = 15, height = 15, units = "in", res = 110)
ggplot(MatrizMetabolismo2, aes(x=MatrizMetabolismo2$A4GALT)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#666666") + 
  geom_density(fill = "#669900", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(MatrizMetabolismo2$A4GALT)), color="black", linetype="dashed", size=2)+
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,face="bold"), 
        axis.text=element_text(size=35, color="black", face = "bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="log2 TPM", y="Densidad")+
  ggtitle("DISTRIBUCIÓN EXPRESIÓN GEN A4GALT")
dev.off()

png("EjemGrupo6_1Ceres_new.png", width = 15, height = 15, units = "in", res = 110)
ggplot(data=Ceres_metabol, aes(x=Ceres_metabol$A4GALT)) +
  geom_density(color="darkblue", fill="#6699FF", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=2) +
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,color="black", face="bold"), 
        axis.text=element_text(size=35, color="black", face="bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN ESENCIABILIDAD GEN A4GALT")
dev.off()

#~~~~~~~~~~~~~ grupo 7 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
png("EjemGrupo7_new.png", width = 15, height = 15, units = "in", res = 110)
ggplot(MatrizMetabolismo2, aes(x=MatrizMetabolismo2$TM7SF2)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#666666") + 
  geom_density(fill = "#669900", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(MatrizMetabolismo2$TM7SF2)), color="black", linetype="dashed", size=2)+
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,face="bold"), 
        axis.text=element_text(size=35, color="black", face = "bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="log2 TPM", y="Densidad")+
  ggtitle("DISTRIBUCIÓN EXPRESIÓN GEN TM7SF2")
dev.off()

png("EjemGrupo7_Ceres_new.png", width = 15, height = 15, units = "in", res = 110)
ggplot(data=Ceres_metabol, aes(x=Ceres_metabol$TM7SF2)) +
  geom_density(color="darkblue", fill="#6699FF", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=2) +
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,color="black", face="bold"), 
        axis.text=element_text(size=35, color="black", face="bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN ESENCIABILIDAD GEN TM7SF2")
dev.off()

png("EjemGrupo7_1_new.png", width = 15, height = 15, units = "in", res = 110)
ggplot(MatrizMetabolismo2, aes(x=MatrizMetabolismo2$PIK3R3)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#666666") + 
  geom_density(fill = "#669900", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(MatrizMetabolismo2$PIK3R3)), color="black", linetype="dashed", size=2)+
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,face="bold"), 
        axis.text=element_text(size=35, color="black", face = "bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="log2 TPM", y="Densidad")+
  ggtitle("DISTRIBUCIÓN EXPRESIÓN GEN PIK3R3")
dev.off()

png("EjemGrupo7_1Ceres_new.png", width = 15, height = 15, units = "in", res = 110)
ggplot(data=Ceres_metabol, aes(x=Ceres_metabol$PIK3R3)) +
  geom_density(color="darkblue", fill="#6699FF", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=2) +
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,color="black", face="bold"), 
        axis.text=element_text(size=35, color="black", face="bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN ESENCIABILIDAD GEN PIK3R3")
dev.off()

#~~~~~~~~~~~~~ grupo 8 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
png("EjemGrupo8_new.png", width = 15, height = 15, units = "in", res = 110)
ggplot(MatrizMetabolismo2, aes(x=MatrizMetabolismo2$TXNRD1)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#666666") + 
  geom_density(fill = "#669900", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(MatrizMetabolismo2$TXNRD1)), color="black", linetype="dashed", size=2)+
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,face="bold"), 
        axis.text=element_text(size=35, color="black", face = "bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="log2 TPM", y="Densidad")+
  ggtitle("DISTRIBUCIÓN EXPRESIÓN GEN TXNRD1")
dev.off()

png("EjemGrupo8_Ceres_new.png", width = 15, height = 15, units = "in", res = 110)
ggplot(data=Ceres_metabol, aes(x=Ceres_metabol$TXNRD1)) +
  geom_density(color="darkblue", fill="#6699FF", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=2) +
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,color="black", face="bold"), 
        axis.text=element_text(size=35, color="black", face="bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN ESENCIABILIDAD GEN TXNRD1")
dev.off()

png("EjemGrupo8_1_new.png", width = 15, height = 15, units = "in", res = 110)
ggplot(MatrizMetabolismo2, aes(x=MatrizMetabolismo2$SCD)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#666666") + 
  geom_density(fill = "#669900", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(MatrizMetabolismo2$SCD)), color="black", linetype="dashed", size=2)+
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,face="bold"), 
        axis.text=element_text(size=35, color="black", face = "bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="log2 TPM", y="Densidad")+
  ggtitle("DISTRIBUCIÓN EXPRESIÓN GEN SCD")
dev.off()

png("EjemGrupo8_1Ceres_new.png", width = 15, height = 15, units = "in", res = 110)
ggplot(data=Ceres_metabol, aes(x=Ceres_metabol$SCD)) +
  geom_density(color="darkblue", fill="#6699FF", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=2) +
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,color="black", face="bold"), 
        axis.text=element_text(size=35, color="black", face="bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN ESENCIABILIDAD GEN SCD")
dev.off()
#~~~~~~~~~~~~~ grupo 9 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Ningun dato cae en ese grupo
#~~~~~~~~~~~~~ grupo 10 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
png("EjemGrupo10_new.png", width = 18, height = 18, units = "in", res = 110)
ggplot(MatrizMetabolismo2, aes(x=MatrizMetabolismo2$GFPT2)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#666666") + 
  geom_density(fill = "#669900", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(MatrizMetabolismo2$GFPT2)), color="black", linetype="dashed", size=2)+
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,face="bold"), 
        axis.text=element_text(size=35, color="black", face = "bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="log2 TPM", y="Densidad")+
  ggtitle("DISTRIBUCIÓN EXPRESIÓN GEN GFPT2")
dev.off()

png("EjemGrupo10_Ceres_new.png", width = 15, height = 15, units = "in", res = 110)
ggplot(data=Ceres_metabol, aes(x=Ceres_metabol$GFPT2)) +
  geom_density(color="darkblue", fill="#6699FF", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=2) +
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,color="black", face="bold"), 
        axis.text=element_text(size=35, color="black", face="bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN ESENCIABILIDAD GEN GFPT2")
dev.off()

png("EjemGrupo10_1_new.png", width = 15, height = 15, units = "in", res = 110)
ggplot(MatrizMetabolismo2, aes(x=MatrizMetabolismo2$ALDH3A1)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#666666") + 
  geom_density(fill = "#669900", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(MatrizMetabolismo2$ALDH3A1)), color="black", linetype="dashed", size=2)+
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,face="bold"), 
        axis.text=element_text(size=35, color="black", face = "bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="log2 TPM", y="Densidad")+
  ggtitle("DISTRIBUCIÓN EXPRESIÓN GEN ALDH3A1")
dev.off()

png("EjemGrupo10_1Ceres_new.png", width = 15, height = 15, units = "in", res = 110)
ggplot(data=Ceres_metabol, aes(x=Ceres_metabol$ALDH3A1)) +
  geom_density(color="darkblue", fill="#6699FF", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=2) +
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,color="black", face="bold"), 
        axis.text=element_text(size=35, color="black", face="bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN ESENCIABILIDAD GEN ALDH3A1")
dev.off()

#~~~~~~~~~~~~~ grupo 11 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
png("EjemGrupo11_new.png", width = 15, height = 15, units = "in", res = 110)
ggplot(MatrizMetabolismo2, aes(x=MatrizMetabolismo2$GGT1)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#666666") + 
  geom_density(fill = "#669900", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(MatrizMetabolismo2$GGT1)), color="black", linetype="dashed", size=2)+
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,face="bold"), 
        axis.text=element_text(size=35, color="black", face = "bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="log2 TPM", y="Densidad")+
  ggtitle("DISTRIBUCIÓN EXPRESIÓN GEN GGT1")
dev.off()

png("EjemGrupo11_Ceres_new.png", width = 15, height = 15, units = "in", res = 110)
ggplot(data=Ceres_metabol, aes(x=Ceres_metabol$GGT1)) +
  geom_density(color="darkblue", fill="#6699FF", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=2) +
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,color="black", face="bold"), 
        axis.text=element_text(size=35, color="black", face="bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN ESENCIABILIDAD GEN GGT1")
dev.off()

png("EjemGrupo11_1_new.png", width = 15, height = 15, units = "in", res = 110)
ggplot(MatrizMetabolismo2, aes(x=MatrizMetabolismo2$ABCC3)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#666666") + 
  geom_density(fill = "#669900", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(MatrizMetabolismo2$ABCC3)), color="black", linetype="dashed", size=2)+
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,face="bold"), 
        axis.text=element_text(size=35, color="black", face = "bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="log2 TPM", y="Densidad")+
  ggtitle("DISTRIBUCIÓN EXPRESIÓN GEN ABCC3")
dev.off()

png("EjemGrupo11_1Ceres_new.png", width = 15, height = 15, units = "in", res = 110)
ggplot(data=Ceres_metabol, aes(x=Ceres_metabol$ABCC3)) +
  geom_density(color="darkblue", fill="#6699FF", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=2) +
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,color="black", face="bold"), 
        axis.text=element_text(size=35, color="black", face="bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN ESENCIABILIDAD GEN ABCC3")
dev.off()

#~~~~~~~~~~~~~ grupo 12 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
png("EjemGrupo12_new.png", width = 18, height = 18, units = "in", res = 110)
ggplot(MatrizMetabolismo2, aes(x=MatrizMetabolismo2$ASS1)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#666666") + 
  geom_density(fill = "#669900", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(MatrizMetabolismo2$ASS1)), color="black", linetype="dashed", size=2)+
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,face="bold"), 
        axis.text=element_text(size=35, color="black", face = "bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="log2 TPM", y="Densidad")+
  ggtitle("DISTRIBUCIÓN EXPRESIÓN GEN ASS1")
dev.off()

png("EjemGrupo12_Ceres_new.png", width = 15, height = 15, units = "in", res = 110)
ggplot(data=Ceres_metabol, aes(x=Ceres_metabol$ASS1)) +
  geom_density(color="darkblue", fill="#6699FF", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=2) +
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,color="black", face="bold"), 
        axis.text=element_text(size=35, color="black", face="bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN ESENCIABILIDAD GEN ASS1")
dev.off()

png("EjemGrupo12_1_new.png", width = 15, height = 15, units = "in", res = 110)
ggplot(MatrizMetabolismo2, aes(x=MatrizMetabolismo2$PYGL)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "#666666") + 
  geom_density(fill = "#669900", alpha = 0.3) + 
  geom_vline(aes(xintercept=mean(MatrizMetabolismo2$PYGL)), color="black", linetype="dashed", size=2)+
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,face="bold"), 
        axis.text=element_text(size=35, color="black", face = "bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="log2 TPM", y="Densidad")+
  ggtitle("DISTRIBUCIÓN EXPRESIÓN GEN PYGL")
dev.off()

png("EjemGrupo12_1Ceres_new.png", width = 15, height = 15, units = "in", res = 110)
ggplot(data=Ceres_metabol, aes(x=Ceres_metabol$PYGL)) +
  geom_density(color="darkblue", fill="#6699FF", alpha=.5)+
  geom_vline(aes(xintercept=(-1)), color="red", linetype="dashed", size=2) +
  theme(panel.background = element_rect(fill='white', colour='white'),
        axis.title=element_text(size=35,color="black", face="bold"), 
        axis.text=element_text(size=35, color="black", face="bold"), 
        plot.title = element_text(size=20, face="bold"))+
  labs(x="Valores de CERES", Y="Densidad")+
  ggtitle("DISTRIBUCIÓN ESENCIABILIDAD GEN PYGL")
dev.off()

