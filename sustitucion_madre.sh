#!/bin/bash
#Coges solo la columna de las lineas celulares
cut -f2 -d $'\t' /home/estrella/TFM/Datos/depmap/M19Q2_ceres_metabol.csv | sed '1d' | sed 's/"//g' > /home/estrella/TFM/lista_cell_line.txt

#Carpeta donde se van a guardar el script para cada linea
mkdir /home/estrella/TFM/Para_LaPalma/script_per_line_Lower_locales

for cell_line in $(cat /home/estrella/TFM/lista_cell_line.txt);
do
sed "s/Clave/${cell_line}/g" /home/estrella/TFM/Para_LaPalma/Combinando_LowerGlobal_locales_LaPalma.py > /home/estrella/TFM/Para_LaPalma/script_per_line_Lower_locales/modelo_${cell_line}.py
done
