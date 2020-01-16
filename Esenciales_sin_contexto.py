#!/usr/bin/env python
from corda import reaction_confidence
from corda import test_model
from corda import CORDA
from cobra.io import read_sbml_model
import pandas as pd
import re  
import cobra.test
from cobra.flux_analysis import (single_gene_deletion)
from scipy.stats import mannwhitneyu
import seaborn as sns
import os
import matplotlib.pyplot as plt

#Cargando datos
print('--------CARGANDO DATOS----------')
sbml_fname = './Datos/models/Recon2.2.1_RPMI_trimed_gene_symbols.xml'
reference_model = read_sbml_model(sbml_fname)

df_ceres = pd.read_csv("./Datos/depmap/M19Q2_ceres_metabol_python.csv", sep = '\t')
df_ceres = df_ceres.set_index('cell_line')


# Identificando los genes esenciales, como los genes cuya deleccion individual no producen crecimiento (Ecuación biomasa < 0.01).
print('------SINGLE GENE DELETION-----')
lista_genes = reference_model.genes
resultado_knocked_out = single_gene_deletion(reference_model, lista_genes)

rename_dict = {i:list(i)[0] for i in resultado_knocked_out.index} #Renonbro la frozen set para sea una lista
df_deletion_renamed = resultado_knocked_out.rename(rename_dict, axis=0)

threshold = 0.01 * df_deletion_renamed.growth.max() #Este es el threshold para luego poder seleccionar los que esten por debajo de el.

mask = df_deletion_renamed.growth < threshold 
essential = df_deletion_renamed.index[mask]

mask2 = df_deletion_renamed.growth >= threshold
non_essential = df_deletion_renamed.index[mask2]

df_ceres2 = df_ceres.T
essential_in_ceres = set(essential) & set(df_ceres2.index) 
non_essential_in_ceres = set(non_essential) & set(df_ceres2.index)
total_rx = len(reference_model.reactions)

print('Numero de genes esenciales predichos:', str(len(essential)))
print('Numero de genes no esenciales predichos:', str(len(non_essential)))
print('Numero de reacciones en el modelo', str(total_rx))

#Guardadndo la lista de genes esenciales sin aplicar contexto
print('-----------------GUARDADNDO LISTA GENES----------------')
with open('Esenciales_sin_contexto.txt', 'w') as f:
    for gen in essential_in_ceres:
        f.write("%s\n" % gen)

