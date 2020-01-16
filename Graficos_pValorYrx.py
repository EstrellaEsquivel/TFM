#!/usr/bin/env python
import pandas as pd
import re 
import seaborn as sns
import os
import matplotlib.pyplot as plt
import math
import numpy as np
from matplotlib.pyplot import figure
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from statsmodels.sandbox.stats.multicomp import multipletests

####### Juntas los dos data frames  para los p-valores #############################

print('Juntando los dos data frames: locales y globales')

df_pVal_Mceres = pd.read_csv('/home/cns87/cns87433/thresholding/Th_upper_y_lower/df_pVal_Mceres.csv', sep = '\t')
df_pVal_Mceres_local = pd.read_csv('/home/cns87/cns87433/thresholding/con_th_locales/df_pVal_Mceres_locales.csv', sep = '\t')

df_todo = pd.concat([df_pVal_Mceres, df_pVal_Mceres_local])

p_valor_ajustado_B=multipletests(df_todo['p_valor'].tolist(), alpha=.05, method='bonferroni') #Correcion Bonferroni
df_todo['pValue_corregido_B'] = p_valor_ajustado_B[1]

p_valor_ajustado_FDR=multipletests(df_todo['p_valor'].tolist(), alpha=.05, method='fdr_bh') #Correccion Benjamini/Hochberg
df_todo['pValue_corregido_FDR'] = p_valor_ajustado_FDR[1]


df_todo['log_pValue_corregido_B'] = (np.log(df_todo['pValue_corregido_B'].astype('float64')))*(-1)
df_todo['log_pValue_corregido_FDR'] = (np.log(df_todo['pValue_corregido_FDR'].astype('float64')))*(-1)

df_todo.columns = ['Cell_Lines','Models', 'pValue', 'MeanCERES', 'MeanExpr', 'logPValue', 'p_adjusted_B', 'p_adjusted_FDR', 'LogP_adjusted_B', 'LogP_adjusted_FDR']


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####### Represento el heatmap clustering con los p-valores ################################


print('Represento el heatmap clustering con los p-valores')
#Seleccionas los valores que quieres para el heatmap
heatmap_data = pd.pivot_table(df_todo, values='LogP_adjusted_B', index=['Cell_Lines'], columns='Models')

#Cambias la columna de cell_line, por solo el nombre de su tejido.
index = heatmap_data.index.tolist()
tejidos = []
for cadena in index:
    s = '_'
    a= cadena.split('_')[1:]
    b= s.join(a)
    tejidos.append(b)       
    
    
heatmap_data = pd.pivot_table(df_todo, values='LogP_adjusted_B', index=['Cell_Lines'], columns='Models')
heatmap_data['Cell_Lines'] = tejidos
heatmap_data = heatmap_data.set_index('Cell_Lines')

# Preparas un vector de colores para mapearlo con los tejidos
color = sns.color_palette("bright", 25)
my_palette = dict(zip(heatmap_data.index.unique(), color))
#row_colors = heatmap_data.index.map(my_palette)
row_colors = [my_palette[x] for x in heatmap_data.index]

sns.set(font_scale=4)            
bar_color={'label': 'Exponent value of p Value adjusted by Bonferroni'}

#Significativo p_palor_corregido <= alpha | en log: log_p_valor_corregido <= log(alpha)
g = sns.clustermap(heatmap_data, cmap=sns.cm.rocket_r, row_colors=row_colors, xticklabels=True,
                   yticklabels=False, method='average', cbar_kws=bar_color,
                   mask=(heatmap_data<1.30), vmin=1.30) #La mascara es para que lo no significativos no tengan color

g.cax.set_position((0.94,0.19,0.03,0.5))

g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = 40)

for a in g.ax_col_dendrogram.collections: #Ancho de linea del dendrograma
    a.set_linewidth(7)

g.fig.set_size_inches((30,25))

g.savefig("./thresholding/Cluster_heat_Map_corregido_mascara.png", dpi=50)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####### Juntas los dos data frames  para las reacciones #############################

print('Juntas los dos data frames  para las reacciones: locales y globales')
df_Nrx_globales = pd.read_csv('./thresholding/Th_upper_y_lower/df_Nrx_globales.csv', sep = '\t')
df_Nrx_locales = pd.read_csv('./thresholding/con_th_locales/df_Nrx_locales.csv', sep = '\t')


df_Nrx = pd.concat([df_Nrx_globales, df_Nrx_locales])
df_Nrx.columns = ['cell_lines','Models', 'Rx_total', 'Genes_essential', 'SumCERES']
df_Nrx = df_Nrx.set_index('cell_lines')

df_Nrx['meanCERES'] = df_Nrx['SumCERES']/df_Nrx['Genes_essential']

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####### Represento el heatmap clustering con los valores del total de Rx ################################

print('Represento el heatmap clustering con los valores del total de Rx')
heatmap_data = pd.pivot_table(df_Nrx, values='Rx_total', index=['cell_lines'], columns='Models')

#Cambias la columna de cell_line, por solo el nombre de su tejido.
index = heatmap_data.index.tolist()
tejidos = []
for cadena in index:
    s = '_'
    a= cadena.split('_')[1:]
    b= s.join(a)
    tejidos.append(b)       
    
    
heatmap_data = pd.pivot_table(df_Nrx, values='Rx_total', index=['cell_lines'], columns='Models')
heatmap_data['cell_lines'] = tejidos
heatmap_data = heatmap_data.set_index('cell_lines')

# Preparas un vector de colores para mapearlo con los tejidos
color = sns.color_palette("bright", 25)
my_palette = dict(zip(heatmap_data.index.unique(), color))
#row_colors = heatmap_data.index.map(my_palette)
row_colors = [my_palette[x] for x in heatmap_data.index]
          
sns.set(font_scale=6)            
bar_color={'label': 'Number of total reaction'}
g = sns.clustermap(heatmap_data,  yticklabels=False, cmap=sns.cm.rocket_r, xticklabels=True,
                   method='average', row_colors=row_colors, cbar_kws=bar_color)
g.cax.set_position((0.94,0.19,0.03,0.5))
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = 40)

for a in g.ax_col_dendrogram.collections: #Ancho de linea del dendrograma
    a.set_linewidth(7)

g.fig.set_size_inches((40,40))
g.fig.set_size_inches((30,30))
g.savefig("./thresholding/HeatMap_Nrx.png", bbox_inches='tight', dpi=50)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Numero total de Rx vs Numero total de genes predichos como esenciales ######

print('Represento el Numero total de Rx vs Numero total de genes predichos como esenciales')
sns.set(font_scale=2)
g = sns.jointplot(x="Rx_total", y="Genes_essential", data=df_Nrx, alpha=0.15, color="royalblue") 

g.fig.set_size_inches(18, 12)
fig_ruta = './thresholding/Rxtotal_vs_predictGenes2.0.png'
g.savefig(fig_ruta, bbox_inches='tight')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Numero total de Rx vs Media CERES score de los genes predichos ######

print('Numero total de Rx vs Media CERES score de los genes predichos')
sns.set(font_scale=2)
g = sns.jointplot(x="Rx_total", y="meanCERES", data=df_Nrx, alpha=0.15, color="royalblue")

g.fig.set_size_inches(18, 12)
fig_ruta = './thresholding/Rxtotal_vs_SumCERES2.0.png'
g.savefig(fig_ruta, bbox_inches='tight')


