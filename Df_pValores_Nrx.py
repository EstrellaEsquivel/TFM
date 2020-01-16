#!/usr/bin/env python
import pandas as pd
import re  
from scipy.stats import mannwhitneyu
import seaborn as sns
import os
import matplotlib.pyplot as plt
import math
import numpy as np
from matplotlib.pyplot import figure
from statsmodels.sandbox.stats.multicomp import multipletests

df_expresion = pd.read_csv('/home/cns87/cns87433/Para_Mare/M19Q2_metabol_modelo.csv', sep = '\t')
df_expresion = df_expresion.set_index('cell_line')
df_ceres = pd.read_csv('/home/cns87/cns87433/Para_Mare/M19Q2_ceres_metabol_python.csv', sep = '\t')
df_ceres = df_ceres.set_index('cell_line')
lista_cell_lines = set(df_expresion.index) & set(df_ceres.index)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### Generas el DF de los p-valores con th globales ########
df_pVal_Mceres = pd.DataFrame(columns=['cell_line', 'Modelo', 'p_valor', 'M_ceres', 'M_expr'], index=None)
lista_modelos = ['P05_P90', 'P05_P80', 'P05_P70', 'P05_P65', 'P10_P90', 'P10_P80', 'P10_P70', 'P10_P65', 'P15_P90', 'P15_P80', 'P15_P70', 'P15_P65']


lista_modelo = []
lista_p_valores = []
lista_mean_Sceres = []
lista_mean_Expr = []
lista_cell_line = []

print('Generando el data frame con p-valor, Mean_Ceres, Mean_Exp, -log_p_valor')
for cell_line in lista_cell_lines:
    fold_path = './thresholding/Th_upper_y_lower/' + cell_line
    if not os.path.exists(fold_path):
        continue
    else:
        fname = fold_path + '/p_valor.txt'
        print(cell_line)
        with open(fname, 'r') as file:
            for line in file:
                a = line.split('\t')[0] #dentro de la primera fila el primer elemento, modelo
                b = line.split('\t')[1] # p-valor
                c = line.split('\t')[2] # ceres
                d = line.split('\t')[3] #expr
                e = cell_line
                lista_modelo.append(a)
                lista_p_valores.append(b)
                lista_mean_Sceres.append(c)
                lista_mean_Expr.append(d)
                lista_cell_line.append(e)

df_pVal_Mceres['Modelo'] = lista_modelo
df_pVal_Mceres['p_valor'] = lista_p_valores
df_pVal_Mceres['M_ceres'] = lista_mean_Sceres
df_pVal_Mceres['M_expr'] = lista_mean_Expr
df_pVal_Mceres['cell_line'] = lista_cell_line


df_pVal_Mceres['M_expr'].replace(regex=True,inplace=True,to_replace=r'\n',value=r'')
df_pVal_Mceres['log_p_valor'] = (np.log(df_pVal_Mceres['p_valor'].astype('float64')))*(-1)
df_pVal_Mceres.to_csv('./thresholding/Th_upper_y_lower/df_pVal_Mceres.csv', sep = '\t', index=None, header=True)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### Generas el DF con los p-valores con th locales ########
df_pVal_Mceres_locales = pd.DataFrame(columns=['cell_line', 'Modelo', 'p_valor', 'M_ceres', 'M_expr'], index=None)


print('Generando el data frame con p-valor,Mean_Ceres, Mean_Exp, -log_p_valor')
lista_modelo = []
lista_p_valores = []
lista_mean_Sceres = []
lista_mean_Expr = []
lista_cell_line = []

for cell_line in lista_cell_lines:
    folder_path = '/home/cns87/cns87433/thresholding/con_th_locales/' + str(cell_line)
    if os.path.exists(folder_path):
        fname = folder_path +'/p_valor.txt'
        with open(fname, 'r') as file:
            for line in file:
                a = line.split('\t')[0] #dentro de la primera fila el primer elemento, modelo
                b = line.split('\t')[1] # p-valor
                c = line.split('\t')[2] # ceres
                d = line.split('\t')[3] # expr
                e = cell_line
                lista_modelo.append(a)
                lista_p_valores.append(b)
                lista_mean_Sceres.append(c)
                lista_mean_Expr.append(d)
                lista_cell_line.append(e)
                print(cell_line)
    else:
        continue

df_pVal_Mceres_locales['Modelo'] = lista_modelo
df_pVal_Mceres_locales['p_valor'] = lista_p_valores
df_pVal_Mceres_locales['M_ceres'] = lista_mean_Sceres
df_pVal_Mceres_locales['M_expr'] = lista_mean_Expr
df_pVal_Mceres_locales['cell_line'] = lista_cell_line
    
    
df_pVal_Mceres_locales['M_expr'].replace(regex=True,inplace=True,to_replace=r'\n',value=r'')
df_pVal_Mceres_locales['log_p_valor'] = (np.log(df_pVal_Mceres_locales['p_valor'].astype('float64')))*(-1)
   
    
df_pVal_Mceres_locales['M_expr'].replace(regex=True,inplace=True,to_replace=r'\n',value=r'')
df_pVal_Mceres_locales['log_p_valor'] = (np.log(df_pVal_Mceres_locales['p_valor'].astype('float64')))*(-1)
df_pVal_Mceres_locales.to_csv('./thresholding/con_th_locales/df_pVal_Mceres_locales.csv', sep = '\t', index=None, header=True)

#~ #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ ##### Generas el DF con el numero de Rx con th globales ########
df_Nrx_globales = pd.DataFrame(columns=['cell_line', 'Modelo', 'Rx_totales', 'Rx_essential', 'Suma_CERES'], index=None)

lista_modelo = []
lista_Rx_totales = []
lista_Rx_essential = []
lista_cell_line = []
lista_suma = []

for cell_line in lista_cell_lines:
    fold_path = './thresholding/Th_upper_y_lower/' + cell_line
    fname = fold_path + '/Number_reactions.txt'
    print(cell_line)
    with open(fname, 'r') as file:
        for line in file:
            a = line.split('\t')[0] #dentro de la primera fila el primer elemento, modelo
            b = line.split('\t')[1] # nRX
            c = line.split('\t')[2] # nEssential
            d = line.split('\t')[3] # sumaCeresPredichosEssential
            e = cell_line
            lista_modelo.append(a)
            lista_Rx_totales.append(b)
            lista_Rx_essential.append(c)
            lista_suma.append(d)
            lista_cell_line.append(e)

df_Nrx_globales['Modelo'] = lista_modelo
df_Nrx_globales['Rx_totales'] = lista_Rx_totales
df_Nrx_globales['Rx_essential'] = lista_Rx_essential
df_Nrx_globales['cell_line'] = lista_cell_line
df_Nrx_globales['Suma_CERES'] = lista_suma

df_Nrx_globales.to_csv('./thresholding/Th_upper_y_lower/df_Nrx_globales.csv', sep = '\t', index=None, header=True)

#~ #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ ##### Generas el DF con el numero de Rx con th locales ########

df_Nrx_locales = pd.DataFrame(columns=['cell_line', 'Modelo', 'Rx_totales', 'Rx_essential', 'Suma_CERES'], index=None)

lista_modelo = []
lista_Rx_totales = []
lista_Rx_essential = []
lista_cell_line = []
lista_suma = []

for cell_line in lista_cell_lines:
    fold_path = './thresholding/con_th_locales/' + cell_line
    fname = fold_path + '/Number_reactions.txt'
    print(cell_line)
    with open(fname, 'r') as file:
        for line in file:
            a = line.split('\t')[0] #dentro de la primera fila el primer elemento, modelo
            b = line.split('\t')[1] # nRX
            c = line.split('\t')[2] # nEssential
            d = line.split('\t')[3] # sumaCeresPredichosEssential
            e = cell_line
            lista_modelo.append(a)
            lista_Rx_totales.append(b)
            lista_Rx_essential.append(c)
            lista_suma.append(d)
            lista_cell_line.append(e)

df_Nrx_locales['Modelo'] = lista_modelo
df_Nrx_locales['Rx_totales'] = lista_Rx_totales
df_Nrx_locales['Rx_essential'] = lista_Rx_essential
df_Nrx_locales['cell_line'] = lista_cell_line
df_Nrx_locales['Suma_CERES'] = lista_suma

df_Nrx_locales.to_csv('./thresholding/con_th_locales/df_Nrx_locales.csv', sep = '\t', index=None, header=True)

