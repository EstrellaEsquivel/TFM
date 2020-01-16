from corda import reaction_confidence
from corda import test_model
from corda import CORDA
from cobra.io import read_sbml_model
from cobra.io import write_sbml_model
import pandas as pd
from csm4cobra.context_specific_deletions import reaction_confidence
import re  
import cobra.test
from cobra.flux_analysis import (single_gene_deletion)
from scipy.stats import mannwhitneyu
import seaborn as sns
import os
import matplotlib.pyplot as plt
plt.switch_backend('agg')

print('--------CARGANDO DATOS----------')
sbml_fname = '/home/cns87/cns87433/Para_Mare/Recon2.2.1_RPMI_trimed_gene_symbols.xml'
reference_model = read_sbml_model(sbml_fname)
df_expresion = pd.read_csv('/home/cns87/cns87433/Para_Mare/M19Q2_metabol_modelo.csv', sep = '\t')
df_expresion = df_expresion.set_index('cell_line')
df_gpr = pd.read_csv("/home/cns87/cns87433/Para_Mare/Merge_GP_expression.csv",sep='\t')
df_gpr = df_gpr.set_index('reaction_id')
df_ceres = pd.read_csv('/home/cns87/cns87433/Para_Mare/M19Q2_ceres_metabol_python.csv', sep = '\t')
df_ceres = df_ceres.set_index('cell_line')


print('--------GENE MAPPING-------')
lista_cell_lines = ['Clave']
#lista_cell_lines = ['U118MG_CENTRAL_NERVOUS_SYSTEM']

for cell_line in lista_cell_lines:
    rxns_mapped_expr = {}
    df_ceres_t = df_ceres.T
    ceres_genes_linea = df_ceres_t[cell_line]
    lista_modelos = ['P05_P90_LocalMEDIA','P05_P80_LocalMEDIA', 'P05_P70_LocalMEDIA', 'P05_P65_LocalMEDIA', 'P10_P90_LocalMEDIA', 'P10_P80_LocalMEDIA', 'P10_P70_LocalMEDIA', 
                     'P10_P65_LocalMEDIA', 'P15_P90_LocalMEDIA', 'P15_P80_LocalMEDIA', 'P15_P70_LocalMEDIA', 'P15_P65_LocalMEDIA', 'Frecuencia']
    # df_frecuencia_esenciales = Df con todos los modelos y genes. Los genes esenciales en cada modelo es 1 y los no esenciales son 0. 
    df_frecuencia_esenciales = pd.DataFrame(index=ceres_genes_linea.index, columns=lista_modelos)


    for rxn in reference_model.reactions:
        gene_conf = {}
        for g in rxn.genes:
            if g.id in df_expresion.columns: 
                gene_conf[g.id] = df_expresion.loc[cell_line, g.id] 
            else:
                gene_conf[g.id] = -1
        rxns_mapped_expr[rxn.id] = reaction_confidence(rxn.gene_reaction_rule, gene_conf)

    print('---------DANDO CONFIDENCIAS-------')
    df_gpr_media = df_gpr.mean(axis=1)
    lista_th_lower = {'P05': df_gpr_media.quantile(0.05)}#, 'P15': df_gpr_media.quantile(0.15),'P10': df_gpr_media.quantile(0.10)}
    lista_th_upper = {'P90': df_gpr_media.quantile(0.90)}#, 'P65': df_gpr_media.quantile(0.65),'P80': df_gpr_media.quantile(0.80), 'P70': df_gpr_media.quantile(0.70)}
    


    folder_output = './thresholding/con_th_locales2/' + cell_line
    if not os.path.exists(folder_output):
        os.mkdir(folder_output)

        
    fname = folder_output + '/p_valor_locales.txt'
    fname2 = folder_output + '/Number_reactions_locales.txt'
    with open(fname, 'w') as p_valor:
        with open(fname2, 'w') as Number_reactions:
            for th_l in lista_th_lower:
                for th_u in lista_th_upper:
                    # Dando confidencias a los Rxs [high confidence = 3, Medium confidence = 2, Low confidence = 1, Not_expressed = -1, Unknow = 0]
                    rxns_conf = {}
                    for r in reference_model.reactions:
                        if r.id in rxns_mapped_expr and r.id in df_gpr_media.index:
                            #Para evitar un error. La rxn debe estar en la lista de rx mapeadas y en la matriz de las medias.
                            if rxns_mapped_expr[r.id] >= lista_th_upper[th_u] and df_gpr_media[r.id] >= lista_th_upper[th_u]:
                                if rxns_mapped_expr[r.id] >= df_gpr_media[r.id]:
                                    rxns_conf[r.id] = 3
                                else:
                                    rxns_conf[r.id] = 2
                            elif rxns_mapped_expr[r.id] >= lista_th_upper[th_u] and df_gpr_media[r.id] < lista_th_upper[th_u]:
                                rxns_conf[r.id] = 3


                            elif rxns_mapped_expr[r.id] <= lista_th_lower[th_l] and df_gpr_media[r.id] <= lista_th_lower[th_l]:
                                if rxns_mapped_expr[r.id] >= df_gpr_media[r.id]:
                                    rxns_conf[r.id] = 1
                                else:
                                    rxns_conf[r.id] = -1
                            elif rxns_mapped_expr[r.id] <= lista_th_lower[th_l] and df_gpr_media[r.id] > lista_th_lower[th_l]:
                                rxns_conf[r.id] = -1


                            elif lista_th_lower[th_l] <= rxns_mapped_expr[r.id] <= lista_th_upper[th_u] and df_gpr_media[r.id] >= lista_th_upper[th_u]:
                                rxns_conf[r.id] = 1
                            elif lista_th_lower[th_l] <= rxns_mapped_expr[r.id] <= lista_th_upper[th_u] and df_gpr_media[r.id] <= lista_th_lower[th_l]:
                                rxns_conf[r.id] = 2
                            elif lista_th_lower[th_l] <= rxns_mapped_expr[r.id] <= lista_th_upper[th_u] and lista_th_lower[th_l] <= df_gpr_media[r.id] <= lista_th_upper[th_u]:
                                if rxns_mapped_expr[r.id] > df_gpr_media[r.id]:
                                    rxns_conf[r.id] = 2
                                else:
                                    rxns_conf[r.id] = 1


                            elif rxns_mapped_expr[r.id] == -1:
                                rxns_conf[r.id] = 0
                        else:
                            rxns_conf[r.id] = 0
                    rxns_conf["biomass_reaction"] = 3  # Esta el 'objetive funtion' del FBA.

                    print('Total de reacciones del modelo de referencia Recon2.2:', str(len(reference_model.reactions)))

                    print('--------RECONSTRUYENDO MODELO--------')
                    CORDA_builder = CORDA(reference_model, rxns_conf)
                    CORDA_builder.build()
                                    
                    print('Resumen de CORDA:/n'+ str(CORDA_builder))
                    csm2 = CORDA_builder.cobra_model()
                    write_sbml_model(csm2, './thresholding/con_th_locales2/' + cell_line + '/csm2_%s_%s_LocalMEDIA.sbml' % (th_l, th_u))
                    print('Numero de reacciones de modelo reconstruido:', str(len(csm2.reactions)))

                    # Identificando los genes esenciales, como los genes cuya deleccion individual no producen crecimiento (Ecuación biomasa < 0.01).
                    print('------SINGLE GENE DELETION-----')
                    resultado_knocked_out = single_gene_deletion(csm2, proprocesses=1)
                    rename_dict = {i:list(i)[0] for i in resultado_knocked_out.index}

                    print('Seleccionando los genes esenciales y no esenciales')
                    df_deletion_renamed = resultado_knocked_out.rename(rename_dict, axis=0)

                    threshold = 0.01 * df_deletion_renamed.growth.max()

                    mask = df_deletion_renamed.growth < threshold 
                    essential = df_deletion_renamed.index[mask]

                    mask2 = df_deletion_renamed.growth >= threshold
                    non_essential = df_deletion_renamed.index[mask2]

                    print('Numero de genes esenciales predichos:', str(len(essential)))
                    print('Numero de genes no esenciales predichos:', str(len(non_essential)))
                    
                    total_rx = len(csm2.reactions)
                    essential = len(essential)
                    
                    line = th_l + '_' + th_u+'_LocalMEDIA'+ '\t' + str(total_rx) + '\t'+ str(essential) + '\t' + str(valor_ceres_suma)+ '\n'
                    Number_reactions.write(line)
                    

                    # Compruebas cuantos genes esenciales coinciden con CERES.
                    print('Comprobando que genes predichos coinciden con los de CERES')
                    df_ceres_t = df_ceres.T
                    ceres_genes_linea = df_ceres_t[cell_line]
                    essential_in_ceres = set(essential) & set(ceres_genes_linea.index) 
                    non_essential_in_ceres = set(non_essential) & set(ceres_genes_linea.index)
                    print('Numero de genes esenciales que estan en CERES:', str(len(ceres_genes_linea.loc[essential_in_ceres][ceres_genes_linea>-0.2])))
                    print('Numero de genes no esenciales que estan en CERES:', str(len(ceres_genes_linea.loc[essential_in_ceres])))
                    
                    
                    #Añades al data frame de la frecuencia de los genes esenciales, si el gen esta como esencial o no esta.
                    print('Añadiendo al df de la frecuencia de esenciales, presencia o ausencia de gen esencial en el modelo')
                    for gene in ceres_genes_linea.index:
                        if gene in essential_in_ceres:
                            df_frecuencia_esenciales.at[gene,th_l+'_'+th_u] = 1
                        else:
                            df_frecuencia_esenciales.at[gene,th_l+'_'+th_u] = 0
                    
                                                    
                    print('Mann Whitney')
                    x = ceres_genes_linea.loc[essential_in_ceres]
                    y = ceres_genes_linea.loc[non_essential_in_ceres]
                    U, p = mannwhitneyu(x, y, use_continuity=True)
                    
                            
                    #Calculas media de la expresion y del Score Ceres de los genes predichos como esenciales y su suma.
                    df_expresion_t = df_expresion.T
                    expr_genes_linea = df_expresion_t[cell_line]
                            
                    ScoreCeres_genes_predict = ceres_genes_linea.loc[essential_in_ceres]
                    Expr_genes_predict = expr_genes_linea.loc[essential_in_ceres]

                    mean_ScoreCeres_predict = ScoreCeres_genes_predict.mean()
                    mean_Expr_predict = Expr_genes_predict.mean()
                    valor_ceres_suma = ScoreCeres_genes_predict.sum()
                    
                    total_rx = len(csm2.reactions)
                    Number_essential = len(essential)
                    
                    print('Guardado el numero de reacciones totales, numero de genes predichos como esenciales y su suma en ceres')
                    line = th_l + '_' + th_u+'_LocalMEDIA' + '\t' + str(total_rx) + '\t'+ str(Number_essential) + '\t' + str(valor_ceres_suma)+ '\n'
                    Number_reactions.write(line)
                            
                    print('Guardado p-valor')
                    line = th_l + '_' + th_u+'_LocalMEDIA'+ '\t' + str("%.4e" % p) + '\t'+ str(mean_ScoreCeres_predict) + '\t'+ str(mean_Expr_predict) + '\n'
                    p_valor.write(line)
                    
                    print('Guardando grafico')
                    n_essential = str(len(ceres_genes_linea.loc[essential_in_ceres]))
                    n_non_essential = str(len(ceres_genes_linea.loc[non_essential_in_ceres]))
                    fig = plt.figure()
                    ax1 = fig.add_subplot(1,1,1) #(fila, columna, index)

                    #Dibujo con seaborn en el subplot(ax=ax1)
                    sns.despine()
                    sns.distplot(ceres_genes_linea.loc[essential_in_ceres], kde=True, color='red', ax=ax1).set_title('P05_P90_ThLOCAL_media', size= 18, loc='center', fontweight="bold")
                    sns.distplot(ceres_genes_linea.loc[non_essential_in_ceres], kde=True, color='blue', ax=ax1)

                    #Parametros del subplot
                    ax1.set_xlabel('Esenciabilidad\n'+ cell_line, size = 16)
                    ax1.set_ylabel('Densidad', size = 16)
                    ax1.set_xlim([-3.0, 1.0])
                    ax1.set_ylim([0.0, 2.5])
                    ax1.text(-2.8,2.2,'No_esenciales= '+n_non_essential, color = 'black', fontsize=16)
                    ax1.text(-2.8,2,'Esenciales= '+n_essential, color = 'black', fontsize=16)
                    ax1.text(-2.8,1.8,'p-valor= '+str("%.4e" % p), color= 'black', fontsize=16)
                                                    
                    #Guardando la figura
                    fig.set_size_inches(8, 7)
                    fig_ruta = './thresholding/con_th_locales2/' + cell_line + '/grafico_%s_%s.png' % (th_l, th_u) 
                    fig.savefig(fig_ruta)
                df_frecuencia_esenciales.to_csv('./thresholding/con_th_locales2/' + cell_line + '/df_frecuencia_esenciales.csv', index=True, header=True)

