MÁSTER EN BIOINFORMÁTICA APLICADA A MEDICINA PERSONALIZADA Y SALUD (2018-2019)

TITULO: Reconstrucción de modelos específicos de contexto en líneas celulares de cáncer para identificar genes esenciales metabólicos y predecir nuevas dianas terapéuticas.

ANÁLISIS PRELIMINARES
Se debe seguir el suguiente orden:
   -Etiquetado
   -Filtrado_datos_Y_analisis_estadistico
   -Solo_para_genes_metabolicos
   -Analisis_genes
   -Grupos_genes

RECONSTRUCCION DE MODELOS ESPECÍFICOS DE CONTEXTO
Se debe seguir el siguiente orden. 
   -sustitucion_madre
   -lista_tareas
   -bash_greasy
   -Th_globales
   -Th_globales_y_localMEAN
   
 Estos scripts se ejecutaron en el supercomputador LaPalma. Como el proceso de reconstrucción no es paralelizable, para reconstruir los 13.417 modelos se sustituye en el script generico la palabra  "cell_line" por el nombre específico de cada línea célular, generandose un script por línea. Despues generamos un listado de tareas, para que cada tarea se ejecute en un core, que es lo que indicamos en el último script.

ÁNALISIS DE LOS MODELOS RECONSTRUIDOS ESPECÍFICOS DE CONTEXTO
   -Esenciales_sin_contexto
   -Df_pValores_Nrx
   -Graficos_pValorYrx
   -Analisis_genes
   -Analisis de los genes en los modelos
 
