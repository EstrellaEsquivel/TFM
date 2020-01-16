#!/bin/bash
# Haces un listado del nombre de los scripts (un scrip por linea celular)
ls /home/estrella/TFM/Para_LaPalma/script_per_line_Lower_locales | sed 's/^/python /g' > lista_tareas_total.txt 

