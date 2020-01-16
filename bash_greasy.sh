#!/usr/bin/env bash

#SBATCH -n 160
#SBATCH --cpus-per-task 1
##SBATCH -n 48
##SBATCH --ntasks-per-core 1
#SBATCH --ntasks-per-node 16
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH -J "doble_th"
#SBATCH -t 03-00:00:00
##SBATCH -t 02:00:00
##SBATCH -queue=debug

SECONDS=0

module purge
module load python/3.6.4
module load gnu/4.8.5 openmpi/gnu/3.0.1 greasy/2.2/gnu/4.8.5/openmpi/3.0.1

SECONDS=0

#/apps/GREASY/2.2/GNU/4.8.5/OPENMPI/3.0.1/bin/greasy lista_tareas_total.txt-92309.rst
/apps/GREASY/2.2/GNU/4.8.5/OPENMPI/3.0.1/bin/greasy lista_tareas_total.txt


duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
