#!/bin/bash
#PBS -A iwg-952-aa  # Identifiant Rap; ID
#PBS -l walltime=6:00:00    # Walltime
#PBS -l nodes=1:ppn=1  # Nombre de noeuds.
#PBS -l pmem=7900m

cd "${PBS_O_WORKDIR}"

# Modules
module load ifort_icc/15.0
module load openmpi/1.8.3-intel
module load python/2.6.7

# Run
# mpirun -np 1 ./bin/epoch3d <<< $dir/pulse
# python $dir/interpulse.py $n $nb_snapshots $dir
# qsub $dir/job_pulse$((n+1))

