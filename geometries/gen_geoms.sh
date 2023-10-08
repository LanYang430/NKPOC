#!/bin/bash

n_max=192

module purge
module load gromacs/intel/ips2018/u1/2019.5

gmx_mpi editconf -f em.gro -o cage_000.pdb

i=187
while [ $i -lt $n_max ]; do
    label_curr=`./padding.py $i`
    label_next=`./padding.py $((i+1))`
    gmx_mpi insert-molecules -f cage_${label_curr}.pdb -ci co2.gro -nmol 1 -o cage_${label_next}.pdb -try 1000 >& log_$i
    i=$((i+1))
done

