#!/bin/bash

for i in cage_[0-9][0-9][0-9].pdb; do
    cat $i | sed "s/C   CO2/MC  CO2/" | sed "s/O1  CO2/MO1 CO2/" | sed "s/O2  CO2/MO2 CO2/" > tmp
    cp tmp $i
done
