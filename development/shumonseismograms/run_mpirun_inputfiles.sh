#!/bin/bash
conda activate moose
make -j 8

for i in {1..18}; do

    mpirun -np 8 ../../farms-opt -i ./inputfiles/explicitdynamic_${i}.i