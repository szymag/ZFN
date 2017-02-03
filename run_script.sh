#!/usr/bin/env bash

export PYTHONPATH=`pwd`

for i in {1..400}
    do
        python ./src/eig_problem/ZagadnienieWlasne.py 1 c_coef_500.txt 0.5cdys_${i}.dat ${i}/1000
    done
