#!/usr/bin/env bash

export PYTHONPATH=`pwd`

for i in {50..200}
    do
        python ./src/eig_problem/ZagadnienieWlasne.py 1 c_coef_100.txt 0_${i}.dat ${i}/2000
    done
#python ./src/eig_problem/ZagadnienieWlslasne.py 40 p_coef_500*2.txt test.dat 0.2