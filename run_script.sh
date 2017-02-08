#!/usr/bin/env bash

export PYTHONPATH=`pwd`

for i in {1..4}
    do
        python ./src/eig_problem/ZagadnienieWlasne.py 1 c_coef_500.txt test_${i}.dat ${i}/1000
    done
#python ./src/eig_problem/ZagadnienieWlasne.py 40 p_coef_500*2.txt test.dat 0.2