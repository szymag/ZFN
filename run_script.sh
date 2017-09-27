#!/usr/bin/env bash

export PYTHONPATH=`pwd`
for j in {0..90..2}
    do
    for i in {1..200}
        do
            python ./src/eig_problem/ZagadnienieWlasne.py 1 c_coef_100.txt ${j}_${i}.dat ${i}/1000 ${j}
        done
    done
#python ./src/eig_problem/ZagadnienieWlslasne.py 40 p_coef_500*2.txt test.dat 0.2