#!/usr/bin/env bash

export PYTHONPATH=`pwd`

for i in {1..9}
    do
        python ./src/eig_problem/ZagadnienieWlasne.py 1 p_coef_200*2.txt "dys_$i.txt" $i/100
    done