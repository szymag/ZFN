#!/usr/bin/env bash

export PYTHONPATH=`pwd`

python ./src/eig_problem/ZagadnienieWlasne.py 1e-9 tm_coef_100*5.txt tm5.txt 32
python ./src/eig_problem/ZagadnienieWlasne.py 1e-9 tm_coef_50*7.txt tm7.txt 128
python ./src/eig_problem/ZagadnienieWlasne.py 1e-9 tm_coef_10*9.txt tm9.txt 512
#python ./src/eig_problem/ZagadnienieWlasne.py 1e-9 tm_coef_5*11.txt tm11.txt 2048
#python ./src/eig_problem/ZagadnienieWlasne.py 1e-9 tm_coef_50*33.txt tm13.txt 8192
#python ./src/eig_problem/ZagadnienieWlasne.py 0.99 tm_coef_100*5.txt tm5.txt 32
#python ./src/eig_problem/ZagadnienieWlasne.py 0.99 tm_coef_50*7.txt tm7.txt 64
#python ./src/eig_problem/ZagadnienieWlasne.py 0.99 tm_coef_10*9.txt tm9.txt 512
#python ./src/eig_problem/ZagadnienieWlasne.py 0.99 tm_coef_5*11.txt tm11.txt 2048
#python ./src/eig_problem/ZagadnienieWlasne.py 0.99 tm_coef_50*33.txt tm13.txt 8192