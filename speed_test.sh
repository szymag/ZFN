#!/bin/bash

NS=(11 23 27 29 35 37 41 43 47 51 57 65 87 103)
MAIN_DIR="speed_test_results"

mkdir -v -p "${MAIN_DIR}"

for i in `seq 0 $((${#NS[*]} - 1))` ; do
    vec_num=${NS[${i}]}
    vec_num_sq=$((${vec_num} * ${vec_num}))
    test_dir="${MAIN_DIR}/${vec_num_sq}"
    mkdir -v -p "${test_dir}"
    echo -e -n "Running test for ${vec_num_sq} ... "
    time_file="${test_dir}/${vec_num_sq}_time.log"
    python ./src/eig_problem/ZagadnienieWlasne.py ${vec_num_sq} &> "${time_file}"
    mv ${vec_num_sq}.txt ${test_dir}/
    echo -e -n "done!"
    res=`grep "function calls" ${time_file}`
    echo "${res}"
done
