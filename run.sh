#!/bin/bash

for x in $(seq -w 1); do
  for s in ${@:1}; do

    if [[ $x -eq 0 ]]; then ./test.py -s ${s} -i split -l0 --print_init -w ${s}/split/s000 --shading target district; fi

    ./test.py -s ${s} -i power:100000 -t 0.01 -x${x} -l0     --print_init -w ${s}/power/s00${x} --shading target district 

    ./test.py -s ${s} -m dist_a       -t 0.01 -x${x} -n10000 -c 5 --conv_iter 1000 --destrand_min 5 --destrand_max 50 --tabu_length 2   --shading target district
    ./test.py -s ${s} -m dist_p       -t 0.01 -x${x} -n10000 -c 5 --conv_iter 1000 --destrand_min 5 --destrand_max 50 --tabu_length 2   --shading target district
    ./test.py -s ${s} -m dyn_radius   -t 0.01 -x${x} -n10000 -c 5 --conv_iter 500  --destrand_min 5 --destrand_max 50                   --shading target district
    ./test.py -s ${s} -m harm_radius  -t 0.01 -x${x} -n10000 -c 5 --conv_iter 500  --destrand_min 5 --destrand_max 50 --allow_trades    --shading target district
    ./test.py -s ${s} -m mean_radius  -t 0.01 -x${x} -n10000 -c 5 --conv_iter 500  --destrand_min 5 --destrand_max 50                   --shading target district
    ./test.py -s ${s} -m polsby       -t 0.01 -x${x} -n10000 -c 5 --conv_iter 1000 --destrand_min 5 --destrand_max 50 --tabu_length 2   --shading target district   
    ./test.py -s ${s} -m rohrbach     -t 0.02 -x${x} -n10000 -c 5 --conv_iter 500  --destrand_min 5 --destrand_max 50                   --shading target district     
    ./test.py -s ${s} -m reock        -t 0.01 -x${x} -n10000 -c 5 --conv_iter 500  --destrand_min 5 --destrand_max 50 --allow_trades    --shading target district  -o scc -p scc
    ./test.py -s ${s} -m ehrenburg    -t 0.01 -x${x} -n10000 -c 5 --conv_iter 500  --destrand_min 3 --destrand_max 50 --allow_trades    --shading target district  -o lic -p lic -r
    ./test.py -s ${s} -m hull_a       -t 0.01 -x${x} -n10000 -c 5 --conv_iter 500  --destrand_min 5 --destrand_max 50 --allow_trades    --shading target district 
    ./test.py -s ${s} -m hull_p       -t 0.01 -x${x} -n10000 -c 5 --conv_iter 500  --destrand_min 3 --destrand_max 50 --allow_trades    --shading target district 
    ./test.py -s ${s} -m exchange     -t 0.01 -x${x} -n10000 -c 5 --conv_iter 500  --destrand_min 3 --destrand_max 50                   --shading target district 
    ./test.py -s ${s} -m inertia_a    -t 0.01 -x${x} -n10000 -c 5 --conv_iter 500  --destrand_min 5 --destrand_max 50                   --shading target district
    ./test.py -s ${s} -m inertia_p    -t 0.01 -x${x} -n10000 -c 5 --conv_iter 500  --destrand_min 5 --destrand_max 50                   --shading target district
    ./test.py -s ${s} -m axis_ratio   -t 0.01 -x${x} -n10000 -c 5 --conv_iter 500  --destrand_min 5 --destrand_max 50                   --shading target district

    # if [[ $x -le 1 ]]; then ./test.py -s ${s} -m path_frac -t 0.05 -n2500 -x${x} --destrand_min 20 --destrand_max 5000 --shading target district; fi

  done
done

