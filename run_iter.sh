#!/bin/bash

s=$1
x=$2

cycles=10
max_iter=10000
# shading="target district"
shading="none"

if [[ "$METHOD" == "SPLIT" ]]; then 
  ./run.py -s ${s} -i split -l0 --print_init -w ${s}/split/s001 --shading target district
elif [[ "$METHOD" == "PATH_FRAC" ]]; then
  ./run.py -s ${s} -m path_frac -t 0.05 -n2000 --conv_iter 100 -x${x} --destrand_min 20 --destrand_max 5000 --shading target district
else
  ./run.py -s ${s} -i power:100000 -t 0.01 -x${x} -l0 -x0 --power_restart --print_init -c 100 -w ${s}/power/$(printf "s%03d" $x)/ -m power --shading $shading
  ./run.py -s ${s} -m dist_a       -t 0.01 -x${x} -n$max_iter -c $cycles --conv_iter 1000 --destrand_min 5 --destrand_max 50 --tabu_length 2   --shading $shading 
  ./run.py -s ${s} -m dist_p       -t 0.01 -x${x} -n$max_iter -c $cycles --conv_iter 1000 --destrand_min 5 --destrand_max 50 --tabu_length 2   --shading $shading 
  ./run.py -s ${s} -m dyn_radius   -t 0.01 -x${x} -n$max_iter -c $cycles --conv_iter 500  --destrand_min 5 --destrand_max 50                   --shading $shading 
  ./run.py -s ${s} -m harm_radius  -t 0.01 -x${x} -n$max_iter -c $cycles --conv_iter 500  --destrand_min 5 --destrand_max 50 --allow_trades    --shading $shading 
  ./run.py -s ${s} -m mean_radius  -t 0.01 -x${x} -n$max_iter -c $cycles --conv_iter 500  --destrand_min 5 --destrand_max 50                   --shading $shading 
  ./run.py -s ${s} -m polsby       -t 0.01 -x${x} -n$max_iter -c $cycles --conv_iter 1000 --destrand_min 5 --destrand_max 50 --tabu_length 2   --shading $shading    
  ./run.py -s ${s} -m rohrbach     -t 0.02 -x${x} -n$max_iter -c $cycles --conv_iter 500  --destrand_min 5 --destrand_max 50                   --shading $shading      
  ./run.py -s ${s} -m reock        -t 0.01 -x${x} -n$max_iter -c $cycles --conv_iter 500  --destrand_min 5 --destrand_max 50 --allow_trades    --shading $shading  -o scc -p scc
  ./run.py -s ${s} -m ehrenburg    -t 0.01 -x${x} -n$max_iter -c $cycles --conv_iter 500  --destrand_min 3 --destrand_max 50 --allow_trades    --shading $shading  -o lic -p lic -r
  ./run.py -s ${s} -m hull_a       -t 0.01 -x${x} -n$max_iter -c $cycles --conv_iter 500  --destrand_min 5 --destrand_max 50 --allow_trades    --shading $shading  
  ./run.py -s ${s} -m hull_p       -t 0.01 -x${x} -n$max_iter -c $cycles --conv_iter 500  --destrand_min 3 --destrand_max 50 --allow_trades    --shading $shading  
  ./run.py -s ${s} -m exchange     -t 0.01 -x${x} -n$max_iter -c $cycles --conv_iter 500  --destrand_min 3 --destrand_max 50                   --shading $shading  
  ./run.py -s ${s} -m inertia_a    -t 0.01 -x${x} -n$max_iter -c $cycles --conv_iter 500  --destrand_min 5 --destrand_max 50                   --shading $shading 
  ./run.py -s ${s} -m inertia_p    -t 0.01 -x${x} -n$max_iter -c $cycles --conv_iter 500  --destrand_min 5 --destrand_max 50                   --shading $shading 
  ./run.py -s ${s} -m axis_ratio   -t 0.01 -x${x} -n$max_iter -c $cycles --conv_iter 500  --destrand_min 5 --destrand_max 50                   --shading $shading 
fi

