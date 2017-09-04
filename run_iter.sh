#!/bin/bash

s=${STATE}
x=${SEED}
m=${METHOD}

echo RUNNING :: $x $s $m

C=20
N=10000
shading="none"
#shading="target district"

# POWER DIST RADII IPQ CIRCLES HULL INERTIA AXIS SPLIT PATH_FRAC

if [[ "$METHOD" == "SPLIT" ]]; then 
  ./run.py -s ${s} -i split -l0 --print_init -w ${s}/split/s001 --shading $shading
elif [[ "$METHOD" == "PATH_FRAC" ]]; then
  ./run.py -s ${s} -m path_frac    -t 0.01 -n$N -c$C --conv_iter 100 -x${x} --destrand_min 20 --destrand_max 50 --shading $shading -v 1 --tabu_length 3
elif [[ "$METHOD" == "DIST" ]]; then
  ./run.py -s ${s} -m dist_a       -t 0.01 -x${x} -n$N -c$C --conv_iter 1000 --destrand_min 5 --destrand_max 50 --tabu_length 2                  --shading $shading 
  ./run.py -s ${s} -m dist_p       -t 0.01 -x${x} -n$N -c$C --conv_iter 1000 --destrand_min 5 --destrand_max 50 --tabu_length 2                  --shading $shading 
elif [[ "$METHOD" == "RADII" ]]; then
  ./run.py -s ${s} -m dyn_radius   -t 0.01 -x${x} -n$N -c$C --conv_iter 500  --destrand_min 5 --destrand_max 50                                  --shading $shading 
  ./run.py -s ${s} -m harm_radius  -t 0.01 -x${x} -n$N -c$C --conv_iter 500  --destrand_min 5 --destrand_max 50 --allow_trades                   --shading $shading 
  ./run.py -s ${s} -m mean_radius  -t 0.01 -x${x} -n$N -c$C --conv_iter 500  --destrand_min 5 --destrand_max 50                                  --shading $shading 
  ./run.py -s ${s} -m rohrbach     -t 0.02 -x${x} -n$N -c$C --conv_iter 500  --destrand_min 5 --destrand_max 50                                  --shading $shading      
elif [[ "$METHOD" == "IPQ" ]]; then
  ./run.py -s ${s} -m polsby       -t 0.01 -x${x} -n$N -c$C --conv_iter 1000 --destrand_min 5 --destrand_max 50 --tabu_length 10 --allow_trades  --shading $shading    
elif [[ "$METHOD" == "CIRCLES" ]]; then
  ./run.py -s ${s} -m exchange     -t 0.01 -x${x} -n$N -c$C --conv_iter 500  --destrand_min 3 --destrand_max 50                                  --shading $shading  
  ./run.py -s ${s} -m reock        -t 0.01 -x${x} -n$N -c$C --conv_iter 500  --destrand_min 5 --destrand_max 50 --allow_trades                   --shading $shading  -o scc -p scc
  ./run.py -s ${s} -m ehrenburg    -t 0.01 -x${x} -n$N -c$C --conv_iter 500  --destrand_min 3 --destrand_max 50 --allow_trades                   --shading $shading  -o lic -p lic -r
elif [[ "$METHOD" == "HULL" ]]; then
  ./run.py -s ${s} -m hull_a       -t 0.01 -x${x} -n$N -c$C --conv_iter 500  --destrand_min 5 --destrand_max 50 --allow_trades                   --shading $shading  
  ./run.py -s ${s} -m hull_p       -t 0.01 -x${x} -n$N -c$C --conv_iter 500  --destrand_min 3 --destrand_max 50 --allow_trades                   --shading $shading  
elif [[ "$METHOD" == "INERTIA" ]]; then
  ./run.py -s ${s} -m inertia_a    -t 0.01 -x${x} -n$N -c$C --conv_iter 500  --destrand_min 5 --destrand_max 50                                  --shading $shading 
  ./run.py -s ${s} -m inertia_p    -t 0.01 -x${x} -n$N -c$C --conv_iter 500  --destrand_min 5 --destrand_max 50                                  --shading $shading 
elif [[ "$METHOD" == "AXIS" ]]; then
  ./run.py -s ${s} -m axis_ratio   -t 0.01 -x${x} -n$N -c$C --conv_iter 500  --destrand_min 5 --destrand_max 50                                  --shading $shading 
else
  ./run.py -s ${s} -i power:100000 -t 0.01 -x${x} -l0 -x0 --power_restart --print_init -c$C -w ${s}/power/$(printf "s%03d" $x) -m power --shading $shading
fi

