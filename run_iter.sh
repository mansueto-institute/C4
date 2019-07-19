#!/bin/bash

echo RUNNING :: $STATE $METHOD $SEED

s=${STATE}
x=${SEED}
m=${METHOD}

C="${CYCLES:=20}"
N="${ITERATIONS:=10000}"
shading="${SHADING:=none}"

# POWER DIST RADII IPQ CIRCLES HULL_P HULL_A INERTIA AXIS SPLIT PATH_FRAC

if [[ "$METHOD" == "SPLIT" ]]; then 
  ./run.py -s ${s} -i split -l0 --print_init -w ${s}/split/s001 --shading $shading
elif [[ "$METHOD" == "PATH_FRAC" ]]; then
  ./run.py -s ${s} -m path_frac    -t 0.01  -x${x} -n$N -c5  --conv_iter 100  --destrand_min 20 --destrand_max 50 --shading $shading -v 1 --tabu_length 3
elif [[ "$METHOD" == "DIST" ]]; then
  ./run.py -s ${s} -m dist_a       -t 0.01  -x${x} -n$N -c$C --conv_iter 1000 --destrand_min 5 --destrand_max 50 --tabu_length 2                  --shading $shading 
  ./run.py -s ${s} -m dist_p       -t 0.01  -x${x} -n$N -c$C --conv_iter 1000 --destrand_min 5 --destrand_max 50 --tabu_length 2                  --shading $shading 
elif [[ "$METHOD" == "RADII" ]]; then
  ./run.py -s ${s} -m rohrbach     -t 0.005 -x${x} -n$N -c$C --conv_iter 1000 --destrand_min 5 --destrand_max 50                                  --shading $shading      
  ./run.py -s ${s} -m harm_radius  -t 0.005 -x${x} -n$N -c$C --conv_iter 1000 --destrand_min 5 --destrand_max 50 --allow_trades                   --shading $shading 
  ./run.py -s ${s} -m mean_radius  -t 0.005 -x${x} -n$N -c$C --conv_iter 1000 --destrand_min 5 --destrand_max 50                                  --shading $shading 
  ./run.py -s ${s} -m dyn_radius   -t 0.005 -x${x} -n$N -c$C --conv_iter 1000 --destrand_min 5 --destrand_max 50                                  --shading $shading 
elif [[ "$METHOD" == "IPQ" ]]; then
  ./run.py -s ${s} -m polsby       -t 0.01  -x${x} -n$N -c$C --conv_iter 1000 --destrand_min 5 --destrand_max 50 --tabu_length 10 --allow_trades  --shading $shading 
# ./run.py -s ${s} -m polsby_w     -t 0.01  -x${x} -n$N -c$C --conv_iter 1000 --destrand_min 5 --destrand_max 50 --tabu_length 10 --allow_trades  --shading $shading --ctol 0.02
elif [[ "$METHOD" == "CIRCLES" ]]; then

  ./run.py -s ${s} -m exchange     -t 0.005 -x${x} -n$N -c$C --conv_iter 500  --destrand_min 5 --destrand_max 50                                  --shading $shading  

  # These both allow trades to move cells closer to the center...
  ./run.py -s ${s} -m ehrenburg    -t 0.005 -x${x} -n$N -c$C --conv_iter 500  --ctol 0.02 \
                                                                              --destrand_min 10 --destrand_max 50 --allow_trades --shading $shading -o lic -p lic -r -v 1
  ./run.py -s ${s} -m reock        -t 0.01  -x${x} -n$N -c$C --conv_iter 500  --destrand_min 5 --destrand_max 50 --allow_trades                   --shading $shading  -o scc -p scc
elif [[ "$METHOD" == "HULL_P" ]]; then
  ./run.py -s ${s} -m hull_p       -t 0.01  -x${x} -n$N -c$C --conv_iter 500  --destrand_min 3 --destrand_max 50 --allow_trades                   --shading $shading  -o hull 
elif [[ "$METHOD" == "HULL_A" ]]; then
  ./run.py -s ${s} -m hull_a       -t 0.01  -x${x} -n$N -c$C --conv_iter 500  --destrand_min 5 --destrand_max 50 --allow_trades                   --shading $shading  -o hull 
elif [[ "$METHOD" == "INERTIA" ]]; then
  ./run.py -s ${s} -m inertia_a    -t 0.01  -x${x} -n$N -c$C --conv_iter 500  --destrand_min 5 --destrand_max 50                                  --shading $shading 
  ./run.py -s ${s} -m inertia_p    -t 0.01  -x${x} -n$N -c$C --conv_iter 500  --destrand_min 5 --destrand_max 50                                  --shading $shading 
elif [[ "$METHOD" == "AXIS" ]]; then
  ./run.py -s ${s} -m axis_ratio   -t 0.01  -x${x} -n$N -c$C --conv_iter 500  --destrand_min 5 --destrand_max 50                                  --shading $shading 
else
  ./run.py -s ${s} -i power:100000 -t 0.01  -x${x} -l0 -c100 --power_restart --print_init -w ${s}/power/$(printf "s%03d" $x) -m power --shading $shading 
fi

