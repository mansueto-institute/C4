#!/bin/bash

echo RUNNING :: $STATE $METHOD $SEED $SCALING

s=${STATE}
x=${SEED}
m=${METHOD}

C="${CYCLES:=20}"
N="${ITERATIONS:=10000}"
shading="${SHADING:=none}"
scaling="${SCALING:=1.0}"

# POWER DIST RADII IPQ CIRCLES HULL_P HULL_A INERTIA AXIS SPLIT PATH_FRAC


if [[ "$METHOD" == "POWER" ]]; then 
  ./run.py -s ${s} -i power:100000 -t 0.01 -x ${x} -l0  -c100 --power_restart -w ${s}/power_regimes/s${x}/    --scale_regimes           $scaling --shading $shading --print_init -m power -v 1
fi


if [[ "$METHOD" == "HULL" ]]; then 
  ./run.py -s ${s} -m hull_a       -t 0.01 -x ${x} -n$N -c$C  --conv_iter 250 -w ${s}/hull_a_regimes/s${x}/   --scale_regimes           $scaling --shading $shading --destrand_min 5 --destrand_max 50 --allow_trades -o hull -v 1
fi

if [[ "$METHOD" == "IPQ" ]]; then 
  ./run.py -s ${s} -m polsby_w     -t 0.01 -x ${x} -n$N -c$C  --conv_iter 500 -w ${s}/polsby_w_regimes/s${x}/ --scale_regime_perimeters $scaling --shading $shading --destrand_min 5 --destrand_max 50 --tabu_length 10 --allow_trades --ctol 0.02 -v 1
fi


