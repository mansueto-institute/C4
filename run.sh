#!/bin/bash

for x in $(seq -w 100); do
  for s in pa tn md nj oh wi tx fl; do 

    if [[ $s -eq 1 ]];  then ./test.py -s ${s} -i split -l0 --print_init -w pa/split/s${x} --shading target district; fi
    if [[ $x -le 10 ]]; then ./test.py -s ${s} -m path_frac -t 0.05 -n500 -x${x} --destrand_min 20 --destrand_max 5000 --shading target district; fi

    ./test.py -s ${s} -i power:100000 -t 0.005 -l0    -x${x} --print_init -w pa/power/s${x}                       --shading target district 
    ./test.py -s ${s} -m dyn_radius   -t 0.01  -n2000 -x${x} --destrand_min 5 --destrand_max 20                   --shading target district
    ./test.py -s ${s} -m polsby       -t 0.005 -n2500 -x${x} --destrand_min 3 --destrand_max 50 --tabu_length 12  --shading target district   
    ./test.py -s ${s} -m rohrbach     -t 0.02  -l1500 -x${x} --destrand_min 3 --destrand_max 50                   --shading target district     
    ./test.py -s ${s} -m reock        -t 0.02  -n1000 -x${x} --destrand_min 5 --destrand_max 50 --allow_trades    --shading target district  -c scc -p scc
    ./test.py -s ${s} -m ehrenburg    -t 0.02  -n2000 -x${x} --destrand_min 5 --destrand_max 50 --allow_trades    --shading target district  -c lic -p lic
    ./test.py -s ${s} -m hull_p       -t 0.05  -n200  -x${x} --destrand_min 3 --destrand_max 30 --allow_trades    --shading target district 
    ./test.py -s ${s} -m exchange     -t 0.005 -n5000 -x${x} --destrand_min 3 --destrand_max 30                   --shading target district 

  done
done

