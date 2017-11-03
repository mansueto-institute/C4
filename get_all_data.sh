#!/bin/bash

states="al ak az ar ca co ct de dc fl ga hi id il in ia ks ky la me md ma mi mn ms mo mt ne nv nh nj nm ny nc nd oh ok or pa ri sc sd tn tx ut vt va wa wv wi wy"

for s in $states; do
  ./run.py -s $s -i power:100000 --power_restart --print_init -c 100 -t 0.01 -l0 -x1 -w $s/power/s001 -m power
done

