#!/bin/bash

# 1 2 4 5 6 8 9 10 11 12 13 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 44 45 46 47 48 49 50 51 53 54 55 56 60 66 69 72 78

psql -d census -U jsaxon < simplify_fn.sql

for x in 42; do
  echo STATE ${x}
  cat simplify_states.sql | sed "s/XXSTATEXX/$x/" > temp.sql
  psql -d census -U jsaxon < temp.sql
done


