#!/bin/bash

for s in va; do
for x in $(seq 300 349); do
for m in POWER DIST RADII IPQ CIRCLES HULL INERTIA AXIS SPLIT PATH_FRAC; do

  aws batch submit-job \
      --job-name c4-${s}-${x}-${m}-03a --job-queue cluscious-queue \
      --job-definition arn:aws:batch:us-east-1:808035620362:job-definition/cluscious-def:6 \
      --container-overrides '{"environment" : [{"name" : "STATE", "value" : "'${s}'"}, {"name" : "SEED", "value" : "'${x}'"}, {"name" : "METHOD", "value" : "'${m}'"}]}'

done
done
done

