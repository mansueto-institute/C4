#!/bin/bash

sql_vars=""
cpy_vars=""
for v in DP02_0066 DB00001; do
  sql_vars=${sql_vars}", "$v" float"
  cpy_vars=${cpy_vars}" "$v
done

echo $sql_vars
echo $cpy_vars
./add_acs_tract_vars.py $cpy_vars

psql -d census -U jsaxon << EOD

DROP TABLE acs2015;
CREATE TABLE acs2015 (state smallint, county smallint, tract int ${sql_vars});
\\copy acs2015 FROM '/media/jsaxon/brobdingnag/data/db/census_tract_vars.csv' DELIMITER ',' CSV

EOD

