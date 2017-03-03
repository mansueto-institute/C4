#!/bin/bash

# mkdir tmp
cd tmp
wget -r ftp://ftp2.census.gov/geo/tiger/GENZ2015/shp/cb_2015_*_tract_500k.zip 
mv ftp2.census.gov/geo/tiger/GENZ2015/shp/cb_2015_* .
for x in `ls *zip`; do unzip $x; done
 
 
psql -d census -U jsaxon << EOD
  -- Shapefile type: Polygon
  -- Postgis type: MULTIPOLYGON[2]
  SET CLIENT_ENCODING TO UTF8;
  SET STANDARD_CONFORMING_STRINGS TO ON;
  BEGIN;
  CREATE TABLE "public"."census_tracts_2015" (gid serial,
      "statefp" smallint,
      "countyfp" smallint,
      "tractce" int,
      "affgeoid" varchar(20),
      "geoid" varchar(11),
      "name" varchar(100),
      "lsad" varchar(2),
      "aland" float8,
      "awater" float8);
  SELECT AddGeometryColumn('public','census_tracts_2015','geom','2163','MULTIPOLYGON',2);
  CREATE INDEX ON "public"."census_tracts_2015" USING GIST ("geom");
  COMMIT;
  ANALYZE "public"."census_tracts_2015";
EOD


### shp2pgsql -I -s 4269:2163 -p -W "latin1" cb_2015_23_tract_500k public.census_tracts_2015

for x in $(ls cb_2015_*shp | sed "s/.shp//"); do 
  echo $x
  shp2pgsql -I -s 4269:2163 -a -W "latin1" $x public.census_tracts_2015 | psql -d census -U jsaxon
done

psql -d census -U jsaxon << EOD

  ALTER TABLE census_tracts_2015 DROP COLUMN gid,
                                 DROP COLUMN geoid,
                                 DROP COLUMN affgeoid,
                                 DROP COLUMN name,
                                 DROP COLUMN lsad;

  ALTER TABLE census_tracts_2015 RENAME COLUMN statefp  to state;
  ALTER TABLE census_tracts_2015 RENAME COLUMN countyfp to county;
  ALTER TABLE census_tracts_2015 RENAME COLUMN tractce  to tract;
  
  ALTER TABLE census_tracts_2015 ADD PRIMARY KEY (state, county, tract);
  
  ALTER TABLE census_tracts_2015 ADD centroid GEOMETRY;
  UPDATE census_tracts_2015 SET centroid = ST_Centroid(geom);

EOD


cd ../

## rm -rf tmp


# 1 2 4 5 6 8 9 10 11 12 13 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 44 45 46 47 48 49 50 51 53 54 55 56 60 66 69 72 78
# 8 15 16 17 23 33 37 42 47 48

psql -d census -U jsaxon < simplify_fn.sql

for x in \
  1 2 4 5 6 8 9 10 11 12 13 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 44 45 46 47 48 49 50 51 53 54 55 56 60 66 69 72 78
  do
  echo STATE ${x}
  cat simplify_states.sql | sed "s/XXSTATEXX/$x/" > temp.sql
  psql -d census -U jsaxon < temp.sql
done


