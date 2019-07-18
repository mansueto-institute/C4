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

  ALTER TABLE census_tracts_2015 RENAME COLUMN statefp  TO state;
  ALTER TABLE census_tracts_2015 RENAME COLUMN countyfp TO county;
  ALTER TABLE census_tracts_2015 RENAME COLUMN tractce  TO tract;

  ALTER TABLE census_tracts_2015 ADD COLUMN pop INT DEFAULT 0;
  ALTER TABLE census_tracts_2015 ADD COLUMN black INT DEFAULT 0;
  ALTER TABLE census_tracts_2015 ADD COLUMN hispanic INT DEFAULT 0;
  ALTER TABLE census_tracts_2015 ADD COLUMN vap  INT DEFAULT 0;
  ALTER TABLE census_tracts_2015 ADD COLUMN bvap INT DEFAULT 0;
  ALTER TABLE census_tracts_2015 ADD COLUMN hvap INT DEFAULT 0;
  
  ALTER TABLE census_tracts_2015 ADD PRIMARY KEY (state, county, tract);
  
  ALTER TABLE census_tracts_2015 ADD centroid GEOMETRY;
  UPDATE census_tracts_2015 SET centroid = ST_Centroid(geom);

  ALTER TABLE census_tracts_2015 ADD COLUMN area float;
  UPDATE census_tracts_2015 SET area = ST_Area(ST_Transform(census_tracts_2015.geom, epsg)) FROM states WHERE state = fips;

  ALTER TABLE census_tracts_2015 ADD COLUMN geoid BIGINT;
  UPDATE census_tracts_2015 SET geoid = state::bigint * 1000000000 + county * 1000000 + tract;

  ALTER TABLE census_tracts_2015 ADD geomsimp GEOMETRY;

EOD


cd ../

## rm -rf tmp


