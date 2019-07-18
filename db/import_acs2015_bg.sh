#!/bin/bash

# mkdir tmp
cd tmp
# wget -r ftp://ftp2.census.gov/geo/tiger/GENZ2015/shp/cb_2015_*_bg_500k.zip 
# mv ftp2.census.gov/geo/tiger/GENZ2015/shp/cb_2015_*bg_* .
# for x in `ls *zip`; do unzip $x; done


psql -d census -U jsaxon << EOD
  -- Shapefile type: Polygon
  -- Postgis type: MULTIPOLYGON[2]
  SET CLIENT_ENCODING TO UTF8;
  SET STANDARD_CONFORMING_STRINGS TO ON;
  BEGIN;
  DROP TABLE census_bg_2015;
  CREATE TABLE "public"."census_bg_2015" (gid serial,
      "statefp" smallint,
      "countyfp" smallint,
      "tractce" int,
      "blkgrpce" smallint,
      "affgeoid" varchar(21),
      "geoid" varchar(12),
      "name" varchar(100),
      "lsad" varchar(2),
      "aland" float8,
      "awater" float8);
  SELECT AddGeometryColumn('public','census_bg_2015','geom','2163','MULTIPOLYGON',2);
  COMMIT;
EOD


### shp2pgsql -I -s 4269:2163 -p -W "latin1" cb_2015_42_bg_500k public.census_bg_2015

for x in $(ls cb_2015_*bg_*shp | sed "s/.shp//"); do 
  echo $x
  shp2pgsql -I -s 4269:2163 -a -W "latin1" $x public.census_bg_2015 | grep -v "GIST\|ANALYZE" | psql -d census -U jsaxon
done

psql -d census -U jsaxon << EOD

  ALTER TABLE census_bg_2015 DROP COLUMN gid,
                                 DROP COLUMN geoid,
                                 DROP COLUMN affgeoid,
                                 DROP COLUMN name,
                                 DROP COLUMN lsad;

  ALTER TABLE census_bg_2015 RENAME COLUMN statefp  TO state;
  ALTER TABLE census_bg_2015 RENAME COLUMN countyfp TO county;
  ALTER TABLE census_bg_2015 RENAME COLUMN tractce  TO tract;
  ALTER TABLE census_bg_2015 RENAME COLUMN blkgrpce TO bgroup;
  ALTER TABLE census_bg_2015 RENAME COLUMN aland    TO area;

  ALTER TABLE census_bg_2015 ADD COLUMN pop INTEGER;
  
  ALTER TABLE census_bg_2015 ADD PRIMARY KEY (state, county, tract, bgroup);
  
  ALTER TABLE census_bg_2015 ADD centroid GEOMETRY;
  UPDATE census_bg_2015 SET centroid = ST_Centroid(geom);

  CREATE INDEX ON "public"."census_bg_2015" USING GIST ("geom");
  -- ALTER TABLE census_bg_2015 ADD geomsimp GEOMETRY;
  -- CREATE INDEX ON "public"."census_bg_2015" USING GIST ("geomsimp");
  ANALYZE "public"."census_bg_2015";

EOD


cd ../

## rm -rf tmp


