#!/bin/bash

states="01 02 04 05 06 08 09 10 11 12 13 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 44 45 46 47 48 49 50 51 53 54 55 56"

mkdir tb2000 tb2010

# shp2pgsql -I -s 4269:2163 -p -n -W "latin1" tl_2010_01_tabblock00.shp public.census_blocks_2000
# shp2pgsql -I -s 4269:2163 -p -n -W "latin1" tl_2010_01_tabblock10.shp public.census_blocks_2010

psql -d census -U jsaxon << EOD
  -- Shapefile type: Polygon
  -- Postgis type: MULTIPOLYGON[2]
  SET CLIENT_ENCODING TO UTF8;
  SET STANDARD_CONFORMING_STRINGS TO ON;
  BEGIN;

  DROP TABLE IF EXISTS census_blocks_2000;
  DROP TABLE IF EXISTS census_blocks_2010;

  CREATE TABLE "public"."census_blocks_2000" (gid serial,
      "statefp00" smallint,
      "countyfp00" smallint,
      "tractce00" int,
      "blockce00" int, 
      "blkidfp00" varchar(15),
      "name00" varchar(10),
      "mtfcc00" varchar(5),
      "ur00" varchar(1),
      "uace00" varchar(5),
      "funcstat00" varchar(1),
      "aland00" float8,
      "awater00" float8,
      "intptlat00" float8,
      "intptlon00" float8);

  CREATE TABLE "public"."census_blocks_2010" (gid serial,
      "statefp10" smallint,
      "countyfp10" smallint,
      "tractce10" int,
      "blockce10" int,
      "geoid10" varchar(15),
      "name10" varchar(10),
      "mtfcc10" varchar(5),
      "ur10" varchar(1),
      "uace10" varchar(5),
      "uatyp10" varchar(1),
      "funcstat10" varchar(1),
      "aland10" float8,
      "awater10" float8,
      "intptlat10" float8,
      "intptlon10" float8);

  COMMIT;
EOD


cd tb2000
for fips in $states; do 
  # wget ftp://ftp2.census.gov/geo/tiger/TIGER2010/TABBLOCK/2000/tl_2010_${fips}_tabblock00.zip
  # unzip tl_2010_${fips}_tabblock00.zip

  shp2pgsql -I -s 4269:2163 -a -n -W "latin1" tl_2010_${fips}_tabblock00.shp public.census_blocks_2000 | psql -d census -U jsaxon
done
cd -


cd tb2010
for fips in $states; do 
  # wget ftp://ftp2.census.gov/geo/tiger/TIGER2010/TABBLOCK/2010/tl_2010_${fips}_tabblock10.zip
  # unzip tl_2010_${fips}_tabblock10.zip

  shp2pgsql -I -s 4269:2163 -a -n -W "latin1" tl_2010_${fips}_tabblock10.shp public.census_blocks_2010| psql -d census -U jsaxon
done
cd -


psql -d census -U jsaxon << EOD

  ALTER TABLE census_blocks_2000 RENAME COLUMN statefp00  TO state;
  ALTER TABLE census_blocks_2000 RENAME COLUMN countyfp00 TO county;
  ALTER TABLE census_blocks_2000 RENAME COLUMN tractce00  TO tract;
  ALTER TABLE census_blocks_2000 RENAME COLUMN blockce00  TO block;
  ALTER TABLE census_blocks_2000 DROP COLUMN gid,
                                 DROP COLUMN blkidfp00,
                                 DROP COLUMN name00,
                                 DROP COLUMN mtfcc00,
                                 DROP COLUMN ur00,
                                 DROP COLUMN uace00,
                                 DROP COLUMN funcstat00;
  ALTER TABLE census_blocks_2000 RENAME COLUMN aland00  TO area;
  ALTER TABLE census_blocks_2000 RENAME COLUMN awater00 TO awater;
  ALTER TABLE census_blocks_2000 ADD PRIMARY KEY (state, county, tract, block);

  ALTER TABLE census_blocks_2010 RENAME COLUMN statefp10  TO state;
  ALTER TABLE census_blocks_2010 RENAME COLUMN countyfp10 TO county;
  ALTER TABLE census_blocks_2010 RENAME COLUMN tractce10  TO tract;
  ALTER TABLE census_blocks_2010 RENAME COLUMN blockce10  TO block;
  ALTER TABLE census_blocks_2010 DROP COLUMN gid,
                                 DROP COLUMN name10,
                                 DROP COLUMN mtfcc10,
                                 DROP COLUMN ur10,
                                 DROP COLUMN uace10,
                                 DROP COLUMN funcstat10;
  ALTER TABLE census_blocks_2010 RENAME COLUMN aland10  TO area;
  ALTER TABLE census_blocks_2010 RENAME COLUMN awater10 TO awater;
  ALTER TABLE census_blocks_2010 ADD PRIMARY KEY (state, county, tract, block);

  SELECT AddGeometryColumn('public','census_blocks_2000','geom','2163','POINT',2);
  UPDATE census_blocks_2000 SET geom = ST_Transform(ST_SetSRID(ST_Point(intptlon00, intptlat00), 4269), 2163);
  ALTER TABLE census_blocks_2000 DROP COLUMN intptlon00, DROP COLUMN intptlat00;
  CREATE INDEX ON "public"."census_blocks_2000" USING GIST ("geom");
  ANALYZE "public"."census_blocks_2000";

  SELECT AddGeometryColumn('public','census_blocks_2010','geom','2163','POINT',2);
  UPDATE census_blocks_2010 SET geom = ST_Transform(ST_SetSRID(ST_Point(intptlon10, intptlat10), 4269), 2163);
  ALTER TABLE census_blocks_2010 DROP COLUMN intptlon10, DROP COLUMN intptlat10;
  CREATE INDEX ON "public"."census_blocks_2010" USING GIST ("geom");
  ANALYZE "public"."census_blocks_2010";

  ALTER TABLE census_blocks_2000 ADD COLUMN pop INTEGER;
  ALTER TABLE census_blocks_2000 ADD COLUMN black INTEGER;
  ALTER TABLE census_blocks_2000 ADD COLUMN hispanic INTEGER;
  ALTER TABLE census_blocks_2000 ADD COLUMN vap INTEGER;
  ALTER TABLE census_blocks_2000 ADD COLUMN bvap INTEGER;
  ALTER TABLE census_blocks_2000 ADD COLUMN hvap INTEGER;

  ALTER TABLE census_blocks_2010 ADD COLUMN pop INTEGER;
  ALTER TABLE census_blocks_2010 ADD COLUMN black INTEGER;
  ALTER TABLE census_blocks_2010 ADD COLUMN hispanic INTEGER;
  ALTER TABLE census_blocks_2010 ADD COLUMN vap INTEGER;
  ALTER TABLE census_blocks_2010 ADD COLUMN bvap INTEGER;
  ALTER TABLE census_blocks_2010 ADD COLUMN hvap INTEGER;

EOD


