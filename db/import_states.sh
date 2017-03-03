#!/bin/bash

# mkdir tmp
cd tmp
# wget ftp://ftp2.census.gov/geo/tiger/GENZ2015/shp/cb_2015_us_state_500k.zip
# unzip cb_2015_us_state_500k.zip

# shp2pgsql -I -s 4269:2163 -p -W "latin1" cb_2015_us_state_500k public.cb_2015_us_state_500k

psql -d census -U jsaxon << EOD
  SET CLIENT_ENCODING TO UTF8;
  SET STANDARD_CONFORMING_STRINGS TO ON;
  BEGIN;
  CREATE TABLE "public"."states" (gid serial,
      "statefp" smallint,
      "statens" varchar(8),
      "affgeoid" varchar(11),
      "geoid" varchar(2),
      "stusps" varchar(2),
      "name" varchar(100),
      "lsad" varchar(2),
      "aland" float8,
      "awater" float8);
  SELECT AddGeometryColumn('public','states','geom','2163','MULTIPOLYGON',2);
  CREATE INDEX ON "public"."states" USING GIST ("geom");
  COMMIT;
  ANALYZE "public"."states";
EOD

shp2pgsql -I -s 4269:2163 -a -W "latin1" cb_2015_us_state_500k public.states | psql -d census -U jsaxon

psql -d census -U jsaxon << EOD

  ALTER TABLE states DROP COLUMN gid,
                     DROP COLUMN geoid,
                     DROP COLUMN affgeoid,
                     DROP COLUMN statens,
                     DROP COLUMN lsad;

  ALTER TABLE states RENAME COLUMN statefp TO fips;
  ALTER TABLE states RENAME COLUMN stusps  TO usps;
  
  ALTER TABLE states ADD PRIMARY KEY(fips);

  ALTER TABLE states ADD centroid GEOMETRY;
  UPDATE states SET centroid = ST_Centroid(geom);

  ALTER TABLE states ADD epsg  SMALLINT;
  ALTER TABLE states ADD seats SMALLINT;

  CREATE TEMP TABLE tmp_x (fips SMALLINT, seats SMALLINT, epsg SMALLINT);
  \\copy tmp_x FROM '/media/jsaxon/brobdingnag/data/db/state_epsg.csv' DELIMITER ',' CSV

  UPDATE states
  SET    epsg  = tmp_x.epsg,
         seats = tmp_x.seats
  FROM   tmp_x
  WHERE  states.fips = tmp_x.fips;


EOD

cd ../
# rm -rf tmp


