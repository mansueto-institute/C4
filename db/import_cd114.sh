#!/bin/bash

psql -d census -U jsaxon <<< "DROP TABLE cd114;"

psql -d census -U jsaxon << EOD
  SET CLIENT_ENCODING TO UTF8;
  SET STANDARD_CONFORMING_STRINGS TO ON;
  BEGIN;
  CREATE TABLE "public"."cd114" (gid serial,
    "statefp" smallint,
    "cd114fp" smallint,
    "affgeoid" varchar(13),
    "geoid" smallint,
    "lsad" varchar(2),
    "cdsessn" smallint,
    "aland" float8,
    "awater" float8);
  SELECT AddGeometryColumn('public','cd114','geom','2163','MULTIPOLYGON',2);
  CREATE INDEX ON "public"."cd114" USING GIST ("geom");
  COMMIT;
  ANALYZE "public"."cd114";
EOD

# mkdir tmp
cd tmp
# wget http://www2.census.gov/geo/tiger/GENZ2015/shp/cb_2015_us_cd114_500k.zip
# unzip cb_2015_us_cd114_500k.zip
shp2pgsql -I -s 4326:2163 -a -W "latin1" cb_2015_us_cd114_500k public.cd114 | psql -d census -U jsaxon
cd ../
# rm -rf tmp

psql -d census -U jsaxon << EOD

  ALTER TABLE cd114 DROP COLUMN geoid,
                         DROP COLUMN gid,            
                         DROP COLUMN affgeoid;

  ALTER TABLE cd114 RENAME COLUMN statefp to fips;
  ALTER TABLE cd114 RENAME COLUMN cd114fp to cd;

  ALTER TABLE cd114 ADD PRIMARY KEY (fips, cd);

  ALTER TABLE cd114 ADD centroid GEOMETRY;
  UPDATE cd114 SET centroid = ST_Centroid(geom);

EOD



