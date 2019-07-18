#!/bin/bash

mkdir tmp
cd tmp
# wget https://www2.census.gov/geo/tiger/GENZ2016/shp/cb_2016_us_zcta510_500k.zip
# unzip cb_2016_us_zcta510_500k.zip

### shp2pgsql -I -s 4269:2163 -p -W "latin1" cb_2016_us_zcta510_500k public.zcta
 
psql -d census -U jsaxon << EOD
  DROP TABLE IF EXISTS zcta;

  -- Shapefile type: Polygon
  -- Postgis type: MULTIPOLYGON[2]
  SET CLIENT_ENCODING TO UTF8;
  SET STANDARD_CONFORMING_STRINGS TO ON;
  BEGIN;
  CREATE TABLE "public"."zcta" (gid serial,
      "zcta5ce10" INTEGER,
      "affgeoid10" varchar(14),
      "geoid10" varchar(5),
      "aland10" float8,
      "awater10" float8);
  SELECT AddGeometryColumn('public','zcta','geom','2163','MULTIPOLYGON',2);
  COMMIT;
EOD

shp2pgsql -I -s 4269:2163 -a -W "latin1" cb_2016_us_zcta510_500k public.zcta | psql -d census -U jsaxon

psql -d census -U jsaxon << EOD

  ALTER TABLE zcta DROP COLUMN affgeoid10,
                   DROP COLUMN geoid10,
                   DROP COLUMN gid;

  ALTER TABLE zcta RENAME COLUMN zcta5ce10 TO zip;
  ALTER TABLE zcta RENAME COLUMN aland10 TO aland;
  ALTER TABLE zcta RENAME COLUMN awater10 TO awater;

  SELECT AddGeometryColumn('public','zcta','centroid','2163','POINT',2);
  UPDATE zcta SET centroid = ST_Centroid(geom);

  ALTER TABLE zcta ADD COLUMN state INTEGER;
  UPDATE zcta SET state = fips FROM states WHERE ST_Within(zcta.centroid, states.geom);

  ALTER TABLE "public"."zcta" ADD PRIMARY KEY (zip);
  CREATE INDEX ON "public"."zcta" USING GIST ("geom");

  ANALYZE "public"."zcta";
  
EOD


cd ../

## rm -rf tmp


