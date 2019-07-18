#!/bin/bash

states="01 02 04 05 06 08 09 10 11 12 13 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 44 45 46 47 48 49 50 51 53 54 55 56"

# mkdir tmp
cd tmp
# for s in $states; do wget http://www2.census.gov/geo/tiger/GENZ2015/shp/cb_2015_${s}_puma10_500k.zip; done
# for x in `ls cb_2015*puma*.zip`; do unzip -o $x; done

psql -d census -U jsaxon << EOD
  DROP TABLE IF EXISTS puma;

  -- Shapefile type: Polygon
  -- Postgis type: MULTIPOLYGON[2]
  SET CLIENT_ENCODING TO UTF8;
  SET STANDARD_CONFORMING_STRINGS TO ON;
  BEGIN;
  CREATE TABLE "public"."puma" (gid serial,
      "statefp10" INTEGER,
      "pumace10" INTEGER,
      "affgeoid10" varchar(16),
      "geoid10" varchar(7),
      "name10" varchar(100),
      "lsad10" varchar(2),
      "aland10" float8,
      "awater10" float8);
  SELECT AddGeometryColumn('public','puma','geom','2163','MULTIPOLYGON',2);
  COMMIT;
EOD


### shp2pgsql -I -s 4269:2163 -p -W "latin1" cb_2015_23_puma10_500k.zip public.puma

for x in $(ls cb_2015_*puma*shp | sed "s/.shp//"); do 
  echo $x
  shp2pgsql -I -s 4269:2163 -a -W "latin1" $x public.puma | grep -v "GIST\|ANALYZE" | psql -d census -U jsaxon
done

psql -d census -U jsaxon << EOD

  ALTER TABLE puma DROP COLUMN affgeoid10,
                   DROP COLUMN geoid10,
                   DROP COLUMN lsad10,
                   DROP COLUMN gid;

  ALTER TABLE puma RENAME COLUMN statefp10 TO state;
  ALTER TABLE puma RENAME COLUMN pumace10  TO puma;
  ALTER TABLE puma RENAME COLUMN name10    TO name;
  ALTER TABLE puma RENAME COLUMN aland10   TO aland;
  ALTER TABLE puma RENAME COLUMN awater10  TO awater;

  ALTER TABLE puma ADD COLUMN pop    INT DEFAULT 0;
  ALTER TABLE puma ADD COLUMN adults INT DEFAULT 0;
  ALTER TABLE puma ADD COLUMN hs     INT DEFAULT 0;
  ALTER TABLE puma ADD COLUMN ba     INT DEFAULT 0;
   
  ALTER TABLE puma ADD PRIMARY KEY (state, puma);
  
  ALTER TABLE puma ADD centroid GEOMETRY;
  UPDATE puma SET centroid = ST_Centroid(geom);

  CREATE INDEX ON "public"."puma" USING GIST ("geom");
  ANALYZE "public"."puma";

EOD


cd ../

## rm -rf tmp


