#!/bin/bash

##  mkdir tmp
cd tmp
##  wget -r ftp://ftp2.census.gov/geo/tiger/TIGER2010/VTD/2010/tl_2010_*_vtd10.zip
##  mv ftp2.census.gov/geo/tiger/TIGER2010/VTD/2010/* .
##  for x in `ls *zip`; do unzip $x; done

####  shp2pgsql -I -s 4269:2163 -p -W "latin1" tl_2010_42001_vtd10 public.vtd_2010

psql -d census -U jsaxon << EOD
  -- Shapefile type: Polygon
  -- Postgis type: MULTIPOLYGON[2]
  SET CLIENT_ENCODING TO UTF8;
  SET STANDARD_CONFORMING_STRINGS TO ON;
  BEGIN;
  DROP TABLE IF EXISTS "public"."vtd_2010";
  CREATE TABLE "public"."vtd_2010" (gid serial,
      "statefp10" smallint,
      "countyfp10" smallint,
      "vtdst10" varchar(6),
      "geoid10" varchar(11),
      "vtdi10" varchar(1),
      "name10" varchar(100),
      "namelsad10" varchar(100),
      "lsad10" varchar(2),
      "mtfcc10" varchar(5),
      "funcstat10" varchar(1),
      "aland10" float8,
      "awater10" float8,
      "intptlat10" float8,
      "intptlon10" float8);
  SELECT AddGeometryColumn('public','vtd_2010','geom','2163','MULTIPOLYGON',2);
  SELECT AddGeometryColumn('public','vtd_2010','centroid','2163','POINT',2);
  COMMIT;
EOD

for x in $(ls tl_2010_*shp | sed "s/.shp//"); do 
  echo $x
  shp2pgsql -I -s 4269:2163 -a -W "latin1" $x public.vtd_2010 | grep -v "GIST\|ANALYZE" | psql -d census -U jsaxon
done


psql -d census -U jsaxon << EOD

  ALTER TABLE vtd_2010 DROP COLUMN intptlat10;
  ALTER TABLE vtd_2010 DROP COLUMN intptlon10;
  ALTER TABLE vtd_2010 DROP COLUMN namelsad10;
  ALTER TABLE vtd_2010 DROP COLUMN mtfcc10;
  ALTER TABLE vtd_2010 DROP COLUMN lsad10;
  ALTER TABLE vtd_2010 DROP COLUMN funcstat10;
  ALTER TABLE vtd_2010 DROP COLUMN vtdi10;
  ALTER TABLE vtd_2010 DROP COLUMN gid;

  ALTER TABLE vtd_2010 RENAME COLUMN geoid10    to geoid;
  ALTER TABLE vtd_2010 RENAME COLUMN name10     to name;
  ALTER TABLE vtd_2010 RENAME COLUMN statefp10  to state;
  ALTER TABLE vtd_2010 RENAME COLUMN countyfp10 to county;
  ALTER TABLE vtd_2010 RENAME COLUMN vtdst10    to vtd;
  ALTER TABLE vtd_2010 RENAME COLUMN aland10    to aland;
  ALTER TABLE vtd_2010 RENAME COLUMN awater10   to awater;
  
  -- ALTER TABLE vtd_2010 ADD geomsimp GEOMETRY;

  UPDATE vtd_2010 SET centroid = ST_Centroid(geom);

  ALTER TABLE vtd_2010 ADD PRIMARY KEY (state, county, vtd);
  CREATE INDEX ON "public"."vtd_2010" USING GIST ("geom");
  ANALYZE "public"."vtd_2010";

EOD


##  cd ../

##  rm -rf tmp


