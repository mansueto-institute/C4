#!/bin/bash

mkdir -p bg1990 bg2000 bg2010

psql -d census -U jsaxon << EOD
  -- Shapefile type: Polygon
  -- Postgis type: MULTIPOLYGON[2]
  SET CLIENT_ENCODING TO UTF8;
  SET STANDARD_CONFORMING_STRINGS TO ON;
  BEGIN;

  DROP TABLE IF EXISTS census_bg_1990;
  DROP TABLE IF EXISTS census_bg_2000;
  DROP TABLE IF EXISTS census_bg_2010;

  CREATE TABLE "public"."census_bg_1990" (gid serial,
      "area" float8,
      "perimeter" float8,
      "N" float8,
      "Ni" float8,
      "poly_" float8,
      "subclass" varchar(13),
      "subclass_" float8,
      "rings_ok" int4,
      "rings_nok" int4,
      "st" smallint,
      "co" smallint,
      "tract" int,
      "bg" smallint,
      "geoid" varchar(12),
      "arealand" float8,
      "areawat" float8,
      "areatot" float8,
      "name" varchar(66));
  
  SELECT AddGeometryColumn('public','census_bg_1990','geom','2163','MULTIPOLYGON',2);

  CREATE TABLE "public"."census_bg_2000" (gid serial,
      "area" numeric,
      "perimeter" numeric,
      "N" float8,
      "Ni" float8,
      "state" smallint,
      "county" smallint,
      "tract" int,
      "blkgroup" smallint,
      "name" varchar(90),
      "lsad" varchar(2),
      "lsad_trans" varchar(50));

  SELECT AddGeometryColumn('public','census_bg_2000','geom','2163','MULTIPOLYGON',2);

  CREATE TABLE "public"."census_bg_2010" (gid serial,
      "geo_id" varchar(60),
      "state" smallint,
      "county" smallint,
      "tract" int,
      "blkgrp" smallint,
      "name" varchar(90),
      "lsad" varchar(7),
      "censusarea" numeric);

  SELECT AddGeometryColumn('public','census_bg_2010','geom','2163','MULTIPOLYGON',2);

  COMMIT;

EOD


cd bg1990
# wget ftp://ftp.census.gov/geo/tiger/PREVGENZ/bg/bg90shp/bg*_d90_shp.zip
# for x in *zip; do unzip $x; done
for x in *shp; do
  v=$(echo $x | sed "s/.shp//")
  shp2pgsql -I -s 4269:2163 -a -W "latin1" $x public.census_bg_1990 | sed "s/${v}_/N/g" | grep -v "ANALYZE\|GIST" | psql -d census -U jsaxon
done 
cd -

cd bg2000
# for x in *zip; do unzip $x; done
# wget ftp://ftp2.census.gov/geo/tiger/PREVGENZ/bg/bg00shp/bg*_d00_shp.zip
for x in *shp; do
  v=$(echo $x | sed "s/.shp//")
  shp2pgsql -I -s 4269:2163 -a -W "latin1" $x public.census_bg_2000 | sed "s/${v}_/N/g" | grep -v "ANALYZE\|GIST" | psql -d census -U jsaxon
done 
cd -

cd bg2010
# for x in *zip; do unzip $x; done
# wget ftp://ftp2.census.gov/geo/tiger/GENZ2010/gz_2010_*_150_00_500k.zip
for x in *shp; do
  shp2pgsql -I -s 4269:2163 -a -W "latin1" $x public.census_bg_2010 | grep -v "ANALYZE\|GIST" | psql -d census -U jsaxon
done 
cd -


psql -d census -U jsaxon << EOD

  BEGIN;

  DELETE FROM census_bg_1990 WHERE st IS NULL; -- two blocks for Ellis Island.

	UPDATE census_bg_1990 SET 
	  geom = g, area = a, perimeter = p
	FROM (
	  SELECT 
	    ST_Multi(ST_Union(geom)) g, SUM(area) a, SUM(perimeter) p, 
	    count(gid) n, st s, co c, tract t, bg b
	  FROM census_bg_1990
	  GROUP BY s, c, t, b HAVING count(gid) > 1) AS gr
	WHERE st = s AND co = c AND tract = t AND bg = b;
	
	DELETE FROM census_bg_1990
	WHERE gid IN 
	 (SELECT gid FROM 
	   (SELECT gid, ROW_NUMBER() OVER (partition BY st, co, tract, bg ORDER BY gid) AS rnum
	    FROM census_bg_1990) t
	    WHERE t.rnum > 1);
	
	UPDATE census_bg_2000 SET 
	  geom = g, area = a, perimeter = p
	FROM (
	  SELECT 
	    ST_Multi(ST_Union(geom)) g, SUM(area) a, SUM(perimeter) p, 
	    count(gid) n, state s, county c, tract t, blkgroup b
	  FROM census_bg_2000 
	  GROUP BY s, c, t, b HAVING count(gid) > 1) AS gr
	WHERE state = s AND county = c AND tract = t AND blkgroup = b;
	
	DELETE FROM census_bg_2000
	WHERE gid IN 
	 (SELECT gid FROM 
	   (SELECT gid, ROW_NUMBER() OVER (partition BY state, county, tract, blkgroup ORDER BY gid) AS rnum
	    FROM census_bg_2000) t
	    WHERE t.rnum > 1);


  ALTER TABLE census_bg_1990 RENAME COLUMN st TO state;
  ALTER TABLE census_bg_1990 RENAME COLUMN co TO county;
  ALTER TABLE census_bg_1990 RENAME COLUMN bg TO bgroup;
  ALTER TABLE census_bg_1990 RENAME COLUMN arealand TO aland;
  ALTER TABLE census_bg_1990 RENAME COLUMN areawat  TO awater;
  ALTER TABLE census_bg_1990 DROP COLUMN gid,
                             DROP COLUMN perimeter,
                             DROP COLUMN "N",
                             DROP COLUMN "Ni",
                             DROP COLUMN subclass,
                             DROP COLUMN subclass_,
                             DROP COLUMN rings_ok,
                             DROP COLUMN rings_nok,
                             DROP COLUMN geoid,
                             DROP COLUMN name;
  ALTER TABLE census_bg_1990 ADD PRIMARY KEY (state, county, tract, bgroup);
  CREATE INDEX ON "public"."census_bg_1990" USING GIST ("geom");
  ANALYZE "public"."census_bg_1990";

  SELECT AddGeometryColumn('public','census_bg_1990','centroid','2163','POINT',2);
  UPDATE census_bg_1990 SET centroid = ST_Centroid(geom);
  UPDATE census_bg_1990 SET area = ST_Area(ST_Transform(census_bg_1990.geom, epsg)) FROM states WHERE state = fips;


  ALTER TABLE census_bg_2000 RENAME COLUMN blkgroup TO bgroup;
  ALTER TABLE census_bg_2000 DROP COLUMN gid,
                             DROP COLUMN perimeter,
                             DROP COLUMN "N",
                             DROP COLUMN "Ni",
                             DROP COLUMN name,
                             DROP COLUMN lsad,
                             DROP COLUMN lsad_trans;
  ALTER TABLE census_bg_2000 ADD PRIMARY KEY (state, county, tract, bgroup);
  CREATE INDEX ON "public"."census_bg_2000" USING GIST ("geom");
  ANALYZE "public"."census_bg_2000";

  SELECT AddGeometryColumn('public','census_bg_2000','centroid','2163','POINT',2);
  UPDATE census_bg_2000 SET centroid = ST_Centroid(geom);
  UPDATE census_bg_2000 SET area = ST_Area(ST_Transform(census_bg_2000.geom, epsg)) FROM states WHERE state = fips;


  ALTER TABLE census_bg_2010 RENAME COLUMN blkgrp     TO bgroup;
  ALTER TABLE census_bg_2010 DROP COLUMN gid,
                             DROP COLUMN geo_id,
                             DROP COLUMN name,
                             DROP COLUMN lsad,
                             DROP COLUMN censusarea;
  ALTER TABLE census_bg_2010 ADD PRIMARY KEY (state, county, tract, bgroup);
  CREATE INDEX ON "public"."census_bg_2010" USING GIST ("geom");
  ANALYZE "public"."census_bg_2010";

  SELECT AddGeometryColumn('public','census_bg_2010','centroid','2163','POINT',2);
  UPDATE census_bg_2010 SET centroid = ST_Centroid(geom);

  ALTER TABLE census_bg_2010 ADD COLUMN area FLOAT;
  UPDATE census_bg_2010 SET area = ST_Area(geom);
  UPDATE census_bg_2010 SET area = ST_Area(ST_Transform(census_bg_2010.geom, epsg)) FROM states WHERE state = fips;


  ALTER TABLE census_bg_1990 ADD COLUMN pop INTEGER;
  ALTER TABLE census_bg_1990 ADD COLUMN black INTEGER;
  ALTER TABLE census_bg_1990 ADD COLUMN hispanic INTEGER;
  ALTER TABLE census_bg_1990 ADD COLUMN vap INTEGER;
  ALTER TABLE census_bg_1990 ADD COLUMN bvap INTEGER;
  ALTER TABLE census_bg_1990 ADD COLUMN hvap INTEGER;

  ALTER TABLE census_bg_2000 ADD COLUMN pop INTEGER;
  ALTER TABLE census_bg_2000 ADD COLUMN black INTEGER;
  ALTER TABLE census_bg_2000 ADD COLUMN hispanic INTEGER;
  ALTER TABLE census_bg_2000 ADD COLUMN vap INTEGER;
  ALTER TABLE census_bg_2000 ADD COLUMN bvap INTEGER;
  ALTER TABLE census_bg_2000 ADD COLUMN hvap INTEGER;

  ALTER TABLE census_bg_2010 ADD COLUMN pop INTEGER;
  ALTER TABLE census_bg_2010 ADD COLUMN black INTEGER;
  ALTER TABLE census_bg_2010 ADD COLUMN hispanic INTEGER;
  ALTER TABLE census_bg_2010 ADD COLUMN vap INTEGER;
  ALTER TABLE census_bg_2010 ADD COLUMN bvap INTEGER;
  ALTER TABLE census_bg_2010 ADD COLUMN hvap INTEGER;

  COMMIT;

EOD

