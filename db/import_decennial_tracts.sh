#!/bin/bash

# mkdir -p tmp/dec/
cd tmp/dec/
# wget -r ftp://ftp2.census.gov/geo/tiger/PREVGENZ/tr/tr90shp/
# wget -r ftp://ftp2.census.gov/geo/tiger/PREVGENZ/tr/tr00shp/
# wget -r ftp://ftp2.census.gov/geo/tiger/GENZ2010/gz_2010_*_140_00_500k.zip
# mv ftp2.census.gov/geo/tiger/PREVGENZ/tr/tr90shp/*zip dec/
# mv ftp2.census.gov/geo/tiger/PREVGENZ/tr/tr00shp/*zip dec/
# mv ftp2.census.gov/geo/tiger/GENZ2010/*zip            dec/
# cd dec
# for x in `ls *zip`; do unzip $x; done

psql -d census -U jsaxon << EOD
  -- Shapefile type: Polygon
  -- Postgis type: MULTIPOLYGON[2]
  SET CLIENT_ENCODING TO UTF8;
  SET STANDARD_CONFORMING_STRINGS TO ON;
  BEGIN;

  DROP TABLE IF EXISTS "public"."census_tracts_1990";
  CREATE TABLE "public"."census_tracts_1990" (gid serial,
    "area" float8,
    "perimeter" float8,
    "N" float8,
    "Ni" float8,
    "st" smallint,
    "co" smallint,
    "tractbase" varchar(4),
    "tractsuf" varchar(2),
    "tract_name" varchar(7));
  SELECT AddGeometryColumn('public','census_tracts_1990','geom','2163','MULTIPOLYGON',2);

  DROP TABLE IF EXISTS "public"."census_tracts_2000";
  CREATE TABLE "public"."census_tracts_2000" (gid serial,
    "area" numeric,
    "perimeter" numeric,
    "N" float8,
    "Ni" float8,
    "state" smallint,
    "county" smallint,
    "tract" int,
    "name" float(8),
    "lsad" varchar(2),
    "lsad_trans" varchar(50));
  SELECT AddGeometryColumn('public','census_tracts_2000','geom','2163','MULTIPOLYGON',2);

  DROP TABLE IF EXISTS "public"."census_tracts_2010";
  CREATE TABLE "public"."census_tracts_2010" (gid serial,
    "geo_id" varchar(60),
    "state" smallint,
    "county" smallint,
    "tract" int,
    "name" varchar(90),
    "lsad" varchar(7),
    "censusarea" numeric);
  SELECT AddGeometryColumn('public','census_tracts_2010','geom','2163','MULTIPOLYGON',2);

  COMMIT;

EOD


for x in $(ls tr*_d90.shp | sed "s/.shp//"); do 
  echo $x
  shp2pgsql -I -s 4269:2163 -a -W "latin1" $x public.census_tracts_1990 | 
    sed "s/${x}_/N/g" | grep -v "GIST\|ANALYZE" | psql -d census -U jsaxon
done

for x in $(ls tr*_d00.shp | sed "s/.shp//"); do 
  echo $x
  shp2pgsql -I -s 4269:2163 -a -W "latin1" $x public.census_tracts_2000 |
    sed "s/${x}_/N/g" | grep -v "GIST\|ANALYZE" | psql -d census -U jsaxon
done

for x in $(ls gz*shp | sed "s/.shp//"); do 
  echo $x
  shp2pgsql -I -s 4269:2163 -a -W "latin1" $x public.census_tracts_2010 | grep -v "GIST\|ANALYZE" | psql -d census -U jsaxon
done



psql -d census -U jsaxon << EOD


  ALTER TABLE census_tracts_1990 ADD COLUMN tract int;
  UPDATE census_tracts_1990 SET tract = CONCAT(tractbase, CASE WHEN tractsuf != '' THEN tractsuf ELSE '00' END)::int;
  ALTER TABLE census_tracts_1990 DROP COLUMN "N", DROP COLUMN "Ni",
                                 DROP COLUMN tractbase, DROP COLUMN tractsuf, DROP COLUMN tract_name;
  ALTER TABLE census_tracts_1990 RENAME COLUMN area TO censusarea;

  ALTER TABLE census_tracts_1990 RENAME COLUMN st TO state;
  ALTER TABLE census_tracts_1990 RENAME COLUMN co TO county;

  UPDATE census_tracts_2000 SET tract = (name * 100)::int;
  ALTER TABLE census_tracts_2000 DROP COLUMN "N", DROP COLUMN "Ni", 
                                 DROP COLUMN name, DROP COLUMN lsad, DROP COLUMN lsad_trans;
  ALTER TABLE census_tracts_2000 RENAME COLUMN area TO censusarea;

  ALTER TABLE census_tracts_2010 DROP COLUMN name, DROP COLUMN lsad;
  ALTER TABLE census_tracts_2010 RENAME COLUMN geo_id TO geoid;


	-- Ellis Island duplicated
	DELETE FROM census_tracts_1990 WHERE state IS NULL;	

	UPDATE census_tracts_1990 SET 
	  geom = g, area = a, perimeter = p
	FROM (
	  SELECT 
	    ST_Multi(ST_Union(geom)) g, SUM(area) a, SUM(perimeter) p, 
	    count(gid) n, state s, county c, tract t 
	  FROM census_tracts_1990
	  GROUP BY s, c, t HAVING count(gid) > 1) AS gr
	WHERE state = s AND county = c AND tract = t;
	
	DELETE FROM census_tracts_1990
	WHERE gid IN 
	 (SELECT gid FROM 
	   (SELECT gid, ROW_NUMBER() OVER (partition BY state, county, tract ORDER BY gid) AS rnum
	    FROM census_tracts_1990) t
	    WHERE t.rnum > 1);
	
	UPDATE census_tracts_2000 SET 
	  geom = g, area = a, perimeter = p
	FROM (
	  SELECT 
	    ST_Multi(ST_Union(geom)) g, SUM(area) a, SUM(perimeter) p, 
	    count(gid) n, state s, county c, tract t 
	  FROM census_tracts_2000 
	  GROUP BY s, c, t HAVING count(gid) > 1) AS gr
	WHERE state = s AND county = c AND tract = t;
	
	DELETE FROM census_tracts_2000
	WHERE gid IN 
	 (SELECT gid FROM 
	   (SELECT gid, ROW_NUMBER() OVER (partition BY state, county, tract ORDER BY gid) AS rnum
	    FROM census_tracts_2000) t
	    WHERE t.rnum > 1);

	ALTER TABLE census_tracts_1990 DROP COLUMN gid;
	ALTER TABLE census_tracts_2000 DROP COLUMN gid;
	ALTER TABLE census_tracts_2010 DROP COLUMN gid;

  ALTER TABLE census_tracts_1990 ADD COLUMN area float;
  UPDATE census_tracts_1990 SET area = ST_Area(ST_Transform(census_tracts_1990.geom, epsg)) FROM states WHERE state = fips;

  ALTER TABLE census_tracts_2000 ADD COLUMN area float;
  UPDATE census_tracts_2000 SET area = ST_Area(ST_Transform(census_tracts_2000.geom, epsg)) FROM states WHERE state = fips;

  ALTER TABLE census_tracts_2010 ADD COLUMN area float;
  UPDATE census_tracts_2010 SET area = ST_Area(ST_Transform(census_tracts_2010.geom, epsg)) FROM states WHERE state = fips;

  ALTER TABLE census_tracts_1990 ADD COLUMN pop INT DEFAULT 0;
  ALTER TABLE census_tracts_1990 ADD COLUMN black INT DEFAULT 0;
  ALTER TABLE census_tracts_1990 ADD COLUMN hispanic INT DEFAULT 0;
  ALTER TABLE census_tracts_1990 ADD COLUMN vap  INT DEFAULT 0;
  ALTER TABLE census_tracts_1990 ADD COLUMN bvap INT DEFAULT 0;
  ALTER TABLE census_tracts_1990 ADD COLUMN hvap INT DEFAULT 0;
  ALTER TABLE census_tracts_1990 ADD centroid GEOMETRY;
  UPDATE census_tracts_1990 SET centroid = ST_Centroid(geom);
  ALTER TABLE "public"."census_tracts_1990" ADD PRIMARY KEY (state, county, tract);
  CREATE INDEX ON "public"."census_tracts_1990" USING GIST ("geom");
  ANALYZE "public"."census_tracts_1990";

  ALTER TABLE census_tracts_2000 ADD COLUMN pop INT DEFAULT 0;
  ALTER TABLE census_tracts_2000 ADD COLUMN black INT DEFAULT 0;
  ALTER TABLE census_tracts_2000 ADD COLUMN hispanic INT DEFAULT 0;
  ALTER TABLE census_tracts_2000 ADD COLUMN vap  INT DEFAULT 0;
  ALTER TABLE census_tracts_2000 ADD COLUMN bvap INT DEFAULT 0;
  ALTER TABLE census_tracts_2000 ADD COLUMN hvap INT DEFAULT 0;
  ALTER TABLE census_tracts_2000 ADD centroid GEOMETRY;
  UPDATE census_tracts_2000 SET centroid = ST_Centroid(geom);
  ALTER TABLE "public"."census_tracts_2000" ADD PRIMARY KEY (state, county, tract);
  CREATE INDEX ON "public"."census_tracts_2000" USING GIST ("geom");
  ANALYZE "public"."census_tracts_2000";

  ALTER TABLE census_tracts_2010 ADD COLUMN pop INT DEFAULT 0;
  ALTER TABLE census_tracts_2010 ADD COLUMN black INT DEFAULT 0;
  ALTER TABLE census_tracts_2010 ADD COLUMN hispanic INT DEFAULT 0;
  ALTER TABLE census_tracts_2010 ADD COLUMN vap  INT DEFAULT 0;
  ALTER TABLE census_tracts_2010 ADD COLUMN bvap INT DEFAULT 0;
  ALTER TABLE census_tracts_2010 ADD COLUMN hvap INT DEFAULT 0;
  ALTER TABLE census_tracts_2010 ADD centroid GEOMETRY;
  UPDATE census_tracts_2010 SET centroid = ST_Centroid(geom);
  ALTER TABLE "public"."census_tracts_2010" ADD PRIMARY KEY (state, county, tract);
  CREATE INDEX ON "public"."census_tracts_2010" USING GIST ("geom");
  ANALYZE "public"."census_tracts_2010";


EOD


cd ../../
./get_dec_populations.py

## rm -rf tmp


