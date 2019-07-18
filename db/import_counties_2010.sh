#!/bin/bash

# mkdir tmp
# cd tmp
# wget https://www2.census.gov/geo/tiger/TIGER2010/COUNTY/2010/tl_2010_us_county10.zip
# unzip tl_2010_us_county10.zip

### shp2pgsql -I -s 4269:2163 -p -W "latin1" tl_2010_us_county10 public.counties_2010
 
psql -d census -U jsaxon << EOD
	-- Shapefile type: Polygon
	-- Postgis type: MULTIPOLYGON[2]
	SET CLIENT_ENCODING TO UTF8;
	SET STANDARD_CONFORMING_STRINGS TO ON;
	BEGIN;
  DROP TABLE IF EXISTS "public"."counties_2010";
	CREATE TABLE "public"."counties_2010" (gid serial,
	    "statefp10" smallint,
	    "countyfp10" smallint,
	    "countyns10" varchar(8),
	    "geoid10" varchar(5),
	    "name10" varchar(100),
	    "namelsad10" varchar(100),
	    "lsad10" varchar(2),
	    "classfp10" varchar(2),
	    "mtfcc10" varchar(5),
	    "csafp10" varchar(3),
	    "cbsafp10" varchar(5),
	    "metdivfp10" varchar(5),
	    "funcstat10" varchar(1),
	    "aland10" float8,
	    "awater10" float8,
	    "intptlat10" float8,
	    "intptlon10" float8);
	ALTER TABLE "public"."counties_2010" ADD PRIMARY KEY (gid);
	SELECT AddGeometryColumn('public','counties_2010','geom','2163','MULTIPOLYGON',2);
	COMMIT;
	ANALYZE "public"."counties_2010";
EOD

shp2pgsql -I -s 4269:2163 -a -W "latin1" tl_2010_us_county10 public.counties_2010 | grep -v "GIST\|ANALYZE" | psql -d census -U jsaxon

psql -d census -U jsaxon << EOD

  ALTER TABLE counties_2010 DROP COLUMN gid,
                            DROP COLUMN countyns10,
                            DROP COLUMN namelsad10,
                            DROP COLUMN lsad10,
                            DROP COLUMN classfp10,
                            DROP COLUMN mtfcc10,
                            DROP COLUMN csafp10,
                            DROP COLUMN cbsafp10,
                            DROP COLUMN metdivfp10,
                            DROP COLUMN funcstat10,
                            DROP COLUMN intptlat10,
                            DROP COLUMN intptlon10;

  ALTER TABLE counties_2010 RENAME COLUMN statefp10  to state;
  ALTER TABLE counties_2010 RENAME COLUMN countyfp10 to county;

  ALTER TABLE counties_2010 RENAME COLUMN geoid10 to geoid;
  ALTER TABLE counties_2010 RENAME COLUMN name10 to name;
  ALTER TABLE counties_2010 RENAME COLUMN aland10 to aland;
  ALTER TABLE counties_2010 RENAME COLUMN awater10 to awater;

  ALTER TABLE counties_2010 ADD PRIMARY KEY (state, county);

	CREATE INDEX ON "public"."counties_2010" USING GIST ("geom");
  ANALYZE "public"."counties_2010";
  
EOD


cd ../

## rm -rf tmp


