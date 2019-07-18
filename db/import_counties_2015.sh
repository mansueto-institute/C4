#!/bin/bash

# mkdir tmp
cd tmp
wget https://www2.census.gov/geo/tiger/TIGER2015/COUNTY/tl_2015_us_county.zip
unzip tl_2015_us_county.zip

### shp2pgsql -I -s 4269:2163 -p -W "latin1" tl_2015_us_county public.counties_2015
 
psql -d census -U jsaxon << EOD
	-- Shapefile type: Polygon
	-- Postgis type: MULTIPOLYGON[2]
	SET CLIENT_ENCODING TO UTF8;
	SET STANDARD_CONFORMING_STRINGS TO ON;
	BEGIN;
  DROP TABLE IF EXISTS "public"."counties_2015";
	CREATE TABLE "public"."counties_2015" (gid serial,
	    "statefp" smallint,
	    "countyfp" smallint,
	    "countyns" varchar(8),
	    "geoid" varchar(5),
	    "name" varchar(100),
	    "namelsad" varchar(100),
	    "lsad" varchar(2),
	    "classfp" varchar(2),
	    "mtfcc" varchar(5),
	    "csafp" varchar(3),
	    "cbsafp" varchar(5),
	    "metdivfp" varchar(5),
	    "funcstat" varchar(1),
	    "aland" float8,
	    "awater" float8,
	    "intptlat" float8,
	    "intptlon" float8);
	ALTER TABLE "public"."counties_2015" ADD PRIMARY KEY (gid);
	SELECT AddGeometryColumn('public','counties_2015','geom','2163','MULTIPOLYGON',2);
	COMMIT;
	ANALYZE "public"."counties_2015";
EOD

shp2pgsql -I -s 4269:2163 -a -W "latin1" tl_2015_us_county public.counties_2015 | grep -v "GIST\|ANALYZE" | psql -d census -U jsaxon

psql -d census -U jsaxon << EOD

  ALTER TABLE counties_2015 DROP COLUMN gid,
                            DROP COLUMN countyns,
                            DROP COLUMN namelsad,
                            DROP COLUMN lsad,
                            DROP COLUMN classfp,
                            DROP COLUMN mtfcc,
                            DROP COLUMN csafp,
                            DROP COLUMN cbsafp,
                            DROP COLUMN metdivfp,
                            DROP COLUMN funcstat,
                            DROP COLUMN intptlat,
                            DROP COLUMN intptlon;

  ALTER TABLE counties_2015 RENAME COLUMN statefp  to state;
  ALTER TABLE counties_2015 RENAME COLUMN countyfp to county;

  ALTER TABLE counties_2015 ADD COLUMN pop INTEGER DEFAULT 0;
  
  ALTER TABLE counties_2015 ADD PRIMARY KEY (state, county);

	CREATE INDEX ON "public"."counties_2015" USING GIST ("geom");
  ANALYZE "public"."counties_2015";
  
EOD


cd ../

## rm -rf tmp


