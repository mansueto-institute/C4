#!/bin/bash

states="01 02 04 05 06 08 09 10 11 12 13 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 44 45 46 47 48 49 50 51 53 54 55 56"

# mkdir tmp
cd tmp
# for s in $states; do
#   wget https://www2.census.gov/geo/tiger/TIGER2017/SLDL/tl_2017_${s}_sldl.zip
#   wget https://www2.census.gov/geo/tiger/TIGER2017/SLDU/tl_2017_${s}_sldu.zip
# done

# for x in tl_2017_*_sld*.zip; do unzip $x; done

### shp2pgsql -I -s 4269:2163 -p -W "latin1" tl_2017_01_sldl public.sldl
### shp2pgsql -I -s 4269:2163 -p -W "latin1" tl_2017_01_sldu public.sldu

psql -d census -U jsaxon << EOD

	SET CLIENT_ENCODING TO UTF8;
	SET STANDARD_CONFORMING_STRINGS TO ON;
	BEGIN;
	CREATE TABLE "public"."sldl" (gid serial,
		"statefp" varchar(2),
		"sldlst" varchar(3),
		"geoid" varchar(5),
		"namelsad" varchar(100),
		"lsad" varchar(2),
		"lsy" int,
		"mtfcc" varchar(5),
		"funcstat" varchar(1),
		"aland" float8,
		"awater" float8,
		"intptlat" varchar(11),
		"intptlon" varchar(12));
	ALTER TABLE "public"."sldl" ADD PRIMARY KEY (gid);
	SELECT AddGeometryColumn('public','sldl','geom','2163','MULTIPOLYGON',2);
	COMMIT;
	ANALYZE "public"."sldl";


	BEGIN;
	CREATE TABLE "public"."sldu" (gid serial,
		"statefp" varchar(2),
		"sldust" varchar(3),
		"geoid" varchar(5),
		"namelsad" varchar(100),
		"lsad" varchar(2),
		"lsy" varchar(4),
		"mtfcc" varchar(5),
		"funcstat" varchar(1),
		"aland" float8,
		"awater" float8,
		"intptlat" varchar(11),
		"intptlon" varchar(12));
	ALTER TABLE "public"."sldu" ADD PRIMARY KEY (gid);
	SELECT AddGeometryColumn('public','sldu','geom','2163','MULTIPOLYGON',2);
	COMMIT;
	ANALYZE "public"."sldu";

  BEGIN;
  CREATE TABLE "public"."sld" ("state" smallint, house varchar(1), id varchar(3), year int);
  SELECT AddGeometryColumn('public','sld','geom','2163','MULTIPOLYGON',2);
  CREATE INDEX ON "public"."sld" USING GIST ("geom");
  ALTER TABLE sld ADD PRIMARY KEY (state, house, id, year);
  COMMIT;
  ANALYZE "public"."sld";

EOD

for s in $(ls tl_2017_*_sldl*shp | sed "s/.shp//"); do 
  shp2pgsql -I -s 4269:2163 -a -W "latin1" $s public.sldl | grep -v "GIST\|ANALYZE" | psql -d census -U jsaxon  
done

for s in $(ls tl_2017_*_sldu*shp | sed "s/.shp//"); do 
  shp2pgsql -I -s 4269:2163 -a -W "latin1" $s public.sldu | grep -v "GIST\|ANALYZE" | psql -d census -U jsaxon  
done

cd ../

# rm -rf tmp

psql -d census -U jsaxon << EOD

  INSERT INTO sld(state, year, id, house, geom)
	SELECT statefp::int, lsy::int, sldlst, 'L', geom
	FROM sldl;

  INSERT INTO sld(state, year, id, house, geom)
	SELECT statefp::int, lsy::int, sldust, 'U', geom
	FROM sldu;

  ANALYZE "public"."sld";

  DROP TABLE sldl;
  DROP TABLE sldu;

EOD



