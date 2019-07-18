#!/bin/bash


mkdir tmp
cd tmp
# wget http://www2.census.gov/geo/tiger/GENZ2015/shp/cb_2015_us_cd114_500k.zip
# wget http://www2.census.gov/geo/tiger/GENZ2016/shp/cb_2016_us_cd115_500k.zip
# wget http://www2.census.gov/geo/tiger/PREVGENZ/cd/cd107shp/cd99_107_shp.zip
# wget ftp://ftp2.census.gov/geo/tiger/GENZ2010/gz_2010_*_500_11_500k.zip
# unzip cd99_107_shp.zip
# unzip cb_2015_us_cd114_500k.zip
# unzip cb_2016_us_cd115_500k.zip
# 
# for x in gz_2010_*_500_11_500k.zip; do unzip $x; done

### shp2pgsql -I -s 4269:2163 -p -W "latin1" cd99_107 public.cd107
### shp2pgsql -I -s 4269:2163 -p -W "latin1" gz_2010_us_500_11_5m public.cd111
### shp2pgsql -I -s 4269:2163 -p -W "latin1" cb_2015_us_cd114_500k public.cd114
### shp2pgsql -I -s 4269:2163 -p -W "latin1" cb_2016_us_cd115_500k public.cd115

psql -d census -U jsaxon << EOD
  DROP TABLE IF EXISTS cd;
  SET CLIENT_ENCODING TO UTF8;
  SET STANDARD_CONFORMING_STRINGS TO ON;

  BEGIN;
  CREATE TABLE "public"."cd107" (gid serial,
      "area" numeric,
      "perimeter" numeric,
      "cd99_107_" float8,
      "cd99_107_i" float8,
      "state" varchar(2),
      "cd" varchar(2),
      "lsad" varchar(2),
      "name" varchar(90),
      "lsad_trans" varchar(50));
  ALTER TABLE "public"."cd107" ADD PRIMARY KEY (gid);
  SELECT AddGeometryColumn('public','cd107','geom','2163','MULTIPOLYGON',2);
  CREATE INDEX ON "public"."cd107" USING GIST ("geom");
  COMMIT;
  ANALYZE "public"."cd107";

  BEGIN;
  CREATE TABLE "public"."cd111" (gid serial,
      "geo_id" varchar(60),
      "state" varchar(2),
      "cd" varchar(2),
      "name" varchar(90),
      "lsad" varchar(7),
      "censusarea" numeric);
  ALTER TABLE "public"."cd111" ADD PRIMARY KEY (gid);
  SELECT AddGeometryColumn('public','cd111','geom','2163','MULTIPOLYGON',2);
  CREATE INDEX ON "public"."cd111" USING GIST ("geom");
  COMMIT;
  ANALYZE "public"."cd111";

  BEGIN;
  CREATE TABLE "public"."cd114" (gid serial,
      "statefp" varchar(2),
      "cd114fp" varchar(2),
      "affgeoid" varchar(13),
      "geoid" varchar(4),
      "lsad" varchar(2),
      "cdsessn" varchar(3),
      "aland" float8,
      "awater" float8);
  ALTER TABLE "public"."cd114" ADD PRIMARY KEY (gid);
  SELECT AddGeometryColumn('public','cd114','geom','2163','MULTIPOLYGON',2);
  CREATE INDEX ON "public"."cd114" USING GIST ("geom");
  COMMIT;
  ANALYZE "public"."cd114";


  BEGIN;
  CREATE TABLE "public"."cd115" (gid serial,
      "statefp" varchar(2),
      "cd115fp" varchar(2),
      "affgeoid" varchar(13),
      "geoid" varchar(4),
      "lsad" varchar(2),
      "cdsessn" varchar(3),
      "aland" float8,
      "awater" float8);
  ALTER TABLE "public"."cd115" ADD PRIMARY KEY (gid);
  SELECT AddGeometryColumn('public','cd115','geom','2163','MULTIPOLYGON',2);
  CREATE INDEX ON "public"."cd115" USING GIST ("geom");
  COMMIT;
  ANALYZE "public"."cd115";


  BEGIN;
  CREATE TABLE "public"."cd" ("state" smallint, "cd" smallint, "sessn" smallint);
  SELECT AddGeometryColumn('public','cd','geom','2163','MULTIPOLYGON',2);
  CREATE INDEX ON "public"."cd" USING GIST ("geom");
  ALTER TABLE cd ADD PRIMARY KEY (sessn, state, cd);
  COMMIT;
  ANALYZE "public"."cd";
EOD

shp2pgsql -I -s 4269:2163 -a -W "latin1" cd99_107 public.cd107              | psql -d census -U jsaxon  
shp2pgsql -I -s 4269:2163 -a -W "latin1" cb_2015_us_cd114_500k public.cd114 | psql -d census -U jsaxon 
shp2pgsql -I -s 4269:2163 -a -W "latin1" cb_2016_us_cd115_500k public.cd115 | psql -d census -U jsaxon 

for s in $(ls gz_2010_*_500_11_500k.shp | sed "s/.shp//"); do 
  shp2pgsql -I -s 4269:2163 -a -W "latin1" $s public.cd111  | psql -d census -U jsaxon  
done

cd ../

# rm -rf tmp

psql -d census -U jsaxon << EOD

  ALTER TABLE cd107 ADD COLUMN sessn smallint DEFAULT 107;
  ALTER TABLE cd111 ADD COLUMN sessn smallint DEFAULT 111;
  ALTER TABLE cd114 ADD COLUMN sessn smallint DEFAULT 114;
  ALTER TABLE cd115 ADD COLUMN sessn smallint DEFAULT 115;

  INSERT INTO cd(sessn, state, cd, geom) SELECT sessn, state::int,   cd::int,      ST_Union(geom) FROM cd107 GROUP BY state, cd, sessn;
  INSERT INTO cd(sessn, state, cd, geom) SELECT sessn, state::int,   cd::int,      geom FROM cd111;
  INSERT INTO cd(sessn, state, cd, geom) SELECT sessn, statefp::int, cd114fp::int, geom FROM cd114;
  INSERT INTO cd(sessn, state, cd, geom) SELECT sessn, statefp::int, cd115fp::int, geom FROM cd115;

  DROP TABLE cd107;
  DROP TABLE cd111;
  DROP TABLE cd114;
  DROP TABLE cd115;

EOD



