SET CLIENT_ENCODING TO UTF8;
SET STANDARD_CONFORMING_STRINGS TO ON;
BEGIN;
CREATE TABLE "public"."census_tracts_2015" (gid serial,
"statefp" integer,
"countyfp" integer,
"tractce" integer,
"geoid" bigint,
"name" varchar(7),
"namelsad" varchar(20),
"mtfcc" varchar(5),
"funcstat" varchar(1),
"aland" float8,
"awater" float8,
"intptlat" float,
"intptlon" float);
ALTER TABLE "public"."census_tracts_2015" ADD PRIMARY KEY (statefp, countyfp, tractce);
SELECT AddGeometryColumn('public','census_tracts_2015','geom','4269','MULTIPOLYGON',2);
CREATE INDEX ON "public"."census_tracts_2015" USING GIST ("geom");
COMMIT;
ANALYZE "public"."census_tracts_2015";
