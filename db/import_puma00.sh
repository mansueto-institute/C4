#!/bin/bash

states="01 02 04 05 06 08 09 10 11 12 13 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 44 45 46 47 48 49 50 51 53 54 55 56"

# mkdir tmp
cd tmp
# for s in $states; do wget https://www2.census.gov/geo/tiger/PREVGENZ/pu/p500shp/p5${s}_d00_shp.zip; done
# for x in p5*_d00_shp.zip; do unzip -o $x; done

psql -d census -U jsaxon << EOD
  DROP TABLE IF EXISTS puma00;

	SET CLIENT_ENCODING TO UTF8;
	SET STANDARD_CONFORMING_STRINGS TO ON;
	BEGIN;
	CREATE TABLE "public"."puma00" (gid serial,
			"area" numeric,
			"perimeter" numeric,
			"p5_d00_" float8,
			"p5_d00_i" float8,
			"puma5" integer,
			"name" varchar(90),
			"lsad" varchar(2),
			"lsad_trans" varchar(50));
	SELECT AddGeometryColumn('public','puma00','geom','2163','MULTIPOLYGON',2);
  COMMIT;
EOD


### shp2pgsql -I -s 4269:2163 -p -W "latin1" cb_2015_23_puma10_500k.zip public.puma

for x in $states; do 
  echo $x
  shp2pgsql -I -s 4269:2163 -a -W "latin1" p5${x}_d00 public.puma00 | sed "s/p5${x}_d00/p5_d00/g" | grep -v "GIST\|ANALYZE" | psql -d census -U jsaxon
done

psql -d census -U jsaxon << EOD

  ALTER TABLE puma00 DROP COLUMN p5_d00_,
                     DROP COLUMN p5_d00_i,
                     DROP COLUMN name,
                     DROP COLUMN lsad,
                     DROP COLUMN lsad_trans;

  ALTER TABLE puma00 RENAME COLUMN puma5 TO puma;

  ALTER TABLE puma00 ADD centroid GEOMETRY;
  UPDATE puma00 SET centroid = ST_Centroid(geom);

  ALTER TABLE puma00 ADD state INTEGER;

  UPDATE puma00 SET state = NULL;
  UPDATE puma00 SET state = states.fips FROM states
  WHERE ST_Within(puma00.centroid, states.geom);

  UPDATE puma00 SET state = states.fips FROM states 
  WHERE ST_Intersects(ST_Buffer(puma00.geom, -100), states.geom) AND state IS NULL;

  UPDATE puma00 SET state = 15 WHERE puma =  301 AND state IS NULL;
  UPDATE puma00 SET state = 15 WHERE puma =  302 AND state IS NULL;
  UPDATE puma00 SET state =  6 WHERE puma = 6121 AND state IS NULL;
  UPDATE puma00 SET state = 12 WHERE puma = 4020 AND state IS NULL;

  UPDATE puma00 SET 
    geom = g, area = a, perimeter = perim
  FROM (
    SELECT 
      ST_Multi(ST_Union(geom)) g, SUM(area) a, SUM(perimeter) perim, 
      count(gid) n, state s, puma p
    FROM puma00
    GROUP BY s, p HAVING count(gid) > 1) AS gr
  WHERE state = s AND puma = p;
  
  DELETE FROM puma00 WHERE gid IN 
   (SELECT gid FROM 
     (SELECT gid, ROW_NUMBER() OVER (partition BY state, puma ORDER BY gid) AS rnum
      FROM puma00) t
      WHERE t.rnum > 1);

	ALTER TABLE puma00 DROP COLUMN gid;

  ALTER TABLE puma00 ADD PRIMARY KEY (state, puma);

	CREATE INDEX ON "public"."puma00" USING GIST ("geom");
	ANALYZE "public"."puma00";

EOD


psql -d census -U jsaxon << EOD

  ALTER TABLE census_tracts_2010 ADD COLUMN puma00 INTEGER DEFAULT NULL;

  UPDATE census_tracts_2010 tr SET puma00 = NULL;

  UPDATE census_tracts_2010 tr SET puma00 = puma FROM puma00 pu
  WHERE ST_Within(tr.centroid, pu.geom) AND
        ST_Within(ST_PointOnSurface(tr.geom), pu.geom) AND
        tr.state = pu.state;

  UPDATE census_tracts_2010 tr SET puma00 = puma FROM puma00 pu
  WHERE ST_Within(ST_PointOnSurface(tr.geom), pu.geom) AND
        tr.state = pu.state AND tr.pop > 0 AND tr.puma00 IS NULL;

  UPDATE census_tracts_2010 tr SET puma00 = puma FROM puma00 pu
  WHERE ST_Intersects(tr.geom, pu.geom) AND
        tr.state = pu.state AND tr.pop > 0 AND tr.puma00 IS NULL;

EOD


cd ../

## rm -rf tmp


