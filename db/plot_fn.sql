DROP FUNCTION IF EXISTS PlotTract10();
CREATE FUNCTION PlotTract10() 
  RETURNS TABLE(state smallint, county smallint, tract integer, geoid BIGINT, geometry geometry(MultiPolygon,2163))
  AS $$
  SELECT
    state, county, tract, SUBSTR(geoid,10,21)::BIGINT geoid,
    CASE 
      WHEN state NOT IN (2, 15) THEN geom
      WHEN state = 2  THEN ST_Translate(ST_Scale(ST_Rotate(geom, 0.63, -3000000, 2500000), 0.38, 0.38, 0.38), -1.5e5, -3.15e6)
      WHEN state = 15 THEN ST_Intersection(ST_Translate(ST_Rotate(geom, 0.79, -5500000, -1000000), 5250000, -1.10e6), ST_Buffer(ST_SetSRID(ST_Point(0, -2e6), 2163), 2e6))
      ELSE NULL END AS geometry
  FROM census_tracts_2010
  WHERE state < 57
  ORDER BY state, county, tract;
  $$
LANGUAGE SQL;

DROP FUNCTION IF EXISTS PlotZip();
CREATE FUNCTION PlotZip() 
  RETURNS TABLE(state int, zip int, geometry geometry(MultiPolygon,2163))
  AS $$
  SELECT
    state, zip, 
    CASE 
      WHEN STATE NOT IN (2, 15) THEN geomsimp
      WHEN STATE = 2  THEN ST_Translate(ST_Scale(ST_Rotate(geomsimp, 0.63, -3000000, 2500000), 0.38, 0.38, 0.38), -1.5e5, -3.15e6)
      WHEN STATE = 15 THEN ST_Intersection(ST_Translate(ST_Rotate(geomsimp, 0.79, -5500000, -1000000), 5250000, -1.10e6), ST_Buffer(ST_SetSRID(ST_Point(0, -2e6), 2163), 2e6))
      ELSE NULL END AS geometry
  FROM zcta
  WHERE state < 57
  ORDER BY state, zcta;
  $$
LANGUAGE SQL;

DROP FUNCTION IF EXISTS PlotState();
CREATE FUNCTION PlotState() 
  RETURNS TABLE(fips smallint, usps character varying(2), epsg smallint, seats smallint, geometry geometry(MultiPolygon,2163))
  AS $$
  SELECT
    fips, usps, epsg, seats,
    CASE 
      WHEN fips NOT IN (2, 15) THEN geom 
      WHEN fips = 2  THEN ST_Translate(ST_Scale(ST_Rotate(geom, 0.63, -3000000, 2500000), 0.38, 0.38, 0.38), -1.5e5, -3.15e6)
      WHEN fips = 15 THEN ST_Intersection(ST_Translate(ST_Rotate(geom, 0.79, -5500000, -1000000), 5250000, -1.10e6), ST_Buffer(ST_SetSRID(ST_Point(0, -2e6), 2163), 2e6))
      ELSE NULL END AS geometry
  FROM states
  WHERE fips < 57
  ORDER BY fips;
  $$
LANGUAGE SQL;

DROP FUNCTION IF EXISTS PlotCounty();
CREATE FUNCTION PlotCounty() 
  RETURNS TABLE(state smallint, county smallint, pop int, geometry geometry(MultiPolygon,2163))
  AS $$
  SELECT
    state, county, pop,
    CASE 
      WHEN state NOT IN (2, 15) THEN geomland
      WHEN state = 2  THEN ST_Translate(ST_Scale(ST_Rotate(geomland, 0.63, -3000000, 2500000), 0.38, 0.38, 0.38), -1.5e5, -3.15e6)
      WHEN state = 15 THEN ST_Intersection(ST_Translate(ST_Rotate(geomland, 0.79, -5500000, -1000000), 5250000, -1.10e6), ST_Buffer(ST_SetSRID(ST_Point(0, -2e6), 2163), 2e6))
      ELSE NULL END AS geometry
  FROM counties_2015
  WHERE state < 57
  ORDER BY state, county;
  $$
LANGUAGE SQL;


DROP FUNCTION IF EXISTS PlotPuma();
CREATE FUNCTION PlotPuma() 
  RETURNS TABLE(state int, puma int, pop int, geometry geometry(MultiPolygon,2163))
  AS $$
  SELECT
    state, puma, pop, 
    CASE 
      WHEN STATE NOT IN (2, 15) THEN geom
      WHEN STATE = 2  THEN ST_Translate(ST_Scale(ST_Rotate(geom, 0.63, -3000000, 2500000), 0.38, 0.38, 0.38), -1.5e5, -3.15e6)
      WHEN STATE = 15 THEN ST_Intersection(ST_Translate(ST_Rotate(geom, 0.79, -5500000, -1000000), 5250000, -1.10e6), ST_Buffer(ST_SetSRID(ST_Point(0, -2e6), 2163), 2e6))
      ELSE NULL END AS geometry
  FROM puma
  WHERE state < 57
  ORDER BY state, puma;
  $$
LANGUAGE SQL;

DROP FUNCTION IF EXISTS PlotPuma2000();
CREATE FUNCTION PlotPuma2000() 
  RETURNS TABLE(state int, puma int, geometry geometry(MultiPolygon,2163))
  AS $$
  SELECT
    state, puma, 
    CASE 
      WHEN STATE NOT IN (2, 15) THEN geom
      WHEN STATE = 2  THEN ST_Translate(ST_Scale(ST_Rotate(geom, 0.63, -3000000, 2500000), 0.38, 0.38, 0.38), -1.5e5, -3.15e6)
      WHEN STATE = 15 THEN ST_Intersection(ST_Translate(ST_Rotate(geom, 0.79, -5500000, -1000000), 5250000, -1.10e6), ST_Buffer(ST_SetSRID(ST_Point(0, -2e6), 2163), 2e6))
      ELSE NULL END AS geometry
  FROM puma00
  WHERE state < 57
  ORDER BY state, puma;
  $$
LANGUAGE SQL;

DROP FUNCTION IF EXISTS PlotDistricts114();
CREATE FUNCTION PlotDistricts114() 
  RETURNS TABLE(state smallint, cd smallint, geometry geometry(MultiPolygon,2163))
  AS $$
  SELECT
    state, cd,
    CASE 
      WHEN state NOT IN (2, 15) THEN geom
      WHEN state = 2  THEN ST_Translate(ST_Scale(ST_Rotate(geom, 0.63, -3000000, 2500000), 0.38, 0.38, 0.38), -1.5e5, -3.15e6)
      WHEN state = 15 THEN ST_Intersection(ST_Translate(ST_Rotate(geom, 0.79, -5500000, -1000000), 5250000, -1.10e6), ST_Buffer(ST_SetSRID(ST_Point(0, -2e6), 2163), 2e6))
      ELSE NULL END AS geometry
  FROM cd
  WHERE state < 57 AND sessn = 114
  ORDER BY state, cd;
  $$
LANGUAGE SQL;


