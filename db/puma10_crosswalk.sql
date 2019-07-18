ALTER TABLE census_tracts_2000 ADD COLUMN puma10 INTEGER DEFAULT NULL;

UPDATE census_tracts_2000 tr SET puma10 = NULL;

UPDATE census_tracts_2000 tr SET puma10 = puma FROM puma pu
WHERE ST_Within(tr.centroid, pu.geom) AND
      ST_Within(ST_PointOnSurface(tr.geom), pu.geom) AND
      tr.state = pu.state;

UPDATE census_tracts_2000 tr SET puma10 = puma FROM puma pu
WHERE ST_Within(ST_PointOnSurface(tr.geom), pu.geom) AND
      tr.state = pu.state AND tr.pop > 0 AND tr.puma10 IS NULL;

UPDATE census_tracts_2000 tr SET puma10 = puma FROM puma pu
WHERE ST_Intersects(tr.geom, pu.geom) AND
      tr.state = pu.state AND tr.pop > 0 AND tr.puma10 IS NULL;

