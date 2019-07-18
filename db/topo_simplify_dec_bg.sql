-- Drop the topology, in case:
SELECT topology.DropTopology('topo_bg_2010');
 
-- Create a topology
SELECT topology.CreateTopology('topo_bg_2010', find_srid('public', 'census_bg_2010', 'geom'));

-- Add a layer to it.
SELECT topology.AddTopoGeometryColumn('topo_bg_2010', 'public', 'census_bg_2010', 'topogeom', 'MULTIPOLYGON');

-- Populate the layer and the topology for each state separately
-- SINGLE::
-- UPDATE census_bg_2010 SET topogeom = toTopoGeom(geom, 'topo_bg_2010', 1) WHERE state = 24; 
-- FOR LOOP:: ref/credit: http://stackoverflow.com/questions/34818875/
DO $$DECLARE r record;
BEGIN
 FOR r IN SELECT DISTINCT state FROM census_bg_2010 LOOP
  BEGIN
		RAISE NOTICE 'state = %', r.state;
    UPDATE census_bg_2010 
		SET topogeom = toTopoGeom(geom,'topo_bg_2010', 1) 
    WHERE state = r.state;
  EXCEPTION
    WHEN OTHERS THEN
     RAISE WARNING 'Loading of % failed: %', r.state, SQLERRM;
  END;
 END LOOP;
END$$;

-- Simplify all edges up to 10 km
SELECT SimplifyEdgeGeom('topo_bg_2010', edge_id, 1000) FROM topo_bg_2010.edge;
SELECT AddGeometryColumn('public','census_bg_2010','geomsimp','2163','MULTIPOLYGON',2);
UPDATE census_bg_2010 SET geomsimp = topogeom::geometry;
CREATE INDEX ON "public"."census_bg_2010" USING GIST ("geomsimp");
ANALYZE "public"."census_bg_2010";

