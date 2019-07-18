-- Drop the topology, in case:
SELECT topology.DropTopology('block_topo');
 
-- Create a topology
SELECT topology.CreateTopology('block_topo', find_srid('public', 'census_bg_2015', 'geom'));

-- Add a layer to it.
SELECT topology.AddTopoGeometryColumn('block_topo', 'public', 'census_bg_2015', 'topogeom', 'MULTIPOLYGON');

-- Populate the layer and the topology for each state separately
-- SINGLE::
-- UPDATE census_bg_2015 SET topogeom = toTopoGeom(geom, 'block_topo', 1) WHERE state = 24; 
-- FOR LOOP:: ref/credit: http://stackoverflow.com/questions/34818875/
DO $$DECLARE r record;
BEGIN
 FOR r IN SELECT DISTINCT state FROM census_bg_2015 LOOP
  BEGIN
		RAISE NOTICE 'state = %', r.state;
    UPDATE census_bg_2015 
		SET topogeom = toTopoGeom(geom,'block_topo', 1) 
    WHERE state = r.state;
  EXCEPTION
    WHEN OTHERS THEN
     RAISE WARNING 'Loading of % failed: %', r.state, SQLERRM;
  END;
 END LOOP;
END$$;

-- Simplify all edges up to 10 km
SELECT SimplifyEdgeGeom('block_topo', edge_id, 10000) FROM block_topo.edge;
SELECT AddGeometryColumn('public','census_bg_2015','geomsimp','2163','MULTIPOLYGON',2);
UPDATE census_bg_2015 SET geomsimp = topogeom::geometry;
CREATE INDEX ON "public"."census_bg_2015" USING GIST ("geomsimp");
ANALYZE "public"."census_bg_2015";

