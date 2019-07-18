-- Drop the topology, in case:
SELECT topology.DropTopology('state_topo');
 
-- Create a topology
SELECT topology.CreateTopology('state_topo', find_srid('public', 'census_tracts_2015', 'geom'));

-- Add a layer to it.
SELECT topology.AddTopoGeometryColumn('state_topo', 'public', 'census_tracts_2015', 'topogeom', 'MULTIPOLYGON');

-- Populate the layer and the topology for each state separately
-- SINGLE::
-- UPDATE census_tracts_2015 SET topogeom = toTopoGeom(geom, 'state_topo', 1) WHERE state = 24; 
-- FOR LOOP:: ref/credit: http://stackoverflow.com/questions/34818875/
DO $$DECLARE r record;
BEGIN
 FOR r IN SELECT DISTINCT state FROM census_tracts_2015 LOOP
  BEGIN
		RAISE NOTICE 'state = %', r.state;
    UPDATE census_tracts_2015 
		SET topogeom = toTopoGeom(geom,'state_topo', 1) 
    WHERE state = r.state;
   EXCEPTION
    WHEN OTHERS THEN
     RAISE WARNING 'Loading of % failed: %', r.state, SQLERRM;
  END;
 END LOOP;
END$$;

-- Simplify all edges up to 1 km
SELECT SimplifyEdgeGeom('state_topo', edge_id, 10000) FROM state_topo.edge;

UPDATE census_tracts_2015 SET geomsimp = topogeom::geometry;

