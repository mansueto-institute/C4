-- Create a temporary table for a single state.
CREATE TABLE temp_table AS SELECT * FROM census_tracts_2015 where state = XXSTATEXX;

-- Create a topology
SELECT topology.CreateTopology('temp_topo', find_srid('public', 'census_tracts_2015', 'geom'));

-- Add a layer to it.
SELECT topology.AddTopoGeometryColumn('temp_topo', 'public', 'temp_table', 'topogeom', 'MULTIPOLYGON');

-- Populate the layer and the topology
UPDATE temp_table SET topogeom = toTopoGeom(geom, 'temp_topo', 1); 

-- Simplify all edges up to 1 km
SELECT SimplifyEdgeGeom('temp_topo', edge_id, 10000) FROM temp_topo.edge;

-- Convert the TopoGeometries to Geometries for visualization
UPDATE census_tracts_2015 
SET geomsimp = temp_table.topogeom::geometry
FROM temp_table
WHERE 
  temp_table.state  = census_tracts_2015.state AND
  temp_table.county = census_tracts_2015.county AND
  temp_table.tract  = census_tracts_2015.tract
;

-- Drop the topology
SELECT topology.DropTopology('temp_topo');

DROP TABLE temp_table;


