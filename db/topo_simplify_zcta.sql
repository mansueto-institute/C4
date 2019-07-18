-- Add a column -- 
SELECT AddGeometryColumn('public','zcta','geomsimp','2163','MULTIPOLYGON',2);

-- Create a topology
SELECT topology.DropTopology('zcta_topo');
SELECT topology.CreateTopology('zcta_topo', find_srid('public', 'zcta', 'geom'));

-- Add a layer to it.
SELECT topology.AddTopoGeometryColumn('zcta_topo', 'public', 'zcta', 'topogeom', 'MULTIPOLYGON');

-- Populate the layer and the topology for each state separately
UPDATE zcta SET topogeom = toTopoGeom(geom, 'zcta_topo', 1);

-- Simplify all edges up to 1 km
SELECT SimplifyEdgeGeom('zcta_topo', edge_id, 1000) FROM zcta_topo.edge;

UPDATE zcta SET geomsimp = topogeom::geometry;

