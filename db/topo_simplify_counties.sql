-- Add a column -- 
SELECT AddGeometryColumn('public','counties_2010','geomsimp','2163','MULTIPOLYGON',2);

-- Create a topology
SELECT topology.DropTopology('county_topo_2010');
SELECT topology.CreateTopology('county_topo_2010', find_srid('public', 'counties_2010', 'geom'));

-- Add a layer to it.
SELECT topology.AddTopoGeometryColumn('county_topo_2010', 'public', 'counties_2010', 'topogeom', 'MULTIPOLYGON');

-- Populate the layer and the topology for each state separately
UPDATE counties_2010 SET topogeom = toTopoGeom(geom, 'county_topo_2010', 1);

-- Simplify all edges up to 1 km
SELECT SimplifyEdgeGeom('county_topo_2010', edge_id, 1000) FROM county_topo_2010.edge;

UPDATE counties_2010 SET geomsimp = topogeom::geometry;

SELECT AddGeometryColumn('public','counties_2010','geomland','2163','MULTIPOLYGON',2);

UPDATE counties_2010 
SET geomland = ST_Multi(ST_CollectionExtract(ST_Intersection(counties_2010.geomsimp, states.geom), 3))
FROM states WHERE states.fips = counties_2010.state;

SELECT AddGeometryColumn('public','counties_2010','centroid','2163','POINT',2);
UPDATE counties_2010 SET centroid = ST_Centroid(geom);


