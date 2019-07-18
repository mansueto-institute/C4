# Database Build Scripts

This directory contains the scripts used to build the postgres database 

## Import

I have used this database for many projects, and so there are scripts for importing many different Census geographies, here.
The ones most relevant to this topic are `import_census_tracts.sh`, which imports 2015 Census 500k Cartographic Boundary File geometries, 
   and `import_cd.sh`, which similarly retrieves the geometries of actual Congressional Districts.
These scripts each have three stages: 
 1. I edit the output of `shp2pgsql -I -s 4269:2163 -p -W "latin1" filename`, to generate the schema of the tables
       in the standard, EPSG 2163 Albers projection of the US.
    I convert certain fields (FIPS codes, for instance) to numeric.
    I also postpone the standard spatial indices (GIST) to run later.
 2. Run shp2pgsl with the -a flag, to append into the just-prepared tables.
 3. Cleanup.  This includes renaming and adding fields, such as centroids and areas.
    I also make "slots" for some ACS data to be added subsequently.

In addition, the `import_census_block.sh` retrieves the block point locations (to which populations are appended),
  which is used to calculate the spatial distribution of popuations for historic districts.
Those distributions are used to calculate compactness scores for historic districts, for the appendix and online maps.

## Topo Simplify

The standard ST_Simplify function simplifies polygons, not edges.  This can affect contiguity between neigbors.
The correct way to deal with this is to convert the geometry to a topology (faces, edges, nodes),
  and simplify the edges. 
The result of this is that two faces along one edge are simplified together.
I follow the approach outlined here:

https://strk.kbt.io/blog/2012/04/13/simplifying-a-map-layer-using-postgis-topology/

In short, I set the threshold of up to 10 km simplification,
  and then halve that repeatedly, until the simplification succeeds (does not reate a new node).
The function is in `simplify_fn.sql` and the script to run it is in `topo_simplify_tracts.sql`.
In short, it is expensive and unnecessary to create the topology for the United States in one swoop,
  so I proceed state by state.

## ACS Data Retrieval

Some ACS and Decennial Census data are appended to the geometries.
The variables most used in this study are in the 2015 5-year estimes, listed in `acssf5y2015.csv`.
These were loaded into the db using `add_acs_tract_vars.py`.
