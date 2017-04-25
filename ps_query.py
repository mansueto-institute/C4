states = """SELECT fips AS id, ST_Simplify(ST_Scale(ST_Transform(geom, epsg), 0.001, 0.001, 0.001), 5) AS state 
            FROM states WHERE usps = upper('{}');"""

tracts = """SELECT census_tracts_2015.county, census_tracts_2015.tract, 
                   ST_Scale(ST_Transform(census_tracts_2015.geomsimp, states.epsg), 0.001, 0.001, 0.001) as geometry, 
                   ST_X(ST_Transform(census_tracts_2015.centroid, states.epsg))/1000 as x, 
                   ST_Y(ST_Transform(census_tracts_2015.centroid, states.epsg))/1000 as y,
                   census_tracts_2015.aland/1e6 as A, b01001_001e as pop,
                   CASE WHEN (ST_NumGeometries(census_tracts_2015.geomsimp) > 1) THEN 1 ELSE 0 END AS split
            FROM census_tracts_2015
            JOIN acssf5y2015 on
              census_tracts_2015.state = acssf5y2015.state AND
              census_tracts_2015.county = acssf5y2015.county AND
              census_tracts_2015.tract = acssf5y2015.tract
            JOIN states ON fips = census_tracts_2015.state 
            WHERE census_tracts_2015.state = {}
            ORDER BY county, tract;
          """

edges  = """
         SELECT
           rn,
           seq, edge_id eid, edge < 0 AS rev,
           CASE WHEN (edge < 0) THEN start_node ELSE end_node END AS nodeA, 
           CASE WHEN (edge > 0) THEN start_node ELSE end_node END AS nodeB,
           ST_Scale(ST_Transform(edge.geom, states.epsg), 0.001, 0.001, 0.001) as lines
         FROM
           census_tracts_2015  AS tab,
           state_topo.relation AS rel,
           state_topo.edge     AS edge,
           states,
           topology.ST_GetFaceEdges('state_topo', rel.element_id) AS t(seq,edge),
           (SELECT
              state, county, tract,
              row_number() over (PARTITION BY state ORDER BY county, tract NULLS LAST) - 1 AS rn
            FROM census_tracts_2015) rn
         WHERE
           states.fips = tab.state AND
           ABS(t.edge) = edge.edge_id AND
           tab.topogeom = (topology.GetTopologyID('state_topo'), 1, rel.topogeo_id, 3) AND
           rn.state = tab.state AND rn.county = tab.county AND rn.tract = tab.tract AND
           tab.state = {}
         ORDER BY
           rn, seq
         ;
         """

nodes = """
        SELECT
          DISTINCT ON(node_id, nseq, eid)
          node_id nid, ne.seq as nseq, abs(ne.edge) as eid,
          ST_X(ST_Transform(node.geom, states.epsg))/1000. x,
          ST_Y(ST_Transform(node.geom, states.epsg))/1000. y
        FROM
          census_tracts_2015  AS tab,
          state_topo.relation AS rel,
          state_topo.edge     AS edge,
          state_topo.node     AS node,
          topology.ST_GetFaceEdges('state_topo', rel.element_id) AS t(seq,edge),
          topology.GetNodeEdges('state_topo', node_id) as ne(seq,edge),
          states
        WHERE
          ABS(t.edge) = edge.edge_id AND
          (node_id = edge.start_node OR node_id = edge.end_node) AND
          tab.topogeom = (topology.GetTopologyID('state_topo'), 1, rel.topogeo_id, 3) AND
          states.fips = tab.state AND
          tab.state = {}
        ORDER BY node_id, nseq
        ;
        """


