states = """SELECT fips AS id, ST_Simplify(ST_Transform(geom, epsg), 5) AS state 
            FROM states WHERE usps = upper('{}');"""

shapes = """SELECT census_tracts_2015.county, census_tracts_2015.tract, 
                   ST_Transform(census_tracts_2015.geomsimp, states.epsg) as geometry, 
                   ST_X(ST_Transform(census_tracts_2015.centroid, states.epsg)) as x, 
                   ST_Y(ST_Transform(census_tracts_2015.centroid, states.epsg)) as y,
                   census_tracts_2015.area as A, b01001_001e as pop,
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
           ST_Transform(edge.geom, states.epsg) as lines
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
          ST_X(ST_Transform(node.geom, states.epsg)) x,
          ST_Y(ST_Transform(node.geom, states.epsg)) y
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

race = """SELECT
            rn.rn, d.state, d.county cid, d.tract,
            b01001_001e pop, b01001a_001e white,
            b01001b_001e black, b01001i_001e hispanic,
            total_vap, black_vap, hispanic_vap
          FROM census_tracts_2015 AS g
          JOIN states AS s ON
            g.state = s.fips
          JOIN acssf5y2015 AS d ON
            d.state  = g.state  AND
            d.county = g.county AND
            d.tract  = g.tract
          JOIN (SELECT state, county, tract,
                       row_number() over (PARTITION BY state ORDER BY county, tract NULLS LAST) - 1 as rn
                FROM census_tracts_2015) rn ON
            g.state  = rn.state  AND
            g.county = rn.county AND
            g.tract  = rn.tract
          WHERE s.usps = UPPER('{}')
          ORDER BY d.state, d.county, d.tract
          ;
          """

