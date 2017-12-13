states = """SELECT fips AS id, ST_Simplify(ST_Transform(geom, epsg), 5) AS state 
            FROM states WHERE usps = upper('{}');"""

shapes = """SELECT census_bg_2010.county, census_bg_2010.tract, census_bg_2010.bgroup,
                   ST_Transform(census_bg_2010.geomsimp, states.epsg) as geometry, 
                   ST_X(ST_Transform(census_bg_2010.centroid, states.epsg)) as x, 
                   ST_Y(ST_Transform(census_bg_2010.centroid, states.epsg)) as y,
                   ST_Area(ST_Transform(census_bg_2010.geom, states.epsg))::int AS A, pop,
                   CASE WHEN (ST_NumGeometries(census_bg_2010.geomsimp) > 1) THEN 1 ELSE 0 END AS split
            FROM census_bg_2010
            JOIN states ON fips = census_bg_2010.state 
            WHERE census_bg_2010.state = {}
            ORDER BY county, tract, bgroup;
         """

edges  = """
         SELECT
           rn,
           seq, edge_id eid, edge < 0 AS rev,
           CASE WHEN (edge < 0) THEN start_node ELSE end_node END AS nodeA, 
           CASE WHEN (edge > 0) THEN start_node ELSE end_node END AS nodeB,
           ST_Transform(edge.geom, states.epsg) as lines
         FROM
           census_bg_2010  AS tab,
           topo_bg_2010.relation AS rel,
           topo_bg_2010.edge     AS edge,
           states,
           topology.ST_GetFaceEdges('topo_bg_2010', rel.element_id) AS t(seq,edge),
           (SELECT
              state, county, tract, bgroup,
              row_number() over (PARTITION BY state ORDER BY county, tract, bgroup NULLS LAST) - 1 AS rn
            FROM census_bg_2010) rn
         WHERE
           states.fips = tab.state AND
           ABS(t.edge) = edge.edge_id AND
           tab.topogeom = (topology.GetTopologyID('topo_bg_2010'), 1, rel.topogeo_id, 3) AND
           rn.state = tab.state AND rn.county = tab.county AND rn.tract = tab.tract AND rn.bgroup = tab.bgroup AND
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
          census_bg_2010  AS tab,
          topo_bg_2010.relation AS rel,
          topo_bg_2010.edge     AS edge,
          topo_bg_2010.node     AS node,
          topology.ST_GetFaceEdges('topo_bg_2010', rel.element_id) AS t(seq,edge),
          topology.GetNodeEdges('topo_bg_2010', node_id) as ne(seq,edge),
          states
        WHERE
          ABS(t.edge) = edge.edge_id AND
          (node_id = edge.start_node OR node_id = edge.end_node) AND
          tab.topogeom = (topology.GetTopologyID('topo_bg_2010'), 1, rel.topogeo_id, 3) AND
          states.fips = tab.state AND
          tab.state = {}
        ORDER BY node_id, nseq
        ;
        """

race = """
       SELECT
         row_number() over
            (PARTITION BY state ORDER BY county, tract NULLS LAST) - 1 AS rn,
         state, county cid, tract, bgroup,
         pop, black, hispanic,
         vap total_vap,
         bvap black_vap,
         hvap hispanic_vap
       FROM census_bg_2010 AS bg
       JOIN states AS st ON
         bg.state  = st.fips
       WHERE
         st.usps = UPPER('{}')
       ;
       """


