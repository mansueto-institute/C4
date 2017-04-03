#!/usr/bin/env python 

import sys, os, re

import pycluscious as pycl
from pycluscious_helper import *


usps, seed = "PA", 1
if len(sys.argv) > 1:
  usps = sys.argv[1]
if len(sys.argv) > 2:
  seed = int(sys.argv[2])


method = methods[seed % len(methods)]

ens_dir("shapes/")
cache_shapefile(usps)
cache_edge_file(usps)
cache_node_file(usps)

gdf = gpd.read_file(shapefile.format(usps))

edges, spw = spw_from_shapefile(shapefile.format(usps))
gdf["ps_n"] = pd.Series(spw.neighbors)
gdf["ps_w"] = pd.Series(spw.weights)
gdf["edge"] = pd.Series(edges) < 0.99


seats = get_seats(usps)
if seats == 1: sys.exit()
u = pycl.universe(seats)


for xi, c in gdf.iterrows():

  u.add_cell(pycl.cell(xi, int(c["pop"]), c.x, c.y, c.a, 
                       {n:w for n, w in zip(c.ps_n, c.ps_w)}, c.edge))


edf = pd.read_csv(edge_file.format(usps))
rn, nb = -1, 0
for ei, e in edf.iterrows():

  if rn != e.rn: rn = e.rn
  # elif nb != e.nodea: continue

  u.add_edge(e.rn, e.eid, e.nodea, e.nodeb)

  nb = e.nodeb


ndf = pd.read_csv(node_file.format(usps))
node = -1
for ni, n in ndf.iterrows():

  if n.nid != node:
    u.add_node(n.nid, n.x, n.y)
    node = n.nid

  u.add_node_edge(n.nid, n.eid)

u.adjacency_to_pointers()
u.node_ids_to_pointers()

u.connect_graph()
u.trim_graph()

u.rand_districts(seed)
u.grow_kmeans(method[0] == "dist_p") # False is population growing

# for r in range(seats):


for i in range(10):

  if not i: continue
  elif i == 1: 
    u.oiterate(method[1], niter = 5000, tol = 0.02, alpha = 4, verbose = False)
    # u.iterate(1000, 0.05)
  else: break

  crm = u.cell_region_map()
  bc = u.border_cells()
  # bc = list(gdf[gdf.edge].index)

  geo_ring = [LineString(u.get_ring(d)) for d in range(seats)]
  ring_df = gpd.GeoDataFrame(geometry=geo_ring)
  
  for col in ["cat"]: # , "qyb"]:
    plot_map(gdf, "results/{}_s{:03d}_i{:03d}_{}_{}".\
                  format(usps, seed, i, method[0], col), crm, bc, col,
             ring = ring_df)

  print("completed iteration ::", i)


