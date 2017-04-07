#!/usr/bin/env python 

import argparse

import sys, os, re
from random import randint

import pycluscious as pycl
from pycluscious_helper import *

import random

def load_data(state):

  ens_data(state)

  u = pycl.universe(get_seats(state))
  
  gdf = gpd.read_file(shapefile.format(state))
  edges, spw = spw_from_shapefile(shapefile.format(state))
  gdf["ps_n"] = pd.Series(spw.neighbors)
  gdf["ps_w"] = pd.Series(spw.weights)
  gdf["edge"] = pd.Series(edges) < 0.99

  ### Add the cells....
  for xi, c in gdf.iterrows():
  
    u.add_cell(pycl.cell(xi, int(c["pop"]), c.x, c.y, c.a, 
                         {n:w for n, w in zip(c.ps_n, c.ps_w)}, c.edge))
  
  
  ### Add the edges....
  edf = pd.read_csv(edge_file.format(state))
  rn, first, skip_rest = -1, -1, False
  for ei, e in edf.iterrows():

    u.add_edge(e.rn, e.eid, e.nodea, e.nodeb)
  
  
  ### Add the nodes....
  ndf = pd.read_csv(node_file.format(state))
  node = -1
  for ni, n in ndf.iterrows():
  
    if n.nid != node:
      u.add_node(n.nid, n.x, n.y)
      node = n.nid
  
    u.add_node_edge(n.nid, n.eid)
  
  ### Now transform everything into pointers.
  u.adjacency_to_pointers()
  u.node_ids_to_pointers()
  
  ### Clean up the graph...
  u.connect_graph()
  u.trim_graph()
  
  u.build_dijkstra_graph()

  return u, gdf


def main(state, seed, method, niter, nloops, load, write, ring, pop, verbose):

  u, gdf = load_data(state)

  if not load:

    u.rand_init(seed)
    u.grow_kmeans(method == "dist_p") # True is population growing
  else:
    u.load_partition(load)
  
  # Print once before looping.
  plot_map(gdf, "results/{}_i000.pdf".format(write),
           crm = u.cell_region_map(), hlt = u.border_cells(), pop = pop)


  for i in range(1, nloops+1):

    u.oiterate(pycl_methods[method], niter = niter, tol = 0.05, alpha = 4, verbose = False)
    
    crm = u.cell_region_map()
    
    if ring: ring_df = gpd.GeoDataFrame(geometry=[LineString(u.get_point_ring(d)) for d in range(u.nregions)])
    else: ring_df = None

    for gi in range(5):
      for r in range(u.nregions):
        u.add_cell_to_region(r, random.choice(u.border_cells(True, r)))

      ring_df = gpd.GeoDataFrame(geometry=[LineString(u.get_point_ring(d)) for d in range(u.nregions)])

      plot_map(gdf, "results/{}_i{:03d}_g{:03d}.pdf".format(write, i, gi),
               crm = u.cell_region_map(), hlt = u.border_cells(), pop = pop, ring = ring_df)
    
    with open ("results/{}_i{:03d}.csv".format(write, i), "w") as out:
      for k, v in crm.items(): out.write("{},{}\n".format(k, v))
    
    print("completed iteration ::", i)



if __name__ == "__main__":


  parser = argparse.ArgumentParser()
  parser.add_argument("-i", "--seed",    default = 0, type = int)
  parser.add_argument("-s", "--state",   default = "pa", type=str.lower, choices = us_states, help='state')
  parser.add_argument("-m", "--method",  default = "dist_a", choices = pycl_methods, type = str)
  parser.add_argument("-n", "--niter",   default = 100, type = int)
  parser.add_argument("-l", "--nloops",  default = 1, type = int)
  parser.add_argument("-f", "--load",    default = "", type = str)
  parser.add_argument("-w", "--write",   default = "", type = str)
  parser.add_argument("-p", "--pop",     action  = "store_true")
  parser.add_argument("-r", "--ring",    action  = "store_true")
  parser.add_argument("-v", "--verbose", default = 0, type = int)
  args = parser.parse_args()

  if not args.write: args.write = "{}_{}_s{:03d}".format(args.state, args.method, args.seed)

  main(**vars(args))



