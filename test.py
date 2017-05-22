#!/usr/bin/env python3

import argparse

import sys, os, re, math
from random import randint, choice

import pycluscious as pycl
from pycluscious_helper import *

import random
import operator

import time

def load_data(state, method):

  ens_data(state)

  u = pycl.universe(get_seats(state))
  
  gdf = gpd.read_file(shapefile.format(state))
  edges, spw = spw_from_shapefile(shapefile.format(state))
  gdf["ps_n"] = pd.Series(spw.neighbors)
  gdf["ps_w"] = pd.Series(spw.weights)
  gdf["edge"] = pd.Series(edges) > 1e-3
  gdf["edge_perim"] = pd.Series(edges)

  ### Add the cells....
  for xi, c in gdf.iterrows():
  
    ##  if c.edge:
    ##    print("Neighbors are ::", {n:w for n, w in zip(c.ps_n, c.ps_w)})
    ##    print("Edge perimeter is ::", c.edge_perim)
    u.add_cell(pycl.cell(xi, int(c["pop"]), c.x, c.y, c.a, 
                         {n:w for n, w in zip(c.ps_n, c.ps_w)},
                         c.edge_perim, c.split))

  
  ### Add the edges....
  edf = pd.read_csv(edge_file.format(state) + ".csv")
  rn, first, skip_rest = -1, -1, False
  for ei, e in edf.iterrows():

    u.add_edge(e.rn, e.eid, e.nodea, e.nodeb)
  
  
  ### Add the nodes....
  ndf = pd.read_csv(node_file.format(state) + ".csv")
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
  print("Connecting and trimming graph.")
  u.connect_graph()
  u.trim_graph()
  u.build_dijkstra_graph()

  print("Finished preparing the data.")

  return u, gdf

def ring_df(u, ring = True):

  if not ring: return None

  return gpd.GeoDataFrame(geometry=[LineString(u.get_point_ring(d)) for d in range(u.nregions)])

def circle_df(u, circle = True):

  if not circle: return None

  return gpd.GeoDataFrame(geometry=[Point(c[0][0], c[0][1]).buffer(c[1] if math.isfinite(c[1]) and c[1] > 0 else 1)
                                    for c in [u.get_circle_coords(r, pycl_circles[circle])
                                    for r in range(u.nregions)]])

def point_df(u, point = True):

  if not point: return None

  return gpd.GeoDataFrame(geometry=[Point(c[0][0], c[0][1])
                                    for c in [u.get_circle_coords(r, pycl_circles[point])
                                    for r in range(u.nregions)]])


def main(state, seed, method, ncycles, niter, nloops, tol, conv_iter, init, write, 
         grasp, allow_trades, destrand_inputs, destrand_min, destrand_max, tabu_length,
         circle, ring, point, print_init, no_plot, shading, borders, verbose):

  u, gdf = load_data(state, method)

  if u.nregions == 1: return

  u.RANDOM       = grasp 
  u.TRADE        = allow_trades
  u.TABU_LENGTH  = u.nregions * tabu_length
  u.DESTRAND_MIN = destrand_min
  u.DESTRAND_MAX = destrand_max

  if not init or "kmeans" in init:
    u.rand_init(seed)
    u.grow_kmeans(method == "dist_p") # True is population growing

  elif "power" in init:
    u.rand_init(seed)
    u.grow_kmeans()

    sinit = init.split(":")
    npiter = int(sinit[1]) if len(init) > 1 else 10000
    u.iterate_power(tol, npiter, 1)

  elif "rand" in init:
    u.rand_init(seed)
    u.grow_random(seed)
    print("Finished random growth initialization.")

  elif "split" in init:
    u.assign_to_zero()
    u.split_line_init()

  elif "csv" in init:
    u.load_partition(init)

    while destrand_inputs:
      destrand_inputs = u.destrand(destrand_min, destrand_max)

  else: 
    print("Tried to initialize with unknown argument,", init)
    print("Exiting....")
    sys.exit()

  ens_dir("res/{}".format(write))

  for c in range(ncycles):

    if ncycles > 1:
      write_cycle = write + "/c{:03d}".format(c)
      ens_dir("res/{}".format(write_cycle))
      if c: u.reboot(seed, pycl_methods[method])
    else: write_cycle = write

    for i in range(0, nloops+1):

      if i: u.oiterate(pycl_methods[method], niter = niter, tol = tol, conv_iter = conv_iter, seed = seed, verbose = verbose)
      elif not print_init: continue
      
      if i == nloops and u.get_best_solution(): u.load_best_solution()

      crm = u.cell_region_map()

      for s in shading:
        style = "" if len(shading) == 1 else "_{}".format(s)
        plot_map(gdf, "res/{}/i{:03d}{}.pdf".format(write_cycle, i, style), label = pycl_formal[method] if i else init.capitalize(),
                 crm = crm, hlt = u.border_cells(True if "ext" in borders else False) if borders else None, shading = s,
                 ring = ring_df(u, (ring or s == "density")),
                 circ = circle_df(u, circle), point = point_df(u, point), legend = verbose)

      with open ("res/{}/i{:03d}.csv".format(write_cycle, i), "w") as out:
        for k, v in crm.items(): out.write("{},{}\n".format(k, v))

      print(write_cycle, ":: completed iteration", i)

    save_geojson(gdf, "res/{}/final.geojson".format(write_cycle), crm,
                 metrics = {v : u.get_objectives(pycl_methods[k]) for k, v in pycl_formal.items()})




if __name__ == "__main__":


  parser = argparse.ArgumentParser()

  # Initialization 
  parser.add_argument("-i", "--init",      default = "", type = str)
  parser.add_argument("-x", "--seed",      default = 0, type = int)
  parser.add_argument("-s", "--state",     default = "pa", type=str.lower, choices = us_states, help='state')
  parser.add_argument("-w", "--write",     default = "", type = str)
  parser.add_argument("-m", "--method",    default = "dist_a", choices = pycl_methods, type = str)

  # Looping parameters.
  parser.add_argument("-c", "--ncycles",   default = 1, type = int, help = "Number of reboots (split/merge)")
  parser.add_argument("-l", "--nloops",    default = 1, type = int, help = "Loops: number of times to run iter")
  parser.add_argument("-n", "--niter",     default = 100, type = int, help = "Max iterations per loop.")
  parser.add_argument("--conv_iter",       default = 0, type = int, help = "Stop after X without improvement.")
  parser.add_argument("-t", "--tol",       default = 0.02, type = float)

  # Metaheuristic and minimum configuration
  parser.add_argument("--grasp",           action  = "store_true")
  parser.add_argument("--tabu_length",     default = 0, type = int)
  parser.add_argument("--destrand_min",    default = 2, type = int)
  parser.add_argument("--destrand_max",    default = 0, type = int)
  parser.add_argument("--destrand_inputs", action  = "store_true")
  parser.add_argument("--allow_trades",    action  = "store_true")

  # Plotting options.
  parser.add_argument("-r", "--ring",      action  = "store_true")
  parser.add_argument("-o", "--circle",    default = "", choices = pycl_circles, type = str)
  parser.add_argument("-p", "--point",     default = None, choices = pycl_circles, type = str)
  parser.add_argument("--shading",         default = ["target"], nargs = "+")
  parser.add_argument("--no_plot",         action  = "store_true")
  parser.add_argument("--borders",         default = "", choices = ["ext", "int"])
  parser.add_argument("--print_init",      action  = "store_true")

  # Verbosity
  parser.add_argument("-v", "--verbose",   default = 0, type = int)
  args = parser.parse_args()

  if not args.write: args.write = "{}/{}/s{:03d}".format(args.state, args.method, args.seed)
  if "all" in args.shading: args.shading = ["district", "target", "density"]

  if args.no_plot: shading = []

  main(**vars(args))



