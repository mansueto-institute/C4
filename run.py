#!/usr/bin/env python3

import argparse

import sys, os, re, math
from random import randint, choice

import pyc4
from pyc4_helper import *

import random
import operator

import time

import warnings
if not sys.warnoptions:
    warnings.simplefilter("ignore")


def load_data(state, method, seats = None, bgroup = False, scale_regimes = 1, scale_regime_perimeters = 1):

  ens_data(state, bgroup)

  if not seats: seats = get_seats(state) 

  if seats == 1:
    print(state, "has a single districting -- returning.")
    sys.exit()

  u = pyc4.universe(seats)

  tag = ""
  if bgroup: tag = "_bgroup"
  
  gdf = gpd.read_file(shapefile.format(state + tag))
  edges, spw = spw_from_shapefile(shapefile.format(state + tag))
  gdf["ps_n"] = pd.Series(spw.neighbors)
  gdf["ps_w"] = pd.Series(spw.weights)
  gdf["edge"] = pd.Series(edges) > 1e-3
  gdf["edge_perim"] = pd.Series(edges)

  # Additional, explicit connections.
  if os.path.exists(conx_file.format(state + tag)):
    for line in open(conx_file.format(state + tag)):

      line = line.split("#")[0].strip()
      x1, x2 = (int(v) for v in line.split(","))
      gdf.loc[x1, "ps_n"].append(x2)
      gdf.loc[x1, "ps_w"].append(1e-6)
      gdf.loc[x2, "ps_n"].append(x1)
      gdf.loc[x2, "ps_w"].append(1e-6)

  ### Add the cells....
  for xi, c in gdf.iterrows():
  
    ##  if c.edge:
    ##    print("Neighbors are ::", {n:w for n, w in zip(c.ps_n, c.ps_w)})
    ##    print("Edge perimeter is ::", c.edge_perim)
    u.add_cell(pyc4.cell(xi, int(c["pop"]), c.x, c.y, c.a, 
                         {n:w for n, w in zip(c.ps_n, c.ps_w)},
                         c.edge_perim, c.split))

  # Regime weight and node assignments, for weighted polsby
  if math.fabs(scale_regimes - 1) > 1e-5 or \
     math.fabs(scale_regime_perimeters -1 ) > 1e-5:
    u.add_regime("county", gdf.county.to_dict(), 
                 scale_regimes, scale_regime_perimeters)
  
  ### Add the edges....
  edf = pd.read_csv(edge_file.format(state + tag) + ".csv")
  rn, first, skip_rest = -1, -1, False
  for ei, e in edf.iterrows():

    u.add_edge(e.rn, e.eid, e.nodea, e.nodeb)
    

  ### Add the nodes....
  ndf = pd.read_csv(node_file.format(state + tag) + ".csv")
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

def circle_df(u, circle = False):

  if not circle: return None

  if circle == "hull":

    return gpd.GeoDataFrame(geometry=[Polygon(u.hull(d, True)) for d in range(u.nregions)])

  return gpd.GeoDataFrame(geometry=[Point(c[0][0], c[0][1]).buffer(c[1] if math.isfinite(c[1]) and c[1] > 0 else 1)
                                    for c in [u.get_circle_coords(r, pyc4_circles[circle])
                                              for r in range(u.nregions)]])

def point_df(u, point = True):

  if not point: return None

  return gpd.GeoDataFrame(geometry=[Point(c[0][0], c[0][1])
                                    for c in [u.get_circle_coords(r, pyc4_circles[point])
                                    for r in range(u.nregions)]])


def main(state, seed, method, seats, ncycles, split_restart, power_restart, niter, nloops, tol, ctol, conv_iter, init, output, write, 
         grasp, allow_trades, destrand_inputs, destrand_min, destrand_max, tabu_length,
         circle, ring, point, print_init, no_plot, shading, borders, verbose, bgroup, 
         scale_regimes, scale_regime_perimeters):

  u, gdf = load_data(state, method, seats, bgroup, scale_regimes, scale_regime_perimeters)

  u.RANDOM       = grasp 
  u.TRADE        = allow_trades
  u.TABU_LENGTH  = u.nregions * tabu_length
  u.DESTRAND_MIN = destrand_min
  u.DESTRAND_MAX = destrand_max

  if "seed" in init:
    u.rand_init(seed)

  elif "kmeans" in init:
    u.rand_init(seed)
    u.grow_kmeans(method == "dist_p") # True is population growing

  elif "power" in init:
    u.rand_init(seed)
    u.grow_kmeans()

    sinit = init.split(":")
    npiter = int(sinit[1]) if len(init) > 1 else 10000
    u.iterate_power(tol, npiter, 1, verbose)

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

  ens_dir("{}/{}".format(output, write))
  ens_dir("{}/json/".format(output))

  for c in range(ncycles):

    if ncycles > 1:
      write_cycle = write + "/c{:03d}".format(c)
      ens_dir("{}/{}".format(output, write_cycle))
      if c:
        if power_restart:
          sinit = init.split(":")
          npiter = int(sinit[1]) if len(init) > 1 else 10000
          u.power_restart(seed + c * 1000, npiter, tol, verbose)

        if split_restart: u.split_restart(seed+c, pyc4_methods[method])

    else: write_cycle = write

    for i in range(0, nloops+1):

      converged = False
      if i:
        converged = u.oiterate(pyc4_methods[method], niter = niter, llh_tol = tol, cut_tol = ctol, conv_iter = conv_iter, seed = seed, verbose = verbose)
      elif not print_init: continue
      
      if (converged or i == nloops) and u.get_best_solution(): u.load_best_solution()

      crm = u.cell_region_map()

      tag = "final" if (i == nloops or converged) else "i{:03d}".format(i)

      for s in shading:

        # Score
        scores = None
        if s == "scores":
            if method not in pyc4_methods: continue
            scores = u.get_objectives(pyc4_methods[method])

        plot_map(gdf, "{}/{}/{}_{}.pdf".format(output, write_cycle, tag, s),
                 label = pyc4_formal[method] if i else init.capitalize(),
                 crm = crm, hlt = u.border_cells(True if "ext" in borders else False) if borders else None, shading = s,
                 ring = ring_df(u, (ring or s == "density")),
                 circ = circle_df(u, circle), point = point_df(u, point), scores = scores, legend = bool(verbose))

      with open ("{}/{}/{}.csv".format(output, write_cycle, tag), "w") as out:
        for k, v in crm.items(): out.write("{},{}\n".format(k, v))

      print(write_cycle, ":: completed iteration", i, flush = True)

      if converged: break

    save_json("{}/json/{}.json".format(output, write_cycle.replace("/", "_")),
              state, pyc4_short[method], write_cycle, gdf, crm = u.cell_region_map(),
              metrics = {pyc4_short[k] : u.get_objectives(v) for k, v in pyc4_methods.items()},
              seats = u.nregions, bgroup = bgroup)

    save_geojson(gdf, "{}/{}/final.geojson".format(output, write_cycle), u.cell_region_map(), state,
                 metrics = {pyc4_short[k] : u.get_objectives(v) for k, v in pyc4_methods.items()},
                 bgroup = bgroup)


if __name__ == "__main__":


  parser = argparse.ArgumentParser()

  # Initialization 
  parser.add_argument("-i", "--init",      default = "seed", type = str)
  parser.add_argument("-x", "--seed",      default = 0, type = int)
  parser.add_argument("-s", "--state",     default = "pa", type=str.lower, choices = us_states, help='state')
  parser.add_argument("-w", "--write",     default = "", type = str)
  parser.add_argument("-m", "--method",    default = "dist_a", choices = pyc4_formal, type = str)
  parser.add_argument("-b", "--bgroup",    action  = "store_true")

  parser.add_argument("--scale_regimes",           default = 1, type = float)
  parser.add_argument("--scale_regime_perimeters", default = 1, type = float)

  # Looping parameters.
  parser.add_argument("-c", "--ncycles",   default = 1, type = int, help = "Number of restarts (split/merge)")
  parser.add_argument("-l", "--nloops",    default = 1, type = int, help = "Loops: number of times to run iter")
  parser.add_argument("-n", "--niter",     default = 100, type = int, help = "Max iterations per loop.")
  parser.add_argument("--conv_iter",       default = 0, type = int, help = "Stop after X without improvement.")
  parser.add_argument("-t", "--tol",       default = 0.01, type = float)
  parser.add_argument("--ctol",            default = 0.02, type = float)
  parser.add_argument("--power_restart",   action  = "store_true", help = "Number of restarts (split/merge)")
  parser.add_argument("--split_restart",   action  = "store_true", help = "Number of restarts (split/merge)")

  # Metaheuristic and minimum configuration
  parser.add_argument("--grasp",           action  = "store_true")
  parser.add_argument("--tabu_length",     default = 0, type = int)
  parser.add_argument("--destrand_min",    default = 2, type = int)
  parser.add_argument("--destrand_max",    default = 0, type = int)
  parser.add_argument("--destrand_inputs", action  = "store_true")
  parser.add_argument("--allow_trades",    action  = "store_true")

  # Plotting options.
  parser.add_argument("-r", "--ring",      action  = "store_true")
  parser.add_argument("-o", "--circle",    default = "", choices = pyc4_circles, type = str)
  parser.add_argument("-p", "--point",     default = None, choices = pyc4_circles, type = str)
  parser.add_argument("--shading",         default = ["target"], nargs = "+")
  parser.add_argument("--no_plot",         action  = "store_true")
  parser.add_argument("--borders",         default = "", choices = ["ext", "int"])
  parser.add_argument("--print_init",      action  = "store_true")
  parser.add_argument("--output",          default = "res")

  # Number of seats -- default is US Cong.
  parser.add_argument("--seats",           default = 0, type = int)

  # Verbosity
  parser.add_argument("-v", "--verbose",   default = 0, type = int)
  args = parser.parse_args()

  if not args.write: args.write = "{}/{}/s{:03d}".format(args.state, args.method, args.seed)
  if "all" in [s.lower() for s in args.shading]: args.shading = ["district", "target", "density", "scores", "counties"]
  if "none" in [s.lower() for s in args.shading]: args.shading = []
  if args.nloops == 0: method = args.init

  if args.ncycles > 1 and not args.power_restart: args.split_restart = True

  if args.no_plot: shading = []

  # print(vars(args))
  main(**vars(args))



