from netrc import netrc
user, acct, passwd = netrc().authenticators("harris")

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import ticker

from fiona.crs import from_epsg
import pysal as ps
import pandas as pd
import geopandas as gpd
import numpy as np
import psycopg2

import ps_query

# Needed to keep the centroids in the boundary.
import shapely
from shapely.wkt import loads, dumps
from shapely.geometry import MultiPolygon, Polygon, Point, MultiLineString, LineString, LinearRing, asShape

import os, glob

con = psycopg2.connect(database = "census", user = user, password = passwd,
                       host = "saxon.harris.uchicago.edu", port = 5432)

shapefile = "shapes/{}.shp"
edge_file = "shapes/{}_edges.csv"
node_file = "shapes/{}_nodes.csv"

def ens_dir(f, quiet = False):
  if not os.path.isdir(f):
    os.makedirs(f)
    print("Remade file", f)


import pycluscious as pycl
methods = [["reock"     , pycl.ObjectiveMethod.REOCK], 
           ["dist_a"    , pycl.ObjectiveMethod.DISTANCE_A],
           ["dist_p"    , pycl.ObjectiveMethod.DISTANCE_P],
           ["inertia_a" , pycl.ObjectiveMethod.INERTIA_A],
           ["inertia_p" , pycl.ObjectiveMethod.INERTIA_P],
           ["polsby"    , pycl.ObjectiveMethod.POLSBY],
           ["hull_a"    , pycl.ObjectiveMethod.HULL_A]]

def get_epsg (usps): return get_state_info(usps)["epsg"]
def get_fips (usps): return get_state_info(usps)["fips"]
def get_seats(usps): return get_state_info(usps)["seats"]

def get_state_info(usps):

   return pd.read_sql("SELECT * FROM states WHERE usps = upper('{}');".format(usps), con).ix[0]



def cache_shapefile(usps, filename = None):

   if not filename:
      filename = shapefile.format(usps.lower())

   if os.path.exists(filename): return

   state = get_state_info(usps)

   geo_df = gpd.GeoDataFrame.from_postgis(ps_query.tracts.format(state["fips"]), con,
                                          geom_col='geometry', crs = from_epsg(state["epsg"]))

   geo_df[["a", "pop", "x", "y", "geometry"]].to_file(filename)



def cache_edge_file(usps, filename = None):

   if not filename:
      filename = edge_file.format(usps.lower())

   if os.path.exists(filename): return

   pd.read_sql(ps_query.edges.format(get_fips(usps)), con).to_csv(filename, index = False)



def cache_node_file(usps, filename = None):

   if not filename:
      filename = node_file.format(usps.lower())

   if os.path.exists(filename): return

   pd.read_sql(ps_query.nodes.format(get_fips(usps)), con).to_csv(filename, index = False)



def plot_map(gdf, filename, crm, bc, col = "qyb", figsize = 10, label = "", ring = None):

    gdf["C"] = pd.Series(crm)
    gdf["B"] = 0
    gdf.loc[bc, "B"] = 1

    dis = gdf.dissolve("C", aggfunc='sum')
    dis.reset_index(inplace = True)

    target = dis["pop"].sum() / dis.shape[0]
    dis["frac"] = dis["pop"] / target

    bounds = gdf.total_bounds
    xr = bounds[2] - bounds[0]
    yr = bounds[3] - bounds[1]

    bins = min(5, dis.shape[0])
    q = ps.Quantiles(dis["frac"], k = bins)

    labels = ["≤ {:.3f}".format(q.bins[i]) for i in range(bins) ]

    dis["qyb"] = q.yb
    dis["qyb"] = pd.Series(dis["qyb"], dtype="category")
    if len(set(labels)) == len(labels):
      dis["qyb"].cat.rename_categories(labels, inplace = True)

    if col == "qyb": 
      alpha, fc = 0.3, "red"
      ax = dis.plot(column = "qyb", alpha = 0.3, categorical = True, cmap = "Greys",
                    legend = True, figsize = (figsize * np.sqrt(xr/yr), figsize * np.sqrt(yr/xr)))

    else: 
      ax = dis.plot("C", alpha = 0.5, categorical = True, cmap = "nipy_spectral", # color = "0.9",
                    legend = False, figsize = (figsize * np.sqrt(xr/yr), figsize * np.sqrt(yr/xr)))

      alpha, fc = 0.3, "grey"

    gdf[gdf["B"] == 1].plot(facecolor = fc, alpha = alpha, linewidth = 0.05, ax = ax)

    ax.set_xlim([bounds[0] - 0.1 * xr, bounds[2] + 0.1 * xr])
    ax.set_ylim([bounds[1] - 0.1 * yr, bounds[3] + 0.1 * yr])

    if label: ax.text(bounds[0] - 0.1*xr, bounds[1] - 0.1*yr, label, fontsize = 10)

    ax.set_axis_off()

    if not filename: return ax

    if ring is not None:
      ring["C"] = ring.index
      ring.plot("C", categorical = True, cmap = "nipy_spectral",  ax = ax, linewidth = 3)
      ring.plot(color = "white", ax = ax, linewidth = 1)

    ax.figure.savefig(filename + ".pdf", bbox_inches='tight', pad_inches=0.05)
    plt.close('all')



def fix_mp(poly):
    
    if type(poly) != shapely.geometry.multipolygon.MultiPolygonAdapter:
        return poly
      
    if poly.is_valid: return poly
    
    mp_out = MultiPolygon()
    for p in list(poly):
        mp_out |= Polygon(p.exterior.coords[:])

    for p in list(poly):
        for ir in p.interiors:
            mp_out -= Polygon(ir.coords[:])

    return mp_out
    

# using pysal and shapely; very slightly modified from the contrib:
# https://github.com/pysal/pysal/blob/master/pysal/contrib/shared_perimeter_weights.py
def spw_from_shapefile(shapefile, norm = False):
    polygons = ps.open(shapefile, 'r').read()
    spolygons = list(map(asShape,polygons))
    spolygons = [fix_mp(p) for p in spolygons]
    perimeters = [p.length if norm else 1. for p in spolygons]
    Wsrc = ps.queen_from_shapefile(shapefile)
    new_weights, edges = {}, {}
    for i in Wsrc.neighbors:
        a = spolygons[i]
        p = perimeters[i]
        new_weights[i] = [] 
        for j in Wsrc.neighbors[i]:

            intersect = a.intersection(spolygons[j])
            new_weights[i].append(intersect.length)

            # if type(intersect) is LineString:
            #   edges[i].append((j, intersect.length))
            # if type(intersect) is MultiLineString:
            #   for line in list(intersect):
            #     edges[i].append((j, line.length))

        edges[i] = sum(new_weights[i])/a.length 

    return edges, ps.W(Wsrc.neighbors, new_weights)



