import os
from netrc import netrc

if os.getenv("GMD_PASSWD") and os.getenv("GMD_USER"):
  user, passwd = os.getenv("GMD_USER"), os.getenv("GMD_PASSWD")
else: user, acct, passwd = netrc().authenticators("harris")

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()

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

shapefile = "shapes/{}.shp"
stateinfo = "shapes/{}_info"
edge_file = "shapes/{}_edges"
node_file = "shapes/{}_nodes"

def ens_dir(f, quiet = False):
  if not os.path.isdir(f):
    os.makedirs(f)
    print("Remade file", f)

def ens_data(usps): 

  ens_dir("shapes/")
  cache_stateinfo(usps)
  cache_shapefile(usps)
  cache_edge_file(usps)
  cache_node_file(usps)


import pycluscious as pycl
pycl_methods = {"reock"       : pycl.ObjectiveMethod.REOCK, 
                "dist_a"      : pycl.ObjectiveMethod.DISTANCE_A,
                "dist_p"      : pycl.ObjectiveMethod.DISTANCE_P,
                "inertia_a"   : pycl.ObjectiveMethod.INERTIA_A,
                "inertia_p"   : pycl.ObjectiveMethod.INERTIA_P,
                "polsby"      : pycl.ObjectiveMethod.POLSBY,
                "hull_a"      : pycl.ObjectiveMethod.HULL_A,
                "hull_p"      : pycl.ObjectiveMethod.HULL_P,
                "path_frac"   : pycl.ObjectiveMethod.PATH_FRAC,
                "ehrenburg"   : pycl.ObjectiveMethod.EHRENBURG,
                "axis_ratio"  : pycl.ObjectiveMethod.AXIS_RATIO,
                "mean_radius" : pycl.ObjectiveMethod.MEAN_RADIUS,
                "dyn_radius"  : pycl.ObjectiveMethod.DYN_RADIUS,
                "harm_radius" : pycl.ObjectiveMethod.HARM_RADIUS,
                "rohrbach"    : pycl.ObjectiveMethod.ROHRBACH,
                "exchange"    : pycl.ObjectiveMethod.EXCHANGE}

pycl_formal  = {
                "dist_a"      : "Distance to Areal Center",
                "dist_p"      : "Distance to Pop. Center",
                "reock"       : "Circumscribing Circles",
                "inertia_a"   : "Moment of Inertia: Area",
                "inertia_p"   : "Moment of Inertia: Pop.",
                "polsby"      : "Isoperimeter Quotient",
                "hull_a"      : "Convex Hull Area Ratio",
                "hull_p"      : "Convex Hull Pop. Ratio",
                "ehrenburg"   : "Inscribed Circles",
                "axis_ratio"  : "Length/Width Ratio",
                "mean_radius" : "Mean Radius",
                "dyn_radius"  : "Dynamic Radius",
                "harm_radius" : "Harmonic Radius",
                "rohrbach"    : "Distance to Perimeter",
                "exchange"    : "Exchange",
                "path_frac"   : "Path Fraction",
               }


pycl_circles = {"area"     : pycl.RadiusType.EQUAL_AREA, 
                "area_pop" : pycl.RadiusType.EQUAL_AREA_POP, 
                "circ"     : pycl.RadiusType.EQUAL_CIRCUMFERENCE, 
                "scc"      : pycl.RadiusType.SCC, 
                "lic"      : pycl.RadiusType.LIC,
                "hull"     : pycl.RadiusType.HULL,
                "power"    : pycl.RadiusType.POWER}


us_states = ["al", "ak", "az", "ar", "ca", "co", "ct", "de", "dc",
             "fl", "ga", "hi", "id", "il", "in", "ia", "ks", "ky",
             "la", "me", "md", "ma", "mi", "mn", "ms", "mo", "mt",
             "ne", "nv", "nh", "nj", "nm", "ny", "nc", "nd", "oh",
             "ok", "or", "pa", "ri", "sc", "sd", "tn", "tx", "ut",
             "vt", "va", "wa", "wv", "wi", "wy"]


def get_epsg (usps): return get_state_info(usps)["epsg"]
def get_fips (usps): return get_state_info(usps)["fips"]
def get_seats(usps): return get_state_info(usps)["seats"]

def get_state_info(usps):
   
   return pd.read_csv(stateinfo.format(usps) + ".csv").ix[0]


def cache_stateinfo(usps, filename = None):

   if not filename:
      filename = stateinfo.format(usps.lower())

   if os.path.exists(filename + ".csv"): return

   con = psycopg2.connect(database = "census", user = user, password = passwd,
                          host = "saxon.harris.uchicago.edu", port = 5432)

   info = pd.read_sql("SELECT epsg, fips, seats FROM states WHERE usps = upper('{}');".format(usps), con)
   info.to_csv(filename + ".csv", index = False)


   state_df = gpd.GeoDataFrame.from_postgis(ps_query.states.format(usps), con,
                                            geom_col='state', crs = from_epsg(get_epsg(usps)))

   state_df[["id", "state"]].to_file(filename + ".shp")


def cache_shapefile(usps, filename = None):

   if not filename:
      filename = shapefile.format(usps.lower())

   if os.path.exists(filename): return

   state = get_state_info(usps)

   con = psycopg2.connect(database = "census", user = user, password = passwd,
                          host = "saxon.harris.uchicago.edu", port = 5432)

   geo_df = gpd.GeoDataFrame.from_postgis(ps_query.tracts.format(state["fips"]), con,
                                          geom_col='geometry', crs = from_epsg(state["epsg"]))

   geo_df["id"] = geo_df.index
   geo_df[["id", "a", "pop", "x", "y", "split", "geometry"]].to_file(filename)



def cache_edge_file(usps, filename = None):

   if not filename:
      filename = edge_file.format(usps.lower())

   if os.path.exists(filename + ".csv"): return

   state = get_state_info(usps)

   con = psycopg2.connect(database = "census", user = user, password = passwd,
                          host = "saxon.harris.uchicago.edu", port = 5432)

   geo_df = gpd.GeoDataFrame.from_postgis(ps_query.edges.format(state["fips"]), con,
                                          geom_col='lines', crs = from_epsg(state["epsg"]))

   geo_df[geo_df["rev"] == False][["eid", "lines"]].to_file(filename + ".shp")

   geo_df[["rn","seq","eid","rev","nodea","nodeb"]].to_csv(filename + ".csv", index = False)


def cache_node_file(usps, filename = None):

   if not filename:
      filename = node_file.format(usps.lower())

   if os.path.exists(filename + ".csv"): return

   state = get_state_info(usps)

   con = psycopg2.connect(database = "census", user = user, password = passwd,
                          host = "saxon.harris.uchicago.edu", port = 5432)

   ndf = pd.read_sql(ps_query.nodes.format(state["fips"]), con)

   ndf[["nid", "x", "y", "nseq", "eid"]].to_csv(filename + ".csv", index = False)

   geometry = [Point(xy) for xy in zip(ndf.x, ndf.y)]
   geo_ndf = gpd.GeoDataFrame(ndf, crs = from_epsg(state["epsg"]), geometry=geometry)

   geo_ndf[ndf["nseq"] == 1][["nid", "geometry"]].to_file(filename + ".shp")




def plot_map(gdf, filename, crm, hlt = None, shading = "district", figsize = 10, label = "", ring = None, circ = None, point = None, legend = False):

    gdf["C"] = pd.Series(crm)

    if hlt:
      gdf["H"] = 0
      gdf.loc[hlt, "H"] = 1

    if shading == "density":
      gdf["density"] = gdf["pop"]/gdf["a"]

    dis = gdf.dissolve("C", aggfunc='sum')
    dis.reset_index(inplace = True)

    target = dis["pop"].sum() / dis.shape[0]
    dis["frac"] = dis["pop"] / target

    bounds = gdf.total_bounds
    xr = bounds[2] - bounds[0]
    yr = bounds[3] - bounds[1]
    fs = (figsize * np.sqrt(xr/yr), figsize * np.sqrt(yr/xr))

    bins = min(5, dis.shape[0])
    q = ps.Quantiles(dis["frac"], k = bins)

    labels = ["≤ {:.3f}".format(q.bins[i]) for i in range(bins) ]
    if len(set(labels)) < len(labels):
      labels = ["≤ {:.5f}".format(q.bins[i]) for i in range(bins) ]

    dis["qyb"] = q.yb
    dis["qyb"] = pd.Series(dis["qyb"], dtype="category")
    if len(set(labels)) == len(labels):
      dis["qyb"].cat.rename_categories(labels, inplace = True)

    if "target" in shading:

      col, alpha, trunc = "coolwarm", 0.7, ""
    
      if dis["frac"].max() > 2: 
          norm = Normalize(vmin = 0, vmax = 2)
          trunc = " (Truncated)"
      elif dis["frac"].max() - 1 < 0.005:
          norm = Normalize(vmin = 0.995, vmax = 1.005)
      else: # regardless, keep it centered
          larger = max(1 - dis["frac"].min(), dis["frac"].max() - 1)
          norm = Normalize(vmin = 1 - larger, vmax = 1 + larger) 
    
      cmap = plt.cm.ScalarMappable(norm=norm, cmap = col)
      
      ax = dis.plot(color = "white", edgecolor = "white", figsize = fs)
      for xi, row in dis.iterrows(): dis[dis.index == xi].plot(ax = ax, alpha = alpha, facecolor = cmap.to_rgba(row["frac"]))

      fig = ax.get_figure()
      cax = fig.add_axes([0.16, 0.13, 0.70, 0.015 * np.sqrt(xr/yr)])
      sm = plt.cm.ScalarMappable(cmap = col, norm=norm)
      sm._A = [] # gross
    
      cb = fig.colorbar(sm, cax = cax, alpha = alpha, # label = "Population / Target" + trunc, labelsize=12,
                        orientation='horizontal', drawedges = True)
      cb.locator = ticker.MaxNLocator(nbins=5)
      cb.formatter.set_useOffset(False)
      cb.set_label("Population / Target" + trunc, size=12)
      cb.ax.tick_params(labelsize=12)
      cb.dividers.set_visible(False)
      cb.update_ticks()

      # ax = dis.plot(column = "qyb", alpha = 0.3, categorical = True, cmap = "Greys", linewidth = 1.5, legend = True, figsize = fs)

      # if hlt: gdf[gdf["H"] == 1].plot(facecolor = "red", alpha = 0.1, linewidth = 0.05, ax = ax)

    elif "density" in shading:
      ax = gdf.plot(column = "density", cmap = "gray", scheme = "quantiles", k = 9, alpha = 0.8, figsize = fs, linewidth = 0)

      dis.plot(color = "blue", alpha = 0.3, linewidth = 1, ax = ax)

    else:
      ax = dis.plot("C", alpha = 0.5, categorical = True, cmap = "nipy_spectral", linewidth = 1, legend = legend, figsize = fs)

      if legend: ax.get_legend().set_bbox_to_anchor((1, 1))

      if hlt: gdf[gdf["H"] == 1].plot(facecolor = "grey", alpha = 0.1, linewidth = 0.05, ax = ax)

    ax.set_xlim([bounds[0] - 0.1 * xr, bounds[2] + 0.1 * xr])
    ax.set_ylim([bounds[1] - 0.1 * yr, bounds[3] + 0.1 * yr])

    if label: ax.text(bounds[0] - 0.16*xr, bounds[3] + 0.12*yr, label, fontsize = 10)

    ax.set_axis_off()

    if ring is not None:
      ring["C"] = ring.index
      if shading == "district":
        ring.plot("C", categorical = True, cmap = "nipy_spectral",  ax = ax, linewidth = 3)
        ring.plot(color = "white", ax = ax, linewidth = 1)
      else:
        ring.plot(color = "black", ax = ax, linewidth = 2.5)
        ring.plot(color = "white", ax = ax, linewidth = 0.7)


    if circ is not None:
      # circ["C"] = circ.index
      circ.plot(color = "white", alpha = 0.2, ax = ax, linewidth = 1)

    if point is not None:
      if "district" in shading:
        point["C"] = point.index
        point.plot("C", categorical = True, cmap = "nipy_spectral", ax = ax, markersize = 3)
      else:
        point.plot(color = "black", ax = ax, markersize = 3)
      point.plot(color = "white", ax = ax, markersize = 1)

    if not filename: return ax


    ax.figure.savefig(filename, bbox_inches='tight', pad_inches=0.05)
    plt.close('all')


def save_geojson(gdf, filename, crm, usps = None, metrics = None):

    gdf["C"] = pd.Series(crm)

    elections, vote_columns = [], []
    if usps and os.path.exists("el/{}_votes.csv".format(usps.lower())):
      votes = pd.read_csv("el/{}_votes.csv".format(usps.lower()), index_col = "rn")
      votes = votes.filter(regex = '[DR][0-9]{2}', axis = 1)

      gdf = gdf.join(votes)
      vote_columns = list(votes.columns)
      elections = set([el[1:] for el in votes.columns])

    dis = gdf.dissolve("C", aggfunc='sum')
    dis.reset_index(inplace = True)

    target = dis["pop"].sum() / dis.shape[0]
    dis["frac"] = (dis["pop"] / target).map('{:.03f}'.format)

    dis = dis[["C", "a", "frac", "geometry", "pop"] + vote_columns]

    dis["a"] *= 3.8610216e-7 # m2 to mi2
    dis["a"] = dis["a"].astype(int)
    dis.rename(columns = {"C" : "ID", "a" : "Area [sq mi]", "pop" : "Population", "frac" : "Pop./Target"}, inplace = True)

    for k, v in metrics.items():
      dis[k] = pd.Series(v).map('{:.03f}'.format)

    for el in elections:
      dis["Party " + el] = "R"
      dis.loc[dis["D" + el] > dis["R" + el], "Party " + el] = "D"

      dis["D" + el + " Share"] = (dis["D" + el] / (dis["D" + el] + dis["R" + el])).map('{:.03f}'.format) 
      dis["R" + el + " Share"] = (dis["R" + el] / (dis["D" + el] + dis["R" + el])).map('{:.03f}'.format) 

    if elections: dis.rename(columns = {v : v + " Votes" for v in vote_columns}, inplace = True)

    dis.crs = gdf.crs
    dis = dis.to_crs(epsg = 4326)
    
    fill = {}
    ndistricts = dis.shape[0]
    for ndi in range(ndistricts):
        color = [int(v*255) for v in plt.get_cmap("nipy_spectral")(ndi/ndistricts)][:3]
        color_hex = "#{0:02X}{1:02X}{2:02X}".format(*color)
        fill[ndi] = color_hex
        
    dis["fill"] = pd.Series(fill)
    dis["stroke-width"] = 2
    dis["stroke"] = "#000000"
    dis["fill-opacity"] = 0.1
    with open(filename, "w") as out: out.write(dis.to_json())
    # dis.to_file(filename, driver='GeoJSON') # crashes.


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

        edges[i] = a.length - sum(new_weights[i]) # /a.length 

    return edges, ps.W(Wsrc.neighbors, new_weights)



