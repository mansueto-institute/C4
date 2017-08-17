import os, sys
from netrc import netrc

if os.getenv("GMD_PASSWD") and os.getenv("GMD_USER"):
  user, passwd = os.getenv("GMD_USER"), os.getenv("GMD_PASSWD")
elif netrc().authenticators("harris"):
  user, acct, passwd = netrc().authenticators("harris")
else: user, acct, passwd = None, None, None

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


# Needed to keep the centroids in the boundary.
import shapely
from shapely.wkt import loads, dumps
from shapely.geometry import MultiPolygon, Polygon, Point, MultiLineString, LineString, LinearRing, asShape

import os, glob

import json

shapefile = "shapes/{}.shp"
stateinfo = "shapes/{}_info"
edge_file = "shapes/{}_edges"
node_file = "shapes/{}_nodes"
race_file = "demographic/{}_race.csv"
vote_file = "demographic/{}_votes.csv"

def ens_dir(f, quiet = False):
  if not os.path.isdir(f):
    os.makedirs(f)
    # print("Remade file", f)

def ens_data(usps, blocks): 

  ens_dir("shapes/")
  ens_dir("demographic/")
  cache_stateinfo(usps)
  cache_shapefile(usps, blocks)
  cache_edge_file(usps, blocks)
  cache_node_file(usps, blocks)
  cache_race_file(usps)


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
                "power"       : "Power Diagram"
               }


pycl_short  = {
               "dist_a"      : "DistArea",
               "dist_p"      : "DistPop",
               "reock"       : "CircCircle",
               "inertia_a"   : "InertiaArea",
               "inertia_p"   : "InertiaPop",
               "polsby"      : "IPQ",
               "hull_a"      : "HullArea",
               "hull_p"      : "HullPop",
               "ehrenburg"   : "InscrCircle",
               "axis_ratio"  : "AxisRatio",
               "mean_radius" : "MeanRadius",
               "dyn_radius"  : "DynamicRadius",
               "harm_radius" : "HarmonicRadius",
               "rohrbach"    : "DistPerimeter",
               "exchange"    : "Exchange",
               "path_frac"   : "PathFraction",
               "power"       : "PowerDiagram",
               "split"       : "SplitLine"
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


def get_epsg (usps): return int(get_state_info(usps)["epsg"])
def get_fips (usps): return int(get_state_info(usps)["fips"])
def get_seats(usps): return int(get_state_info(usps)["seats"])

def get_state_info(usps):
   
   return pd.read_csv(stateinfo.format(usps.lower()) + ".csv").ix[0]


def cache_stateinfo(usps):

   filename = stateinfo.format(usps.lower())

   if os.path.exists(filename + ".csv"): return

   if not passwd:
     print("Failed -- no geo db authentication or cached files.")
     sys.exit(7)

   if blocks: import ps_block_query as query
   else:      import ps_query as query

   con = psycopg2.connect(database = "census", user = user, password = passwd,
                          host = "saxon.harris.uchicago.edu", port = 5432)

   info = pd.read_sql("SELECT epsg, fips, seats FROM states WHERE usps = upper('{}');".format(usps), con)
   info.to_csv(filename + ".csv", index = False)

   state_df = gpd.GeoDataFrame.from_postgis(query.states.format(usps), con,
                                            geom_col='state', crs = from_epsg(get_epsg(usps)))

   state_df[["id", "state"]].to_file(filename + ".shp")


def cache_shapefile(usps, blocks = False):

   tag = "_blocks" if blocks else ""
   filename = shapefile.format(usps.lower() + tag)

   if os.path.exists(filename): return

   if not passwd:
     print("Failed -- no geo db authentication or cached files.")
     sys.exit(7)

   if blocks: import ps_block_query as query
   else:      import ps_query as query

   state = get_state_info(usps)

   con = psycopg2.connect(database = "census", user = user, password = passwd,
                          host = "saxon.harris.uchicago.edu", port = 5432)

   geo_df = gpd.GeoDataFrame.from_postgis(query.shapes.format(state["fips"]), con,
                                          geom_col='geometry', crs = from_epsg(state["epsg"]))

   geo_df["id"] = geo_df.index
   geo_df[["id", "a", "pop", "x", "y", "split", "geometry"]].to_file(filename)



def cache_edge_file(usps, blocks = False):

   tag = "_blocks" if blocks else ""
   filename = edge_file.format(usps.lower() + tag)

   if os.path.exists(filename + ".csv"): return

   if not passwd:
     print("Failed -- no geo db authentication or cached files.")
     sys.exit(7)

   if blocks: import ps_block_query as query
   else:      import ps_query as query

   state = get_state_info(usps)

   con = psycopg2.connect(database = "census", user = user, password = passwd,
                          host = "saxon.harris.uchicago.edu", port = 5432)

   geo_df = gpd.GeoDataFrame.from_postgis(query.edges.format(state["fips"]), con,
                                          geom_col='lines', crs = from_epsg(state["epsg"]))

   geo_df[geo_df["rev"] == False][["eid", "lines"]].to_file(filename + ".shp")

   geo_df[["rn","seq","eid","rev","nodea","nodeb"]].to_csv(filename + ".csv", index = False)


def cache_node_file(usps, blocks = False):

   tag = "_blocks" if blocks else ""
   filename = node_file.format(usps.lower() + tag)

   if os.path.exists(filename + ".csv"): return

   if not passwd:
     print("Failed -- no geo db authentication or cached files.")
     sys.exit(7)

   if blocks: import ps_block_query as query
   else:      import ps_query as query

   state = get_state_info(usps)

   con = psycopg2.connect(database = "census", user = user, password = passwd,
                          host = "saxon.harris.uchicago.edu", port = 5432)

   ndf = pd.read_sql(query.nodes.format(state["fips"]), con)

   ndf[["nid", "x", "y", "nseq", "eid"]].to_csv(filename + ".csv", index = False)

   geometry = [Point(xy) for xy in zip(ndf.x, ndf.y)]
   geo_ndf = gpd.GeoDataFrame(ndf, crs = from_epsg(state["epsg"]), geometry=geometry)

   geo_ndf[ndf["nseq"] == 1][["nid", "geometry"]].to_file(filename + ".shp")


def cache_race_file(usps, tag = None):

  filename = race_file.format(usps.lower())
  if tag: filename = filename.replace(".", "_" + tag + ".")

  if os.path.exists(filename): return

  if not passwd:
    print("Failed -- no geo db authentication or cached files.")
    sys.exit(7)
    
  pd.read_sql("""SELECT
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
                 ORDER BY d.state, d.county, d.tract;
                 """.format(usps),
                 con = psycopg2.connect(database = "census", user = user, password = passwd,
                                        host = "saxon.harris.uchicago.edu", port = 5432),
                 index_col = "rn").to_csv(filename)



def plot_map(gdf, filename, crm, hlt = None, shading = "district", figsize = 10, label = "", ring = None, circ = None, point = None, scores = None, legend = False):

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

    if scores: dis["scores"] = pd.Series(scores)

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

      # if hlt: gdf[gdf["H"] == 1].plot(facecolor = "red", alpha = 0.1, linewidth = 0.05, ax = ax)

    elif "scores" in shading:

      col, alpha, trunc = "cool", 0.7, ""
    
      norm = Normalize(vmin = min(scores.values()), vmax = max([1, max(scores.values())]))
      # print(sum(scores.values()))
    
      cmap = plt.cm.ScalarMappable(norm=norm, cmap = col)
      
      ax = dis.plot(color = "white", edgecolor = "white", figsize = fs)
      for xi, row in dis.iterrows(): dis[dis.index == xi].plot(ax = ax, alpha = alpha, facecolor = cmap.to_rgba(row["scores"]))

      fig = ax.get_figure()
      cax = fig.add_axes([0.16, 0.13, 0.70, 0.015 * np.sqrt(xr/yr)])
      sm = plt.cm.ScalarMappable(cmap = col, norm=norm)
      sm._A = [] # gross
    
      cb = fig.colorbar(sm, cax = cax, alpha = alpha, # label = "Population / Target" + trunc, labelsize=12,
                        orientation='horizontal', drawedges = True)
      cb.locator = ticker.MaxNLocator(nbins=5)
      cb.formatter.set_useOffset(False)
      cb.set_label("Score", size=12)
      cb.ax.tick_params(labelsize=12)
      cb.dividers.set_visible(False)
      cb.update_ticks()

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


def save_json(filename, usps, method, uid, gdf, crm, metrics, tracts = True):

    of = open(filename, 'w')

    js = {"USPS" : usps.upper(), "FIPS" : get_fips(usps), 
          "Method" : method, "Seats" : get_seats(usps), "UID" : uid,
          "Score" : 0 if method is "PowerDiagram" else sum(v for v in metrics[method].values())/get_seats(usps)
         }

    gdf["C"] = pd.Series(crm)

    elections = []
    if usps and os.path.exists(vote_file.format(usps.lower())) and tracts:
      votes = pd.read_csv(vote_file.format(usps.lower()), index_col = "rn")
      votes = votes.filter(regex = '[DR][0-9]{2}', axis = 1)

      gdf = gdf.join(votes)
      vote_columns = list(votes.columns)
      elections = set([int(el[1:]) + (1900 if int(el[1:]) > 50 else 2000) for el in votes.columns])

    race = None
    if usps and os.path.exists(race_file.format(usps.lower())) and tracts:
      race = pd.read_csv(race_file.format(usps.lower()), index_col = "rn")
      race.rename(columns = {"pop" : "TotalPopulation"}, inplace = True)
      gdf = gdf.join(race)


    dis = gdf.dissolve("C", aggfunc='sum')
    dis.reset_index(inplace = True)

    target = dis["pop"].sum() / dis.shape[0]
    dis["frac"] = (dis["pop"] / target).map('{:.03f}'.format)

    dis["a"] *= 3.8610216e-7 # m2 to mi2
    dis["a"] = dis["a"].astype(int)
    dis.rename(columns = {"C" : "ID", "a" : "AreaSqMi", "pop" : "Population", "frac" : "TargetFraction"}, inplace = True)

    js["PopulationDeviation"] = (np.abs(dis["Population"] - target)/target).max()


    js["Elections"] = {}
    for elyr in elections:
      el = "{:02d}".format(elyr % 100)
      dis["Party " + el] = "R"
      dis.loc[dis["D" + el] > dis["R" + el], "Party " + el] = "D"

      dis["D" + el + " Share"] = dis["D" + el] / (dis["D" + el] + dis["R" + el])
      dis["R" + el + " Share"] = dis["R" + el] / (dis["D" + el] + dis["R" + el])

      js["Elections"][elyr] = {}
      js["Elections"][elyr]["DemSeats"] = int(sum(dis["D" + el + " Share"] > 0.5))
      js["Elections"][elyr]["RepSeats"] = int(sum(dis["R" + el + " Share"] > 0.5))
      js["Elections"][elyr]["DemFrac"]  = float(sum(dis["D" + el + " Share"])/js["Seats"])
      js["Elections"][elyr]["RepFrac"]  = float(sum(dis["R" + el + " Share"])/js["Seats"])

    for k, v in metrics.items(): dis[k] = pd.Series(v)

    js["Districts"] = []
    for ri, row in dis.iterrows():

      dist = {"ID" : int(ri), "Populations" : {}, "Elections" : {}, "Spatial" : {}}

      if race:
        dist["Populations"]["Total"] = int(row["Population"])
        dist["Populations"]["Black"] = int(row["black"])
        dist["Populations"]["Hispanic"] = int(row["hispanic"])
        dist["Populations"]["TargetFrac"] = float(row["Population"]/target)
        dist["Populations"]["BlackFrac"] = float(row["black"]/row["Population"])
        dist["Populations"]["HispanicFrac"] = float(row["hispanic"]/row["Population"])

        dist["Populations"]["VAP"] = int(row["total_vap"])
        dist["Populations"]["BlackVAP"] = int(row["black_vap"])
        dist["Populations"]["HispanicVAP"] = int(row["hispanic_vap"])
        dist["Populations"]["BlackVAPFrac"] = float(row["black_vap"]/row["total_vap"])
        dist["Populations"]["HispanicVAPFrac"] = float(row["hispanic_vap"]/row["total_vap"])

      for k in metrics: dist["Spatial"][k] = float(row[k])

      for elyr in elections: 
        el = "{:02d}".format(elyr % 100)
        dist["Elections"][elyr] = {}
        dist["Elections"][elyr]["Party"] = "R" if row["R" + el] > row["D" + el] else "D"
        dist["Elections"][elyr]["DemFrac"] = float(row["D" + el] / (row["D" + el] + row["R" + el]))
        dist["Elections"][elyr]["RepFrac"] = float(row["R" + el] / (row["D" + el] + row["R" + el]))
        dist["Elections"][elyr]["DemVotes"] = int(row["D" + el])
        dist["Elections"][elyr]["RepVotes"] = int(row["R" + el])

      dist["AreaSqMi"] = float(row.AreaSqMi)
      dist["Cells"] = list(int(k) for k, v in crm.items() if v == row.ID)

      js["Districts"].append(dist)


    with open(filename, 'w') as outfile: json.dump(js, outfile)



def save_geojson(gdf, filename, crm, usps = None, metrics = None, tracts = True):

    gdf["C"] = pd.Series(crm)

    elections, vote_columns = [], []
    if usps and os.path.exists(vote_file.format(usps.lower())) and tracts:
      votes = pd.read_csv(vote_file.format(usps.lower()), index_col = "rn")
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



