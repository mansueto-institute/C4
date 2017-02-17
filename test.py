#!/usr/bin/env python 

import sys, os, re

import pycluscious as pycl
from cluster import *

usps, seed = "PA", 0
if len(sys.argv) > 1:
  usps = sys.argv[1]
if len(sys.argv) > 2:
  seed = int(sys.argv[2])


# DISTANCE_A, DISTANCE_P, INERTIA_A, INERTIA_P, HULL_A, POLSBY, REOCK, EHRENBURG, POLSBY_W, PATH_FRAC
methods = [["reock"     , pycl.ObjectiveMethod.REOCK], 
           ["dist_a"    , pycl.ObjectiveMethod.DISTANCE_A],
           ["dist_p"    , pycl.ObjectiveMethod.DISTANCE_P],
           ["inertia_a" , pycl.ObjectiveMethod.INERTIA_A],
           ["polsby"    , pycl.ObjectiveMethod.POLSBY],
           ["hull_a"    , pycl.ObjectiveMethod.HULL_A]]

method = methods[seed % len(methods)]

def get_seats(usps):

   con = psycopg2.connect(database = "census", user = user, password = passwd,
                          host = "saxon.harris.uchicago.edu", port = 5432)
   
   sdf = pd.read_sql("SELECT seats FROM states WHERE usps = upper('{}');".format(usps), con)

   return sdf.ix[0]["seats"]



def plot_map(gdf, filename, crm, bc, col = "qyb", figsize = 10, label = ""):

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
      ax = dis.plot("C", alpha = 0.8, categorical = True, 
                    legend = True, figsize = (figsize * np.sqrt(xr/yr), figsize * np.sqrt(yr/xr)))

      alpha, fc = 0.3, "grey"

    gdf[gdf["B"] == 1].plot(facecolor = fc, alpha = alpha, linewidth = 0.05, ax = ax)

    ax.set_xlim([bounds[0] - 0.1 * xr, bounds[2] + 0.2 * xr])
    ax.set_ylim([bounds[1] - 0.1 * yr, bounds[3] + 0.1 * yr])

    if label: ax.text(bounds[0] - 0.1*xr, bounds[1] - 0.1*yr, label, fontsize = 10)

    ax.set_axis_off()

    if not filename: return ax

    ax.figure.savefig(filename + ".pdf", bbox_inches='tight', pad_inches=0.05)
    plt.close('all')



# using pysal and shapely; very slightly modified from the contrib:
# https://github.com/pysal/pysal/blob/master/pysal/contrib/shared_perimeter_weights.py
def spw_from_shapefile(shapefile, norm = True):
    polygons = ps.open(shapefile, 'r').read()
    polygons = list(map(asShape,polygons))
    perimeters = [p.length if norm else 1. for p in polygons]
    Wsrc = ps.rook_from_shapefile(shapefile)
    new_weights = {}
    for i in Wsrc.neighbors:
        a = polygons[i]
        p = perimeters[i]
        new_weights[i] = [a.intersection(polygons[j]).length for j in Wsrc.neighbors[i]]

        # print(Wsrc.neighbors[i], new_weights[i])
    return ps.W(Wsrc.neighbors, new_weights)


shapefile = "shapes/" + usps.lower() + ".shp"

if not os.path.exists(shapefile):
  state = cluster(usps.upper(), 0, True)
  # state.loop_voronoi(4e3, 2.5e-3)
  state.gdf[["a", "pop", "x", "y", "C", "geometry"]].to_file(shapefile)

gdf = gpd.read_file(shapefile)

spw = spw_from_shapefile(shapefile, False)
gdf["prook_n"] = pd.Series(spw.neighbors)
gdf["prook_w"] = pd.Series(spw.weights)


seats = get_seats(usps)
if seats == 1: sys.exit()
u = pycl.universe(seats)


for xi, c in gdf.iterrows():
  geo = re.sub(r'^POLYGON(.*)', r'MULTIPOLYGON(\1)', str(c.geometry))
  u.add_cell(pycl.cell(xi, int(c["pop"]), c.x, c.y, c.a, 
                       {n:w for n, w in zip(c.prook_n, c.prook_w)}, geo))

u.contiguity_to_neighbors()
u.connect_graph()
u.trim_graph()

u.rand_districts(seed)
u.grow_kmeans(method[0] == "dist_p") # False is population growing

for i in range(10):

  if not i: continue
  elif i == 1: 
    u.oiterate(method[1], niter = 5000, tol = 0.02, alpha = 4, verbose = False)
  else: break

  crm = u.cell_region_map()
  bc = u.border_cells()
  
  for col in ["cat", "qyb"]:
    plot_map(gdf, "results/{}_s{:03d}_i{:03d}_{}_{}".\
                  format(usps, seed, i, method[0], col), crm, bc, col, 
                  label = "{}: {}, seed {}".format(usps, method[0], seed))

  print("completed iteration ::", i)


