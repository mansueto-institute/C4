#!/usr/bin/env python 

import argparse
from sys import exit, stdout

from netrc import netrc
user, acct, passwd = netrc().authenticators("harris")

from fiona.crs import from_epsg
import psycopg2, psycopg2.extras
import geopandas as gpd
import numpy as np
import pandas as pd
import pysal as ps
import requests

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import ticker

from tqdm import tqdm

# Needed to keep the centroids in the boundary.
from shapely.geometry import MultiPolygon, Polygon, Point, LinearRing, asShape

import os, glob

def ens_dir(f, quiet = False):
  try: os.makedirs(f)
  except FileExistsError:
    if not quiet: print("File {} exists.  Cleaning it.".format(f))

    for f in glob.glob(f + "*"): os.remove(f)


def wavg(grp): 
  return grp._get_numeric_data().multiply(grp['a'], axis=0).sum()/grp['a'].sum()


class cluster:

    def __init__(self, usps, seed = 1, quiet = False):
      
        self.usps = usps.upper()
        self.seed = seed
        self.nloops = 0
      
        con = psycopg2.connect(database = "census", user = user, password = passwd,
                               host = "saxon.harris.uchicago.edu", port = 5432)
        
        sql = """SELECT usps, fips, epsg, seats,
                        ST_Scale(ST_Transform(geom, epsg), 0.001, 0.001, 0.001) AS geo 
                 FROM states WHERE usps = upper('{}');"""
        
        sdf = gpd.GeoDataFrame.from_postgis(sql.format(self.usps), con, geom_col='geo')

        sdf.crs = from_epsg(sdf.ix[0]["epsg"])

        state = sdf.ix[0]
        self.fips  = state.fips
        self.epsg  = state.epsg
        self.seats = state.seats
        self.area  = state.geo.area
        self.poly  = state.geo.simplify(5)

        if type(self.poly) is Polygon:
            self.rings = [LinearRing(self.poly.exterior.coords)]
        else: 
            self.rings = [LinearRing(p.exterior.coords) for p in self.poly]
            

        self.initialize_tract_frame()
        self.initialize_districts()

        self.voronoi_classify()

        self.folder = "power/{}/{:02d}/".format(self.usps, self.seed)
        ens_dir(self.folder, quiet)

        if not quiet: self.plot_map("init.pdf".format(usps, seed))
        
        self.print_freq = 500



    def initialize_districts(self):

        # A random sample of tracts, as seeds.
        np.random.seed(self.seed)
        cdf = self.gdf.ix[np.random.choice(self.ntracts, self.seats, replace=False)][["x", "y"]].reset_index(drop = True)
        
        # Initalize to just the largest radius...
        distances = ((cdf[["x", "y"]].as_matrix() - 
                      cdf[["x", "y"]].as_matrix()[:,np.newaxis,:])**2).sum(axis=2)
        cdf["r2"] = np.amin(np.ma.masked_equal(distances, 0), axis = 1).data

        self.cdf = cdf


    def initialize_tract_frame(self):
    
        con = psycopg2.connect(database = "census", user = user, password = passwd,
                               host = "saxon.harris.uchicago.edu", port = 5432)

        sql = """SELECT census_tracts_2015.county, census_tracts_2015.tract, 
                        ST_Scale(ST_Transform(geomsimp, states.epsg), 0.001, 0.001, 0.001) as geometry, 
                        ST_X(ST_Transform(census_tracts_2015.centroid, states.epsg))/1000 as x, 
                        ST_Y(ST_Transform(census_tracts_2015.centroid, states.epsg))/1000 as y,
                        census_tracts_2015.aland/1e6 as A, b01001_001e as pop
                 FROM census_tracts_2015
                 JOIN acssf5y2015 on
                   census_tracts_2015.state = acssf5y2015.state AND
                   census_tracts_2015.county = acssf5y2015.county AND
                   census_tracts_2015.tract = acssf5y2015.tract
                 JOIN states ON fips = census_tracts_2015.state 
                 WHERE census_tracts_2015.state = {}
                 ORDER BY county, tract;
              """
        
        geo_df = gpd.GeoDataFrame.from_postgis(sql.format(self.fips), 
                                            con, geom_col='geometry', crs = from_epsg(self.epsg))
    
        geo_df.to_file("proc/%d.shp" % self.fips)
        geo_df["wrook"] = pd.Series(ps.rook_from_shapefile("proc/%d.shp" % self.fips).neighbors)
        geo_df["points"] = [Point(xy) for xy in zip(geo_df.x, geo_df.y)]


        self.target = geo_df['pop'].sum() / self.seats
        self.ntracts = geo_df.shape[0]

        self.gdf = geo_df
    

    def scale(self, rate = 0.01, clip = 0.01):
    
        groups = self.gdf.groupby("C")
    
        dist_pop = groups["pop"].sum()
        pop = pd.Series(np.zeros(self.seats))
        pop.ix[dist_pop.index] = dist_pop
    
        dfactor = (self.target - pop) / self.target
        
        # Signed squared difference.
        clipped_factor = np.clip(np.abs(dfactor) * dfactor, -clip, +clip) * groups["a"].sum()
        clipped_factor.fillna(clip * self.area / self.seats, inplace = True)
        
        self.cdf["S"]   = clipped_factor
        self.cdf["r2"] += clipped_factor


    def voronoi_classify(self):
        
        # scipy.spatial.distance is equally fast...
        d2 = ((self.gdf[["x", "y"]].as_matrix() - 
               self.cdf[["x", "y"]].as_matrix()[:,np.newaxis,:])**2).sum(axis=2)
        
        self.gdf["C"] = np.argmin(d2 - self.cdf["r2"][:,np.newaxis], axis = 0)



    def plot_map(self, filename = None, annotate = True, figsize = 4):
    
        """
        Somewhat fancy function, for nicely monitoring
        the evolution of district centroids in cdf and
        tract points in gdf.
        """
    
        dis = self.gdf.dissolve("C", aggfunc='sum')
        dis["frac"] = dis["pop"] / self.target
    
        bounds = self.poly.bounds
        xr = bounds[2] - bounds[0]
        yr = bounds[3] - bounds[1]
        
        col, alpha, trunc = "coolwarm", 0.7, ""
    
        if dis["frac"].max() > 2: 
            norm = Normalize(vmin = 0, vmax = 2)
            trunc = " (Truncated)"
        elif dis["frac"].max() - 1 < 1e-5:
            norm = Normalize(vmin = 0.9999, vmax = 1.0001)
        else: # regardless, keep it centered
            larger = max(1 - dis["frac"].min(), dis["frac"].max() - 1)
            norm = Normalize(vmin = 1 - larger, vmax = 1 + larger) 
    
        cmap = plt.cm.ScalarMappable(norm=norm, cmap = col)
        
        ax = dis.plot(color = "white", edgecolor = "white", figsize = (figsize * np.sqrt(xr/yr), figsize * np.sqrt(yr/xr)))
        for xi, row in dis.iterrows():
            dis[dis.index == xi].plot(ax = ax, alpha = alpha, 
                                      facecolor = cmap.to_rgba(row["frac"]))
    
        if annotate: 
            ax.scatter(self.cdf.x, self.cdf.y, c = "k", s = 10, edgecolor = "w", marker = "o", zorder = 3)
    
            for ri, row in self.cdf.iterrows():
                ax.annotate(str(ri), (row.x, row.y), zorder = 4)
                ax.add_artist(plt.Circle((row.x, row.y), np.sqrt(max(row.r2, 1)), color="k", fill=False))    
    
        fig = ax.get_figure()
        cax = fig.add_axes([0.16, 0.13, 0.70, 0.015 * np.sqrt(xr/yr)])
        sm = plt.cm.ScalarMappable(cmap = col, norm=norm)
        sm._A = [] # gross
    
        cb = fig.colorbar(sm, cax = cax, alpha = alpha, label = "Population / Target" + trunc, orientation='horizontal')
        cb.locator = ticker.MaxNLocator(nbins=5)
        cb.formatter.set_useOffset(False)
        cb.update_ticks()
    
        ax.set_xlim([bounds[0] - 0.1 * xr, bounds[2] + 0.1 * xr])
        ax.set_ylim([bounds[1] - 0.1 * yr, bounds[3] + 0.1 * yr])
    
        ax.set_axis_off()
    
        if not filename: return ax
    
        ax.figure.savefig(self.folder + filename, bbox_inches='tight', pad_inches=0.05)
        plt.close('all')



    def closest_point_on_boundary(self, x, y):

        point = Point(x, y)
        
        mini = np.inf
        for ring in self.rings:
        
            d = ring.project(point)
            if d < mini:
                mini = d
                pt_in_poly = ring.interpolate(d)
        
        # return pt_in_poly
        return list(pt_in_poly.coords)[0]



    def go_where_its_hot(self, centers):
        
        # get the directions and distance squared to other cells.
        vectors = centers[["x", "y"]].as_matrix() - centers[["x", "y"]].as_matrix()[:,np.newaxis]
        distance = np.sqrt(np.square(vectors).sum(axis = 2))
        
        distance = np.ma.masked_equal(distance, 0)
        unit_vectors = vectors / distance[:,:,np.newaxis]
    
        # populations and strenths of the pulls.
        pop_diff = (centers["pop"].as_matrix() - centers["pop"].as_matrix()[:,np.newaxis]) / self.target
        pop_diff[pop_diff < 0] = 0
    
        print('\n:::: POP :::', (pop_diff/distance)[:,:,np.newaxis])
        move = (pop_diff / distance).T[:,:,np.newaxis]
    
        step = np.sum(vectors * move, axis = 0).data
        print("\n:::: STEP\n", step)
        return step




    def loop_voronoi(self, nloops, stop = 2.5e-3):

      if self.seats == 1: return # no point!

      nloops = int(nloops)

      ## for i in tqdm(range(self.nloops, self.nloops + nloops), ncols = 80, 
      ##               desc = "{} ({})".format(self.usps, self.seed)): # status bar!

      for i in range(self.nloops, self.nloops + nloops):
                    
          # calculate a scaling from the populations
          self.voronoi_classify()
    
          self.cdf["pop"] = self.gdf.groupby("C")["pop"].sum()

          if np.abs(1 - self.cdf["pop"]/self.target).max() < stop: break

          self.scale(rate = 0.01, clip = 0.01)
      
          # watch out for the nearest (district) neighbors
          # you get "zapped" for bullying smaller districts...
          distance2 = ((self.cdf[["x", "y"]].as_matrix() - # broadcast to seats × seats
                        self.cdf[["x", "y"]].as_matrix()[:,np.newaxis,:])**2).sum(axis=2)
          distance2 = np.ma.masked_equal(distance2, 0) # robust... 
          idx = np.nonzero(distance2 <= self.cdf.r2[np.newaxis,:])
          idx = np.array([np.append(idx[0], idx[1]), np.append(idx[1], idx[0])])
          
          # The zapping is away from the "adversary"
          # with each one stepping in proportion to its population.
          self.cdf["pop"].fillna(1, inplace = True)
          pop_ratio = self.cdf["pop"][:,np.newaxis] / \
                      (self.cdf["pop"][np.newaxis,:] + self.cdf["pop"][:,np.newaxis])
      
          D = pop_ratio[idx[0],idx[1]][:,np.newaxis] * \
              (self.cdf.loc[idx[1],["x", "y"]].values - self.cdf.loc[idx[0],["x", "y"]])
              
          D = D.groupby(D.index).sum()
          self.cdf.loc[D.index, ["x", "y"]] -= 0.05 * D
          
          # Make sure that it hasn't left the state!
          for idx, r in self.cdf.loc[D.index].iterrows():
              if not self.poly.contains(Point(r.x, r.y)):
                  self.cdf.loc[idx, ["x", "y"]] = self.closest_point_on_boundary(r.x, r.y)
          
          # at any rate, the power radius mustn't exceed
          # the distance to nearest neighbor's centroid.
          # recalculate distances post-zapping.
          distance2 = ((self.cdf[["x", "y"]].as_matrix() - # broadcast to seats × seats
                        self.cdf[["x", "y"]].as_matrix()[:,np.newaxis,:])**2).sum(axis=2)
          distance2 = np.ma.masked_equal(distance2, 0) # notes on masking method, below.
          closests = np.amin(distance2, axis = 1).data
          self.cdf["r2"] = np.minimum(self.cdf.r2, closests)
          
          # "regularize" -- can't go below negative avg area.
          # cdf.loc[(cdf["pop"] < 1.1 * self.target) & (cdf["r2"] < 1), "r2"] = 1
          # cdf["r2"] = np.maximum(cdf.r2, - state_area/seats)
      
          ############
          # Also move *slowly* towards the geographic center.
          self.voronoi_classify()
      
          # calculate lines towards the new points
          baryctr = self.gdf.groupby("C").apply(wavg)[["x", "y"]].dropna()
          old_pts = self.cdf.loc[baryctr.index,["x", "y"]]
      
          # how far to step.
          dist_pop = self.gdf.groupby("C")["pop"].sum().loc[baryctr.index]
          mv_rate = np.power(10., - 0.5 - dist_pop / self.target)
          mv_rate = mv_rate.values.reshape(1, mv_rate.size)
      
          self.cdf.loc[baryctr.index,["x", "y"]] = ((1 - mv_rate.T) * old_pts + \
                                                    mv_rate.T  * baryctr)
          
          ### annd.... move towards "hot" districts.
          # if i % 25 == 0: go_where_its_hot(self.cdf)
          # self.cdf[["x", "y"]] -= go_where_its_hot(self.cdf)
      
          # Monitoring and plotting.
          self.voronoi_classify()
    
          if i % self.print_freq == 0: 
              self.plot_map("{:05d}.pdf".format(i))

          
      self.plot_map("final.pdf")
      self.plot_map("final_clear.pdf", annotate = False)

      self.nloops += nloops


def main(usps = "IL", seed = 1, nloop = 4e3, stop = 2.5e-3):
    
    state = cluster(usps, seed)
    state.loop_voronoi(nloop, stop = stop)




if __name__ == "__main__":


  parser = argparse.ArgumentParser()
  parser.add_argument('--state', default = "ks", type=str.upper, help='state')
  parser.add_argument("--seed",  default = 1, type = int)
  parser.add_argument("--nloop", default = int(4e3), type = int)
  parser.add_argument("--stop",  default = 5e-3, type = float)
  args = parser.parse_args()

  args = parser.parse_args()
  main(args.state, args.seed, args.nloop, args.stop)


