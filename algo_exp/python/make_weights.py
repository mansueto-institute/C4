#!/home/jsaxon/anaconda3/envs/gds/bin/python

import seaborn as sns
import pandas as pd
from pandas import IndexSlice as idx
import pysal as ps
import geopandas as gpd
import numpy as np
from sklearn import cluster
import pysal.contrib.clusterpy as cp

from fiona.crs import from_epsg
import psycopg2, psycopg2.extras

import warnings
with warnings.catch_warnings(): # SHUT UP!
    warnings.simplefilter("ignore"); 
    import matplotlib.pyplot as plt
    
from matplotlib.collections import LineCollection

import requests

from copy import deepcopy as dc

import os
import wget

from netrc import netrc
user, acct, passwd = netrc().authenticators("harris")

state_crs = { 6 : 3488, 15 : 2782, 16 : 3522, 17 : 3528,
             23 : 3554, 24 : 3559, 33 : 3613, 37 : 3631, 
             42 : 3364, 44 : 3653, 47 : 3661, 48 : 3665}


def state_seats(stfip):

  con = psycopg2.connect(database = "census", user = user, password = passwd,
                         host = "saxon.harris.uchicago.edu", port = 5432)
  
  sdf = pd.read_sql("SELECT seats FROM states WHERE fips = {};".format(stfip), con)

  return sdf.ix[0].seats


def get_tract_populations(stfip):

    pop_call = "http://api.census.gov/data/2014/acs5?get=NAME,B01001_001E&for=tract:*&in=state:%d"
    cjs = requests.get(pop_call % stfip).json() 
    cdf = pd.DataFrame(cjs[1:], columns = cjs[0])
    cdf["county"] = pd.to_numeric(cdf["county"])
    cdf["tract"] = pd.to_numeric(cdf["tract"])
    cdf["B01001_001E"] = pd.to_numeric(cdf["B01001_001E"])
    cdf.rename(columns = {"B01001_001E" : "pop"}, inplace = True)
    cdf.set_index(["county", "tract"], inplace = True)
    census_df = cdf[["pop"]]

    return census_df



def get_tract_shapes(stfip, FROM_FILE = False):

    if FROM_FILE:

        wget.download("http://www2.census.gov/geo/tiger/GENZ2015/shp/cb_2015_%d_cousub_500k.zip" % stfip, 
                      "census_shp/cb_2015_%d_cousub_500k.zip" % stfip)

        os.system("unzip census_shp/cb_2015_%d_cousub_500k.zip" % stfip)
    
        gdf = gpd.read_file("census_shp/cb_2015_%d_tract_500k.shp" % stfip).to_crs(epsg = state_crs[stfip])
    
        gdf.rename(columns = {"NAME" : "name"}, inplace = True)
        gdf["state"]  = pd.to_numeric(gdf["STATEFP"])
        gdf["county"] = pd.to_numeric(gdf["COUNTYFP"])
        gdf["tract"]  = pd.to_numeric(gdf["TRACTCE"])
    
        gdf = gdf[["county", "tract", "state", "name", "geometry"]]
        gdf["lon"] = [x.coords[0][0] for x in gdf.centroid.geometry]
        gdf["lat"] = [x.coords[0][1] for x in gdf.centroid.geometry]
        
    else:
        
        ## No password for local connections.
        con = psycopg2.connect(database = "census", user = user, password = passwd,
                               host = 'saxon.harris.uchicago.edu', port = 5432)

        sdf = pd.read_sql("SELECT epsg FROM states WHERE fips = {};".format(stfip), con)
        crs = sdf.ix[0].epsg
    
        sql= """SELECT 
                   state, county, tract,
                   ST_Transform(census_tracts_2015.geomsimp, epsg) as geometry, 
                   ST_X(ST_Transform(census_tracts_2015.centroid, epsg)) as lon,
                   ST_Y(ST_Transform(census_tracts_2015.centroid, epsg)) as lat
                FROM census_tracts_2015 
                JOIN states ON fips = census_tracts_2015.state 
                WHERE state = {};"""
    
        gdf = gpd.GeoDataFrame.from_postgis(sql.format(stfip), 
                                            con, geom_col='geometry', crs = from_epsg(crs))


    gdf.sort_values(by = ["county", "tract"], inplace = True)
    # gdf.sort_index(inplace = True)
    # gdf.reset_index(inplace = True)

    ## dissolve into counties, to save their centroids.
    county = gdf.dissolve("county")[["geometry"]]
    county["lon"] = [x.coords[0][0] for x in county.geometry.centroid]
    county["lat"] = [x.coords[0][1] for x in county.geometry.centroid]
    county = county[["lon", "lat"]]

    full_geo = pd.merge(gdf, county, left_on = "county", right_index = True, suffixes = ('', '_county'))

    return full_geo


def retrieve_and_merge_state(stfip, geo_from_file = False):

    census_df = get_tract_populations(stfip)
    geo_df    = get_tract_shapes(stfip, geo_from_file)

    geo_census_df = geo_df.join(census_df, on = ["county", "tract"])
    geo_census_df.sort_values(by = ["county", "tract"], inplace = True)

    geo_census_df.to_file("processed/st_%d_tract_pop.shp" % stfip)



def w2line_graph(w, ctrs):
    
    segments = []
    for i in w.id_order:
        origin = ctrs[i]
        for j in w.neighbors[i]:
            if j < i: continue
            segments.append([origin, ctrs[j]])

    lc_segs = LineCollection(segments)
    lc_segs.set_linewidth(0.5)
    return lc_segs



def get_mod_rook(stfip, VISUALIZE = True):

    '''
    Want to get rook weights, modified so that there
    is at least one path between any pair of tracts in a state.
    In other words: presuming that most islands aren't a 
    perfect multiple of of congressional seats, we need to be able
    to get back to the mainland.

    Naively, we'd do this by inspecting the Shimbel matrix,
    but that takes way too long to run.  So we just loop over 
    the rook weights, and merge all sets.  We then merge
    down the line, finding the point in the smallest subgraph
    that is closest to a foreigner.

    We then iterate on to the next-smallest subgraph until 
    only one remains.
    '''

    # This is the baseline that they'll all get merged onto.
    connected_w = ps.rook_from_shapefile("processed/st_%d_tract_pop.shp" % stfip)
    
    # I thought 'full' was the way to go, here,
    # but that doesn't yield a normal weight mat.
    full_w = ps.knnW_from_shapefile("processed/st_%d_tract_pop.shp" % stfip, 
                                    k = connected_w.n - 1)

    # stored with state-specific CRS;
    # for vincenty convert to plate caree (4326)
    # and then re-derive the centroids.
    gdf = gpd.read_file("processed/st_%d_tract_pop.shp" % stfip)
    centroids = np.stack([gdf.lon, gdf.lat], axis = 1) # use euclidean distance
    

    # identify subgraphs by merging rook weights.
    graphs = []
    for x in range(connected_w.n):
        
        # single star graph
        lone_star = set([x] + connected_w.neighbors[x])
        
        gidx = []
        for gi, g in enumerate(graphs):
            if any(n in g for n in lone_star):
                gidx.append(gi)
        
        if not gidx: # not yet found: new
            graphs.append(lone_star)
    
        else: # add to first
            graphs[gidx[0]] |= lone_star
    
            # loop backwards over any others,
            # merging and popping them off.
            for gi in gidx[:0:-1]:
                graphs[gidx[0]] |= graphs[gi]
                graphs.pop(gi)
    
                
    # Now the loop to merge graphs.
    # Looking for the closest to a larger graph.
    while len(graphs) > 1:
    
        graphs.sort(key = lambda x: len(x), reverse = True)
    
        for local_g in graphs[:0:-1]:
            
            wknn = ps.knnW_from_shapefile("processed/st_%d_tract_pop.shp" % stfip, len(local_g))
            
            # closest point in other subgraphs to this one.
            local, foreigner, min_dist2 = -1, -1, float('inf')
    
            # each point in the small graph -- 
            # usually just a few observations.
            for local_pt in local_g:
                
                # normally only one of these will be foreign.
                for foreign_pt in wknn.neighbors[local_pt]:
                    if foreign_pt in local_g: continue
                      
                    dist2 = pow(centroids[local_pt][0] - centroids[foreign_pt][0], 2) + \
                            pow(centroids[local_pt][1] - centroids[foreign_pt][1], 2)
                                    
                    # so this chooses the optimal point 
                    # in the local star graph.
                    if dist2 < min_dist2:
                        min_dist2 = dist2
                        foreigner = foreign_pt
                        local = local_pt 
                    
            if local < 0 or foreigner < 0: print "shit got fucked"
            bridge = ps.w_subset(full_w, [local, foreigner])
            
            connected_w = ps.w_union(connected_w, bridge)
            
            for foreign_g in graphs[:-1]:
                if foreigner in foreign_g:
                    foreign_g |= local_g
            
            graphs.pop()
    
    # Save the modified file for future use.
    fo = ps.open("processed/st_%d_tract_pop.gal" % stfip, 'w')
    fo.write(connected_w)
    fo.close()


    if VISUALIZE: 
    
        gdf = gpd.read_file("processed/st_%d_tract_pop.shp" % stfip)

        ax = gdf.plot(facecolor = 'white', edgecolor = "0.7", linewidth = 0.4)
        ax.add_collection(w2line_graph(connected_w, np.stack([gdf.lon, gdf.lat], axis = 1)))

        plt.axis('equal')
        ax.set_axis_off()
        plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)
        plt.margins(0,0)
        plt.gca().xaxis.set_major_locator(plt.NullLocator())
        plt.gca().yaxis.set_major_locator(plt.NullLocator())

        ax.figure.savefig("temp.pdf", bbox_inches='tight', pad_inches = 0)

        os.system("pdfcrop temp.pdf img/weights/%d.pdf --margins '20 20 20 20'" % stfip)


def run_clusterpy_districting(stfip, method = "maxp", VISUALIZE = True):

    gcdf = gpd.read_file("processed/st_%d_tract_pop.shp" % stfip)

    varis = ["lat", "lon", "lat_county", "lon_county", "pop"]

    state_pop = gcdf.sum()["pop"] 
    dist_min_pop  = 0.97 * state_pop / state_seats(stfip)

    layer = cp.Layer()
    cp.addArray2Layer(gcdf[varis].values, layer, names = varis)
    cp.addRook2Layer("processed/st_%d_tract_pop.gal" % stfip, layer)


    np.random.seed(1)
    if method == "maxp":
        layer.cluster('maxpTabu', varis, threshold = dist_min_pop,
                      wType = 'rook', tabuLength = 100, maxit = 200)

        # with open("regions", "a") as out:
        #     out.write("STATE == {}\n".format(stfip))
        #     out.write(str(layer.regions) + "\n\n")
        # print(layer.regions)


    if method == "arisel":
        layer.cluster('arisel', varis[:4], state_seats(stfip), wType='rook')
        gcdf['ARiSeL'] = layer.region2areas

    if method == "azp":
        layer.cluster('azp', varis[:4], state_seats(stfip), wType='rook')
        gcdf['AZP'] = layer.region2areas

    print(layer.region2areas)
    gcdf[method] = layer.region2areas

    if VISUALIZE: 

        ax = gcdf.dissolve("state").plot(alpha = 0.0, linewidth = 1.0)

        gcdf.plot(ax = ax, categorical = True, column = method, linewidth=0.25, 
                          cmap = 'Set1', alpha = 0.7, edgecolor = "white", figsize=(5, 5))

        gcdf.dissolve("county").plot(ax = ax, alpha = 0.0, linewidth = 0.5)

        plt.axis('equal')
        ax.set_axis_off()
        plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)
        plt.margins(0,0)
        plt.gca().xaxis.set_major_locator(plt.NullLocator())
        plt.gca().yaxis.set_major_locator(plt.NullLocator())

        ax.figure.savefig("temp.pdf", bbox_inches='tight', pad_inches = 0)

        os.system("pdfcrop temp.pdf img/districts/{}_{}.pdf --margins '20 20 20 20'".format(method, stfip))





for stfip in [53]: # [6, 15, 23, 42, 48]:
    retrieve_and_merge_state(stfip)
    get_mod_rook(stfip, VISUALIZE = True)
    run_clusterpy_districting(stfip, method = "maxp", VISUALIZE = True)
        

