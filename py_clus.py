import math
import random
import pyproj

import psycopg2, psycopg2.extras

from shapely import wkt
from shapely.geometry import MultiPolygon as sh_mpoly

from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

from descartes import PolygonPatch

from netrc import netrc
user, acct, passwd = netrc().authenticators("harris")

import sys, os

def ens_dir(f):

  if not os.path.exists(f):
    os.makedirs(f)

class state_info():

  def __init__(self, census_num, abbrev, seats, epsg):
    self.num    = census_num
    self.abbrev = abbrev
    self.seats  = seats
    self.epsg   = epsg



# Use NSRS2007.
# Prefer full-state to partial, central to sides, South to North,
# unless there are significant features in the north (Oklahoma).
# http://spatialreference.org/ref/?search=nsrs2007

state_dict = {
  "Alabama"              : state_info( 1, "al",  7, 3465), # 3466
  "Alaska"               : state_info( 2, "ak",  1, 3467), 
  "Arizona"              : state_info( 4, "az",  9, 3478), # 3480, 3482
  "Arkansas"             : state_info( 5, "ar",  4, 3484), # 3486
  "California"           : state_info( 6, "ca", 53, 3488), 
  "Colorado"             : state_info( 8, "co",  7, 3501), # 3503, 3505
  "Connecticut"          : state_info( 9, "ct",  5, 3507), 
  "Delaware"             : state_info(10, "de",  1, 3509), 
  "District of Columbia" : state_info(11, "dc",  1, 3559), # (Maryland)
  "Florida"              : state_info(12, "fl", 27, 3513), 
  "Georgia"              : state_info(13, "ga", 14, 3518), # 3520 
  "Hawaii"               : state_info(15, "hi",  2, 2782), 
  "Idaho"                : state_info(16, "id",  2, 3522), 
  "Illinois"             : state_info(17, "il", 18, 2528), 
  "Indiana"              : state_info(18, "in",  9, 3532), 
  "Iowa"                 : state_info(19, "ia",  4, 3536), # 3538
  "Kansas"               : state_info(20, "ks",  4, 3540), # 3542
  "Kentucky"             : state_info(21, "ky",  6, 3546), 
  "Louisiana"            : state_info(22, "la",  6, 3552), # 3550
  "Maine"                : state_info(23, "me",  2, 3554), 
  "Maryland"             : state_info(24, "md",  8, 3559), 
  "Massachusetts"        : state_info(25, "ma",  9, 3585), 
  "Michigan"             : state_info(26, "mi", 14, 3587), # 89, 91
  "Minnesota"            : state_info(27, "mn",  8, 3594), # 95, 96
  "Mississippi"          : state_info(28, "ms",  4, 3597), # 99
  "Missouri"             : state_info(29, "mo",  8, 3601), # 02, 03
  "Montana"              : state_info(30, "mt",  1, 3604), 
  "Nebraska"             : state_info(31, "ne",  3, 3606), 
  "Nevada"               : state_info(32, "nv",  4, 3607),   
  "New Hampshire"        : state_info(33, "nh",  2, 3613), 
  "New Jersey"           : state_info(34, "nj", 12, 3615), 
  "New Mexico"           : state_info(35, "nm",  3, 3617), # 9, 21
  "New York"             : state_info(36, "ny", 27, 3623), # 25, 27, 29
  "North Carolina"       : state_info(37, "nc", 13, 3631),
  "North Dakota"         : state_info(38, "nd",  1, 3635), # 3633
  "Ohio"                 : state_info(39, "oh", 16, 3637), # 3638
  "Oklahoma"             : state_info(40, "ok",  5, 3639), # 3641
  "Oregon"               : state_info(41, "or",  5, 3643), # 45, 47
  "Pennsylvania"         : state_info(42, "pa", 18, 3364), 
  "Rhode Island"         : state_info(44, "ri",  2, 3653), 
  "South Carolina"       : state_info(45, "sc",  7, 3655), 
  "South Dakota"         : state_info(46, "sd",  1, 3659), # 3657
  "Tennessee"            : state_info(47, "tn",  9, 3661), 
  "Texas"                : state_info(48, "tx", 36, 3665), # 63, 66, 67, 71
  "Utah"                 : state_info(49, "ut",  4, 3675), 
  "Vermont"              : state_info(50, "vt",  1, 3684), 
  "Virginia"             : state_info(51, "va", 11, 3687), 
  "Washington"           : state_info(53, "wa", 10, 3689), # 91
  "West Virginia"        : state_info(54, "wv",  3, 3694), 
  "Wisconsin"            : state_info(55, "wi",  8, 3695), 
  "Wyoming"              : state_info(56, "wy",  1, 3703), # 02-05
  "Puerto Rico"          : state_info(72, "pr",  1, 2866)
}


class geop:

    epsg4326   = pyproj.Proj(init="epsg:4326")
    reprformat = "<tract (x,y) = (%.3f, %.3f)>"
    strformat  = "(y, x) = (%.3f, %.3f)"

    def __init__(self, lat, lon):

        self.lon = lon
        self.lat = lat 
        
        self.phi   = math.radians(lat)
        self.theta = math.radians(lon)

        self.x, self.y = geop.epsg4326(lon, lat)
        
        
    def __repr__(self):
        
        return geop.reprformat % (self.x, self.y)


    def __str__(self):
        
        return geop.strformat  % (self.x, self.y)

    def haversine(self, other): # Haversine formula
        
        dphi   = self.phi   - other.phi
        dtheta = self.theta - other.theta
        
        arg = (1 - math.cos(dphi))/2. + math.cos(self.phi) * math.cos(other.phi) * (1 - math.cos(dtheta))/2.   
        arg = math.sqrt(arg)
        
        dist = 2 * math.asin(arg)
        
        return dist

    def fasthav(self, other): # Small angle Haversine
                
        dphi   = self.phi   - other.phi
        dtheta = self.theta - other.theta
        
        arg = dphi * dphi + math.cos(self.phi) * math.cos(other.phi) * dtheta * dtheta   
        arg = math.sqrt(arg) / 2
        
        dist = 2 * math.asin(arg)
        
        return dist
    
    def projdist(self, other):
        ''' Projected Cartesian distance in the local geography.
            3 times faster than fast distance and 6 times faster
            than haversine or "official" haversine
        '''

        dx = self.x - other.x
        dy = self.y - other.y

        return math.sqrt(dx*dx + dy*dy)
     

class tract(geop):
    
    reprformat = "<tract w=%.3f  (x,y) = (%.5f, %.5f)>"
    strformat  = "pop=%d  (lat, lon) = (%.7f, %.7f)   (x, y) = (%.7f, %.7f)"

    def __init__(self, pop, lat, lon, mp = None, county = -1, tract = -1):
        
        geop.__init__(self, lat, lon) 

        self.mp  = mp
        self.pop = pop

        self.tract  = tract
        self.tracts = [tract]

        self.county = county
        
    def __repr__(self):
        
        return tract.reprformat % (self.pop, self.x, self.y)


    def __str__(self):
        
        return tract.strformat  % (self.pop, self.lat, self.lon, self.x, self.y)


    def merge(self, tract): 

        self.mp |= tract.mp

        self.pop += tract.pop
        tract.pop = 0

        self.tracts.append(tract.tract)

    
class district(geop):
    
    reprformat = "<district pop=%.3f  (x,y) = (%.3f, %.3f)>"
    strformat  = "pop=%.1f  w=%.4f (lat, lon) = (%.7f, %.7f)   (x, y) = (%.7f, %.7f)"

    ndistricts = 0
    
    def __init__(self, lat, lon, tracts):
        
        geop.__init__(self, lat, lon) 
        
        self.mp = sh_mpoly()

        self.pop = 0
        self.m = 1.0
        
        self.tracts = tracts # var from states.
        self.tract_idx = []
        self.tract_dist = []

        self.color = colors[district.ndistricts % len(colors)]
        district.ndistricts += 1
        
        
    def __repr__(self):
        
        return district.reprformat % (self.pop, self.x, self.y)

    
    def __str__(self):
        
        return district.strformat  % (self.pop, self.m, self.lat, self.lon, self.x, self.y)

    def add_tract(self, ti, dist = None):

        self.tract_idx .append(ti)
        if dist:
            self.tract_dist.append(dist)
        else:
            self.tract_dist.append(self.fasthav(self.tracts[ti]))

        self.update_centroid(plus = ti)
        self.update_population(plus = ti)
        self.update_multipolygon(plus = ti)
        
    def remove_tract(self, ti):
        
        idx = self.tract_idx.index(ti)
        
        self.tract_dist.pop(idx)
        self.tract_idx .pop(idx)

        self.update_centroid(minus = ti)
        self.update_population(minus = ti)
        self.update_multipolygon(minus = ti)

        
    def clear_tracts(self):
        self.tract_dist.clear()
        self.tract_idx.clear()
        
    def calculate_population(self):
        
        self.pop = 0
        for ti in self.tract_idx:
            self.pop += self.tracts[ti].pop
            

    def update_population(self, plus = -1, minus = -1):
        
        if plus  >= 0: self.pop += self.tracts[plus].pop
        if minus >= 0: self.pop -= self.tracts[minus].pop
            

    def calculate_centroid(self):

        if not self.pop: return (self.x, self.y) # don't update!
        
        xctr, yctr = 0, 0
        
        for ti, di in zip(self.tract_idx, self.tract_dist):
            
            tr = self.tracts[ti]
            
            xctr += tr.pop * tr.x
            yctr += tr.pop * tr.y
            
        self.x = xctr / self.pop
        self.y = yctr / self.pop
        
        self.lon, self.lat = self.epsg4326(self.x, self.y, inverse=True)
        
        return (self.x, self.y)

    def update_centroid(self, plus = -1, minus = -1):

        if plus >= 0:

            tract = self.tracts[plus]

            total_pop = self.pop + tract.pop

            self.x = (self.pop * self.x + tract.pop * tract.x) / total_pop
            self.y = (self.pop * self.y + tract.pop * tract.y) / total_pop

        if minus >= 0:

            tract = self.tracts[minus]

            total_pop = self.pop - tract.pop

            self.x = (self.pop * self.x - tract.pop * tract.x) / total_pop
            self.y = (self.pop * self.y - tract.pop * tract.y) / total_pop


        self.lon, self.lat = self.epsg4326(self.x, self.y, inverse=True)
        
        return (self.x, self.y)


    def calculate_multipolygon(self):
        
        self.mp.empty()
        
        for ti in self.tract_idx:
            self.mp |= self.tracts[ti].mp


    def update_multipolygon(self, plus = -1, minus = -1):
        
        if plus >= 0:
            self.mp |= self.tracts[plus].mp

        if minus >= 0:
            self.mp -= self.tracts[minus].mp
    
                
    def update_weight(self, target, rate = 1.0):
        
        frac = self.pop / target      
        
        if   frac < 0.70: self.m *= 1 + 0.15  * rate
        elif frac < 0.80: self.m *= 1 + 0.10  * rate
        elif frac < 0.90: self.m *= 1 + 0.05  * rate
        elif frac < 0.95: self.m *= 1 + 0.04  * rate
        elif frac < 0.96: self.m *= 1 + 0.03  * rate
        elif frac < 0.97: self.m *= 1 + 0.02  * rate
        elif frac < 0.98: self.m *= 1 + 0.01  * rate
        elif frac < 0.99: self.m *= 1 + 0.005 * rate
            
        if   frac > 1.30: self.m *= 1 - 0.15  * rate
        elif frac > 1.20: self.m *= 1 - 0.10  * rate
        elif frac > 1.10: self.m *= 1 - 0.05  * rate
        elif frac > 1.05: self.m *= 1 - 0.04  * rate
        elif frac > 1.04: self.m *= 1 - 0.03  * rate
        elif frac > 1.03: self.m *= 1 - 0.02  * rate
        elif frac > 1.02: self.m *= 1 - 0.01  * rate
        elif frac > 1.01: self.m *= 1 - 0.005 * rate

            
    def reset_weight(self, val):
        
        self.m = val
            

    def schwartzberg(self, plus = sh_mpoly(), minus = sh_mpoly()):

        mp = self.mp | plus - minus

        return mp.length / (2 * math.sqrt(math.pi * mp.area))


    def hull_ratio(self, plus = sh_mpoly(), minus = sh_mpoly()):

        mp = (self.mp | plus) - minus

        area = mp.area

        if area == 0: return float('inf')

        return mp.convex_hull.area / area


colors = ['red', 'blue', 'limegreen', 'yellow', 'mediumorchid', 
          'darkorange', 'salmon', 'aquamarine', 'deeppink', 'lightsteelblue',
          'darkgray', 'white', 'cornflowerblue', 'sage', 'pink', 'aqua', 'gold', 'teal']

class state():
    
    ## POSSIBLE WAYS FOR INITIALIZING K-MEANS SEEDS.
    CURRENT_CENTROIDS, FORGY, RANDOM_PARTITION = 0, 1, 2

    ## Methods for the trading algorithm
    KMEANS, SCHWARTZBERG, CONVEX_HULL = 0, 1, 2
    
    state_query = "SELECT ST_AsText(ST_Transform(geom, 4269)) boundary FROM states WHERE fips = %d;"

    cd_query = "SELECT ST_x(ST_Transform(centroid, 4269)) x, ST_y(ST_Transform(centroid, 4269)) y FROM cd114 WHERE fips = %d;"
    
    tract_query = '''
                   SELECT
                     "b01001_001e" pop,
                     ST_X(ST_Transform(centroid, {})) x,
                     ST_Y(ST_Transform(centroid, {})) y,
                     ST_AsText(ST_Transform(geomsimp, {})) boundary,
                     census_tracts_2015.county, 
                     census_tracts_2015.tract
                   FROM census_tracts_2015
                   JOIN acssf5y2015 ON (
                     census_tracts_2015.state  = acssf5y2015.state AND
                     census_tracts_2015.county = acssf5y2015.county AND
                     census_tracts_2015.tract  = acssf5y2015.tract
                   )
                   WHERE 
                     census_tracts_2015.state = {} AND 
                     census_tracts_2015.aland > 0
                   ;
                   '''
    
    def __init__(self, name, district_init = 0, seed = 0):

        if name not in state_dict:
            print("%s :: NOT A STATE NAME" % name)
            return
        
        self.name = name
        self.id = state_dict[name].num
        self.abbrev = state_dict[name].abbrev

        epsg = 4269 # state_dict[name].epsg
        self.tract_query = state.tract_query.format(epsg, epsg, epsg, state_dict[name].num)

        self.tag = self.abbrev

        self.init_method = district_init

        self.nfixed_contiguity = 0
        
        self.mp = None
        self.db_load_state_border()
        
        self.ntracts = 0
        self.tracts = []
        self.db_load_census_tracts()
        self.merge_contained_tracts()
        
        self.ndistricts = 0
        self.districts  = []
        if district_init is state.CURRENT_CENTROIDS:
            
            self.tag += "_current"

            self.db_load_district_centroids()
            
        elif district_init is state.FORGY:

            self.tag += "_forgy%d" % seed
            
            random.seed(seed)
            
            self.db_load_ndistricts()
            rand_dist = []
            while len(rand_dist) < self.ndistricts:
                
                t = int(random.random() * self.ntracts)
                if t not in rand_dist:
                    rand_dist.append(t)
                    
                    lat, lon = self.tracts[t].lat, self.tracts[t].lon
                    self.districts.append(district(lat, lon, self.tracts))
            
        else: 
            print("DISTRICT INIT METHOD NOT IMPLEMENTED!")
            
        self.pop = 0
        self.target = 0
        self.calculate_populations()

        self.assign_tracts_k_means()

        
    def db_load_state_border(self):
        
        conn = psycopg2.connect(database = "census", user = user, password = passwd,
                                host = "saxon.harris.uchicago.edu", port = 5432,
                                cursor_factory=psycopg2.extras.DictCursor)
        
        cur  = conn.cursor()
        cur.execute(state.state_query % self.id)
        
        self.mp = wkt.loads(cur.fetchall()[0]['boundary'])
        
        conn.close()
        

    def db_load_census_tracts(self):
        
        conn = psycopg2.connect(database = "census", user = user, password = passwd,
                                host = "saxon.harris.uchicago.edu", port = 5432,
                                cursor_factory=psycopg2.extras.DictCursor)
        
        cur  = conn.cursor()
        cur.execute(self.tract_query)
        
        results = cur.fetchall()
        
        for r in results:
    
            mp = wkt.loads(r['boundary'])
            
            self.tracts.append(tract(r['pop'], r['y'], r['x'], mp,
                                     r['county'], r['tract']))
        
        conn.close()
        
        self.ntracts = len(self.tracts)

    
    def db_load_ndistricts(self):
        
        self.districts.clear()

        conn = psycopg2.connect(database = "census", user = user, password = passwd,
                                host = "saxon.harris.uchicago.edu", port = 5432,
                                cursor_factory=psycopg2.extras.DictCursor)
        
        cur  = conn.cursor()
        cur.execute("SELECT COUNT(cd) N FROM cd114 WHERE fips = %d;" % self.id)
        
        self.ndistricts = cur.fetchone()[0]

        conn.close()
        
        
    def db_load_district_centroids(self):
        
        self.districts.clear()

        conn = psycopg2.connect(database = "census", user = user, password = passwd,
                                host = "saxon.harris.uchicago.edu", port = 5432,
                                cursor_factory=psycopg2.extras.DictCursor)
        
        cur  = conn.cursor()
        cur.execute(state.cd_query % self.id)
        
        cd_loc = cur.fetchall()
        
        for r in cd_loc:
    
            ctr = wkt.loads(r['ctr'])
            self.districts.append(district(r['x'], r['y'], self.tracts))

        conn.close()
        
        self.ndistricts = len(self.districts)
        
        
    def calculate_populations(self):
        
        self.pop = sum([t.pop for t in self.tracts])
        
        if self.ndistricts:
            self.target = 1. * self.pop / self.ndistricts
        else: 
            print("WARNING: NO DISTRICTS LOADED, SETTING TARGET TO 0.")
            self.target = 0


    def merge_contained_tracts(self):

        to_remove = []

        for tr1 in self.tracts:

            poly_list = tr1.mp
            if tr1.mp.geom_type is "Polygon":
                poly_list = [tr1.mp]

            for p in poly_list:

                for int_ring in p.interiors:
                    
                    for tr2 in self.tracts:

                        if tr1 is tr2: continue

                        if int_ring.touches(tr2.mp):
                            tr1.merge(tr2)
                            to_remove.append(tr2)

        to_remove.sort(key = lambda x: x.pop)

        for rm in set(to_remove):

            self.tracts.remove(rm)
            self.ntracts -= 1


        
    def assign_tracts_k_means(self, additive = True):
        
        for d in self.districts:
            d.clear_tracts()
            
        # find the minimum weighted distance
        for ti, t in enumerate(self.tracts):

            idxt = -1
            mini = float('inf')
            for cdi, cd in enumerate(self.districts):

                dist = cd.fasthav(t)

                if additive:
                  dist = dist ** 2 - cd.m
                  if mini > dist:
                    mini = dist
                    idxt = cdi
  
                else: 
                  dist = dist / cd.m
                  if mini > dist:
                    mini = dist
                    idxt = cdi

            self.districts[idxt].add_tract(ti, mini)

        for d in self.districts:
            d.calculate_population()
            d.calculate_centroid()
            d.calculate_multipolygon()

        # self.ensure_contiguity()


    def ensure_contiguity(self):

        split_districts = sum([int(d.mp.geom_type is "MultiPolygon") for d in self.districts])

        if split_districts: self.nfixed_contiguity += 1

        while split_districts:

            # print("%d non-contiguous" % sum([d.mp.geom_type is "MultiPolygon" for d in self.districts]))

            split_tracts = 0

            for di, d in enumerate(self.districts):

                if d.mp.geom_type is "Polygon": continue
                # print("found the multipolygon at", di)

                # kill the splinter with the lowest population.
                splinter_poly = min([poly for poly in d.mp], key = lambda x: self.population_in_polygon(x))
                

                # print("constructed splinter_poly")

                for ti in d.tract_idx:
                    
                    ti_mp = self.tracts[ti].mp

                    if not ti_mp.touches(splinter_poly.boundary): continue
                    if not ti_mp.intersection(splinter_poly).geom_type is "Polygon": continue

                    # edge case: split tract
                    if ti_mp.geom_type == "MultiPolygon" and \
                       not ti_mp.within(splinter_poly):
                        # print("ti %d, county %d, tract %d :: not reassigned -- special case of multipolygon tract." % \
                        #       (ti, self.tracts[ti].county, self.tracts[ti].tract))
                        split_tracts += 1
                        break
                    
                    # print("ti %d, county %d, tract %d, touches %d, intersection %s, within %d" % \
                    #       (ti, self.tracts[ti].county, self.tracts[ti].tract, 
                    #        ti_mp.touches(splinter_poly.boundary), 
                    #        ti_mp.intersection(splinter_poly.boundary).geom_type, 
                    #        ti_mp.within(splinter_poly)))

                    if not ti_mp.within(splinter_poly): continue

                    reassigned = False
                    for dio, do in enumerate(self.districts):

                        if do is d: continue
                        if not ti_mp.touches(do.mp): continue
                        if "Point" in ti_mp.intersection(do.mp).geom_type: continue
                       
                        # print("reassign", di, dio)
                        # print("ensure :: reassigning tr%d to dist%d" % (ti, dio))
                        self.reassign_tract(ti, do)
                        reassigned = True
                        break

                    if reassigned: break

            split_districts = sum([int(d.mp.geom_type is "MultiPolygon") for d in self.districts])
            split_districts -= split_tracts

                
    def reassign_tract(self, tidx, district): 

        for d in self.districts:
            if tidx in d.tract_idx:
                d.remove_tract(tidx)
                break

        district.add_tract(tidx)

        
    def reset_weights(self, val):
        
        for cd in self.districts:
            cd.reset_weight(val = val)
            
    def weighted_kmeans(self, niter = 100, tolerance = 0.05, print_multiple = 0, additive = True):
        
        folder = self.tag + "_weighted_kmeans"

        os.system("rm -rf " + folder)
        ens_dir(folder)

        if additive: self.reset_weights(1e-3)
        else: self.reset_weights(1)

        self.assign_tracts_k_means(additive)
        self.print_iteration(folder + "/init", "Initial")
            
        for i in range(niter):

            self.assign_tracts_k_means(additive)
        
            for cd in self.districts:
                cd.update_weight(self.target, rate = 0.10 * math.pow(0.99, i))  
            
            if all([(math.fabs(1 - cd.pop/self.target) < tolerance) for cd in self.districts]):
                print("\nReached required tolerance of %.3f on iteration %d, terminating." % (tolerance, i))
                self.print_iteration(folder + "/final", "Final: Iteration %d" % i, show = True)
                return

            self.districts.sort(key = lambda x: x.pop)
    
            if not i:
            
                self.print_iteration(folder + "/init", "Initial")
            
            elif print_multiple and not (i % print_multiple):
                print("%03d - Pop" % i, " ".join(["%.3f" % (cd.pop/self.target) for cd in self.districts]))
                self.print_iteration(folder + "/%03d" % i, "Iteration %d" % i)
                

    def population_in_polygon(self, poly):

        pop_in_polygon = 0
        for tr in self.tracts:
            if tr.mp.within(poly):
                pop_in_polygon += tr.pop

        return pop_in_polygon


    def grow_district(self, method, di, tolerance):

        dist = self.districts[di]

        # find the adjacent tracts from larger districts
        adjacent = []
        for dl in self.districts[di:]:

            # if the districts aren't adjacent, continue.
            if not dl.mp.touches(dist.mp): continue
            
            for ti in dl.tract_idx:
                
                ti_mp = self.tracts[ti].mp

                if not ti_mp.touches(dist.mp): continue
                if "Point" in ti_mp.intersection(dist.mp).geom_type: continue

                mini_param = float('inf')
                if method is state.KMEANS:

                    mini_param = dist.fasthav(self.tracts[ti])

                elif method is state.SCHWARTZBERG:
                             
                    nom = dist.schwartzberg() + dl.schwartzberg() 
                    mod = dist.schwartzberg(plus = self.tracts[ti].mp) \
                          + dl.schwartzberg(minus = self.tracts[ti].mp)

                    mini_param = mod / nom


                elif method is state.CONVEX_HULL:

                    nom = dist.hull_ratio() + dl.hull_ratio() 
                    mod = dist.hull_ratio(plus = self.tracts[ti].mp) \
                          + dl.hull_ratio(minus = self.tracts[ti].mp)

                    mini_param = mod / nom

                adjacent.append([ti, mini_param])

        adjacent.sort(key = lambda x: x[1])

        # assign them in order of the impact.
        # potentially multiple tracts, if we're 
        # reasonably far from the target.
        for adj, mp in adjacent:
            self.reassign_tract(adj, dist)
            if dist.pop / self.target > 0.85 or \
               dist.pop / self.target >  1 - 4 * tolerance:
                break

            if self.nfixed_contiguity > 15: break


    def trimming_district(self, method, di, tolerance = 0.0):


        dist = self.districts[di]

        # find the adjacent tracts from smaller districts
        border = []
        for ti in dist.tract_idx:

           for ds in self.districts:

                if ds is dist: continue

                if ds.pop > dist.pop and ds.pop > self.target * (1 + tolerance): continue

                ti_mp = self.tracts[ti].mp

                if not ti_mp.touches(ds.mp): continue
                if "Point" in ti_mp.intersection(ds.mp).geom_type: continue

                mini_param = float('inf')
                if method is state.KMEANS:

                    mini_param = ds.fasthav(self.tracts[ti])  

                elif method is state.SCHWARTZBERG:

                    nom = dist.schwartzberg() + ds.schwartzberg() 
                    mod = dist.schwartzberg(minus = self.tracts[ti].mp) \
                          + ds.schwartzberg(plus = self.tracts[ti].mp)

                    mini_param = mod / nom

                elif method is state.CONVEX_HULL:

                    nom = dist.hull_ratio() + ds.hull_ratio() 
                    mod = dist.hull_ratio(minus = self.tracts[ti].mp) \
                          + ds.hull_ratio(plus = self.tracts[ti].mp)

                    mini_param = mod / nom

                border.append([ti, ds, mini_param])


        border.sort(key = lambda x: x[2])

        # assign them in order of the impact.
        # potentially multiple tracts, if we're 
        # reasonably far from the target.
        for tract, dswap, mp in border:
            self.reassign_tract(tract, dswap)
            if dist.pop / self.target < 1.15 or \
               dist.pop / self.target < 1 + 4 * tolerance:
                break

            if self.nfixed_contiguity > 15: break


    def trading_opt(self, method = 0, niter = 100, tolerance = 0.015, print_multiple = 0, stop_at_tolerance = True):
        
        folder = self.tag + "_trading_opt_%.03f" % tolerance
        if method == state.SCHWARTZBERG: folder += "_schw"
        if method == state.CONVEX_HULL:  folder += "_cvxh"

        os.system("rm -rf " + folder)
        ens_dir(folder)
        
        self.print_iteration(folder + "/init", "Initial")
        
        self.districts.sort(key = lambda x: x.pop)

        for i in range(niter):

            if stop_at_tolerance and \
               self.districts[0].pop > (1 - tolerance) * self.target and \
               self.districts[-1].pop < (1 + tolerance) * self.target:
                break

            # deal first with the smaller districts.
            for di, d in enumerate(self.districts):
                
                if d.pop / self.target > 1 + tolerance: break
                self.grow_district(method, di, tolerance)


            # self.ensure_contiguity()

            # then take on the larger districts
            for di, d in enumerate(self.districts):

                if d.pop / self.target < 1 - tolerance: continue
                self.trimming_district(method, di, tolerance)

            # self.ensure_contiguity()


            self.districts.sort(key = lambda x: x.pop)
            if print_multiple and not (i % print_multiple):
                print("%03d - Pop" % i, " ".join(["%.03f" % (1.0 * cd.pop / self.target) for cd in self.districts]))
                self.print_iteration(folder + "/%03d" % i, "Iteration %d" % i, show = False)


        print("End - Pop", " ".join(["%.03f" % (1.0 * cd.pop/self.target) for cd in self.districts]))
        self.print_iteration(folder + "/final", "Final", show = True)
        
        os.system("convert -loop 0 -delay 300 %s/init.png -delay 50 %s/???.png  -delay 1000 %s/final.png %s/animation.gif" % (folder, folder, folder, folder))

        
    def print_iteration(self, name = "", title = "", size = 10, show = False, img_format = "png"):

        if not title: title = self.name

        fig = plt.figure(figsize=(size, size))
        
        ll_lon, ll_lat, ur_lon, ur_lat = self.mp.bounds
        ll_lon -= 0.10 * (ur_lon - ll_lon)
        ur_lon += 0.10 * (ur_lon - ll_lon)
        ll_lat -= 0.10 * (ur_lat - ll_lat)
        ur_lat += 0.10 * (ur_lat - ll_lat)
    
        m = Basemap(llcrnrlon=ll_lon,llcrnrlat=ll_lat,urcrnrlon=ur_lon,urcrnrlat=ur_lat, 
                    resolution = 'c', epsg=4326)
        m.drawmapboundary(fill_color = 'white')

        for p in self.mp:
            poly = PolygonPatch(p, fc='black', ec = 'k', linewidth = 1, zorder=0)
            plt.gca().add_patch(poly) 

        for cd in self.districts:

            # print(cd)

            mp_district = sh_mpoly()

            for ti in cd.tract_idx:
                
                mp_district |= self.tracts[ti].mp
                
            if mp_district.geom_type is "Polygon":
                
                poly = PolygonPatch(mp_district, fc = cd.color, ec = 'k', linewidth = 2, zorder=0)
                plt.gca().add_patch(poly)
                
            else:
                
                for p_district in mp_district:
                    poly = PolygonPatch(p_district, fc = cd.color, ec = 'k', linewidth = 2, zorder=0)
                    plt.gca().add_patch(poly)    
    
        plt.title(title, fontsize = 20)
    
        if show: plt.show()
        if name: fig.savefig(name + "." + img_format, dpi=fig.dpi, bbox_inches='tight')
        plt.close(fig)

