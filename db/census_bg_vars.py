#!/usr/bin/env python

import sys
import requests
import pandas as pd

import psycopg2
from netrc import netrc
user, acct, passwd = netrc().authenticators("harris")
user, acct, apikey = netrc().authenticators("census")

fips = [2, 1, 5, 4, 6, 8, 9, 11, 10, 12, 13, 15, 19, 16,
        17, 18, 20, 21, 22, 25, 24, 23, 26, 27, 29, 28,
        30, 37, 38, 31, 33, 34, 35, 32, 36, 39, 40, 41,
        42, 44, 45, 46, 47, 48, 49, 51, 50, 53, 55, 54, 56]

fips = [1, 37, 42, 55]

apibase = "http://api.census.gov/data/"
counties90 = apibase + "1990/sf3?key=%s&get=P0010001&for=county:*" % apikey
bgroup90   = apibase + "1990/sf3?key=%s&for=block+group:*&in=state:{}+county:{}&get=" % apikey
vars90     = ['P0010001', 'P0080002', 'P0100001'] # normal
vars90    += sum([[s + str(n) for n in range(13, 32)] 
                  for s in ["P01300", "P014C0", "P014D0", "P015A0", "P015B0"]], [])

counties00 = apibase + "2000/sf1?key=%s&get=P001001&for=county:*" % apikey
bgroup00   = apibase + "2000/sf1?key=%s&get=P001001,P003004,P004002,P005001,P005004,P006002&for=block+group:*&in=state:{}+county:{}" % apikey
counties10 = apibase + "2010/sf1?key=%s&get=P0010001&for=county:*" % apikey 
bgroup10   = apibase + "2010/sf1?key=%s&get=P0010001,P0030003,P0040003,P0100001,P0100004,P0110002&for=block+group:*&in=state:{}+county:{}" % apikey



if False: ## Do 1990.

    ## 1990 is an obnoxious special case.
    # with open ("bgroup_pop1990.csv", "w") as out: 
    #     out.write("s,c,t,bg,population,black,hispanic,vap,bvap,hvap\n")

    print(counties90)
    counties = pd.DataFrame(data = requests.get(counties90).json()[1:],
                            columns = ["p", "s", "c"])[["s", "c"]].astype(int)

    for ni, row in counties.iterrows():
    
        if row.s < 19: continue
        print(row.s, row.c)
    
        dfs = []
        for x in range(0, len(vars90), 50):
            jsre = requests.get(bgroup90.format(row.s, row.c) + ",".join(vars90[x:x+50])).json()
            dfs.append(pd.DataFrame(data = jsre[1:], columns = jsre[0]).fillna(0).astype(int))
            dfs[-1].set_index(["state", "county", "tract", "block group"], inplace = True)
    
        df = pd.concat(dfs, axis = 1).reset_index()
        df["vap"]  = df.filter(regex = "P013", axis = 1).sum(axis = 1)
        df["bvap"] = df.filter(regex = "P014", axis = 1).sum(axis = 1)
        df["hvap"] = df.filter(regex = "P015", axis = 1).sum(axis = 1)
    
        df.rename(columns = {"P0010001" : "population", "P0080002" : "black", "P0100001" : "hispanic",
                             "block group" : "bg"}, inplace = True)
    
        df = df[["state", "county", "tract", "bg", "population", "black", "hispanic", "vap", "bvap", "hvap"]]
    
        with open("bgroup_pop1990.csv", "a") as out: 
            df.to_csv(out, header = False, index = False)


for y, co, bg in []: # [[2000, counties00, bgroup00], [2010, counties10, bgroup10]]:

    with open ("bgroup_pop{}.csv".format(y), "w") as out: 
        out.write("s,c,t,bg,population,black,hispanic,vap,bvap,hvap\n")

    print(co)
    counties = pd.DataFrame(data = requests.get(co).json()[1:],
                            columns = ["p", "s", "c"])[["s", "c"]].astype(int)

    for ni, row in counties.iterrows():
    
        print(row.s, row.c)
        jsre = requests.get(bg.format(row.s, row.c)).json()
        df = pd.DataFrame(data = jsre[1:], columns = ["population", "black", "hispanic", "vap", "bvap", "hvap", "s", "c", "t", "b"])
        df = df[["s", "c", "t", "b", "population", "black", "hispanic", "vap", "bvap", "hvap"]].fillna(0).astype(int)

        with open("bgroup_pop{}.csv".format(y), "a") as out: 
            df.to_csv(out, header = False, index = False)
        



sql_con = psycopg2.connect(database = 'census', user = user, password = passwd,
                           host = "saxon.harris.uchicago.edu", port = 5432)

cur = sql_con.cursor()

for year in [1990]: #, 2000, 2010]:

    ## create a temp table with an appropriate schema
    query = "CREATE TEMP TABLE dpop (s smallint, c smallint, t int, bg int, pop int, black int, hispanic int, vap int, bvap int, hvap int);"
    print(query)
    cur.execute(query)
    
    ## copy the data into the table.
    with open("bgroup_pop{}.csv".format(year)) as f: 
      cur.copy_from(f, "dpop", columns = ["s", "c", "t", "bg", "pop", "black", "hispanic", "vap", "bvap", "hvap"], sep = ",")
      sql_con.commit()

    print(pd.read_sql("select * from dpop limit 100;".format(year), sql_con))
    
    query = """UPDATE census_bg_{} AS bg SET
                 pop = dpop.pop, black = dpop.black, hispanic = dpop.hispanic,
                 vap = dpop.vap, bvap  = dpop.bvap,  hvap     = dpop.hvap     
               FROM dpop 
               WHERE 
                 bg.state  = dpop.s AND bg.county = dpop.c AND
                 (CASE WHEN (bg.tract % 100 = 0) THEN bg.tract/100 ELSE bg.tract END) = dpop.t AND
                 bg.bgroup = dpop.bg;""".format(year)
    print(query)
    cur.execute(query)
    sql_con.commit()
    
    print(pd.read_sql("select state, county, tract, bgroup, pop, black, hispanic, vap, bvap, hvap from census_bg_{} limit 100;".format(year), sql_con))
    
    cur.execute("DROP TABLE dpop;")



