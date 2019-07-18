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

# fips = [1, 37, 42, 55]

apibase = "http://api.census.gov/data/"
tracts00 = apibase + "2000/sf1?key=%s&get=P001001&for=tract:*&in=state:{}" % apikey
blocks00 = apibase + "2000/sf1?key=%s&get=P001001,P003004,P004002,P005001,P005004,P006002&for=block:*&in=state:{}+county:{}+tract:{}" % apikey
tracts10 = apibase + "2010/sf1?key=%s&get=P0010001&for=tract:*&in=state:{}" % apikey 
blocks10 = apibase + "2010/sf1?key=%s&get=P0010001,P0030003,P0040003,P0100001,P0100004,P0110002&for=block:*&in=state:{}+county:{}+tract:{}" % apikey

for y, tr, bl in []: # [[2000, tracts00, blocks00], [2010, tracts10, blocks10]]:

    with open ("block_pop{}.csv".format(y), "w") as out: 
        out.write("s,c,t,b,population,black,hispanic,vap,bvap,hvap\n")

    for s in sorted(fips):
    
        tracts = pd.DataFrame(data = requests.get(tr.format(s)).json()[1:],
                              columns = ["p", "s", "c", "t"])[["s", "c", "t"]].astype(int)
    
        for ni, row in tracts.iterrows():
            print(row.s, row.c, row.t)
            jsre = requests.get(bl.format(row.s, row.c, row.t)).json()
            df = pd.DataFrame(data = jsre[1:], columns = ["population", "black", "hispanic", "vap", "bvap", "hvap", "s", "c", "t", "b"])
            df = df[["s", "c", "t", "b", "population", "black", "hispanic", "vap", "bvap", "hvap"]].fillna(0).astype(int)

            with open("block_pop{}.csv".format(y), "a") as out: 
                df.to_csv(out, header = False, index = False)
        


sql_con = psycopg2.connect(database = 'census', user = user, password = passwd,
                           host = "saxon.harris.uchicago.edu", port = 5432)

cur = sql_con.cursor()

for year in [2000, 2010]:

    ## create a temp table with an appropriate schema
    query = "CREATE TEMP TABLE dpop (s smallint, c smallint, t int, b int, pop int, black int, hispanic int, vap int, bvap int, hvap int);"
    print(query)
    cur.execute(query)
    
    ## copy the data into the table.
    with open("block_pop{}.csv".format(year)) as f: 
      cur.copy_from(f, "dpop", columns = ["s", "c", "t", "b", "pop", "black", "hispanic", "vap", "bvap", "hvap"], sep = ",")
      sql_con.commit()

    print(pd.read_sql("select * from dpop limit 100;".format(year), sql_con))
    
    query = """UPDATE census_blocks_{} AS bl SET
                 pop      = dpop.pop      ,
                 black    = dpop.black    ,
                 hispanic = dpop.hispanic ,
                 vap      = dpop.vap      ,
                 bvap     = dpop.bvap     ,
                 hvap     = dpop.hvap     
               FROM dpop 
               WHERE 
                 bl.state  = dpop.s AND 
                 bl.county = dpop.c AND
                 bl.tract  = dpop.t AND
                 bl.block  = dpop.b;""".format(year)
    print(query)
    cur.execute(query)
    sql_con.commit()
    
    print(pd.read_sql("select state, county, tract, block, pop, black, hispanic, vap, bvap, hvap from census_blocks_{} limit 100;".format(year), sql_con))
    
    cur.execute("DROP TABLE dpop;")




