#!/usr/bin/env python

import sys
import requests
import pandas as pd

import psycopg2
from netrc import netrc
user, acct, passwd = netrc().authenticators("harris")
user, acct, apikey = netrc().authenticators("census")

apibase = "http://api.census.gov/data/"
api2000 = apibase + "2000/sf1?key=%s&get=P001001,P003004,P004002,P005001,P005004,P006002&for=tract:*&in=state:{}"    % apikey 

api2010 = apibase + "2010/sf1?key=%s&get=P0010001,P0030003,P0040003,P0100001,P0100004,P0110002&for=tract:*&in=state:{}" % apikey 

# 2010 endpoint seems to have switched to dec/2010, and the codes have also changed, dropping 0's.
# This now seems right -->> 
# api2010 = apibase + "2010/sf1?key=%s&get=P001001,P003003,P004003,P010001,P010004,P011002&for=tract:*&in=state:{}" % apikey 


api1990 = apibase + "1990/sf3?key=%s&for=tract:*&in=state:{}&get=" % apikey
vars90  = ['P0010001', 'P0080002', 'P0100001'] # normal
vars90 += sum([[s + str(n) for n in range(13, 32)] 
                  for s in ["P01300", "P014C0", "P014D0", "P015A0", "P015B0"]], [])


fips = [2, 1, 5, 4, 6, 8, 9, 11, 10, 12, 13, 15, 19, 16, 
        17, 18, 20, 21, 22, 25, 24, 23, 26, 27, 29, 28, 
        30, 37, 38, 31, 33, 34, 35, 32, 36, 39, 40, 41, 
	42, 44, 45, 46, 47, 48, 49, 51, 50, 53, 55, 54, 56]


if False: ## Do 1990.

    for state in fips:
    
        print(state, end = " ", flush = True)
        dfs = []
        for x in range(0, len(vars90), 50):
            jsre = requests.get(api1990.format(state) + ",".join(vars90[x:x+50])).json()
            dfs.append(pd.DataFrame(data = jsre[1:], columns = jsre[0]).fillna(0).astype(int))
            dfs[-1].set_index(["state", "county", "tract"], inplace = True)
    
        df = pd.concat(dfs, axis = 1).reset_index()
        df["vap"]  = df.filter(regex = "P013", axis = 1).sum(axis = 1)
        df["bvap"] = df.filter(regex = "P014", axis = 1).sum(axis = 1)
        df["hvap"] = df.filter(regex = "P015", axis = 1).sum(axis = 1)
    
        df.rename(columns = {"P0010001" : "population", "P0080002" : "black", "P0100001" : "hispanic"}, inplace = True)
    
        df = df[["state", "county", "tract", "population", "black", "hispanic", "vap", "bvap", "hvap"]]
    
    df.to_csv("tract_pop1990.csv", header = False, index = False)
    print()



for yr, addr in []: #[[2000, api2000], [2010, api2010]]:

  dfs = []
  print(yr, "::", end = " ")
  for st in fips:
    # print(st, end = " ", flush = True)
    print(addr.format(st))
    jsre = requests.get(addr.format(st)).json()
    dfs.append(pd.DataFrame(data = jsre[1:], columns = ["p", "b", "h", "v", "bv", "hv", "s", "c", "t"]))
  print()

  output = pd.concat(dfs).fillna(0)
  output = output[["s", "c", "t", "p", "b", "h", "v", "bv", "hv"]] # change order
  output = output.astype(int)
  output.to_csv("tract_pop{}.csv".format(yr), header = False, index = False)




for yr in [1990, 2000, 2010]:

  sql_con = psycopg2.connect(database = 'census', user = user, password = passwd,
                             host = "saxon.harris.uchicago.edu", port = 5432)

  cur = sql_con.cursor()

  ## create a temp table with an appropriate schema
  cur.execute("CREATE TEMP TABLE dpop (s smallint, c smallint, t int, p int, b int, h int, v int, bv int, hv int);")

  ## copy the data into the table.
  with open("tract_pop{}.csv".format(yr)) as f: 
    cur.copy_from(f, "dpop", columns = ["s", "c", "t", "p", "b", "h", "v", "bv", "hv"], sep = ",")
    sql_con.commit()

  print(pd.read_sql("select * from dpop;", sql_con))
  cur.execute("""UPDATE census_tracts_{} SET
                    pop = p, black = b, hispanic = h,
                    vap = v, bvap = bv, hvap = hv
                 FROM dpop 
                 WHERE state = s AND county = c AND tract = t""".format(yr))
  sql_con.commit()

  print(pd.read_sql("select * from census_tracts_{};".format(yr), sql_con).head())

  cur.execute("DROP TABLE dpop;")



