#!/usr/bin/env python

import sys
import requests
import pandas as pd

import psycopg2
from netrc import netrc
user, acct, passwd = netrc().authenticators("harris")
user, acct, apikey = netrc().authenticators("census")

apibase = "https://api.census.gov/data/"
counties = apibase + "2015/acs5?key=%s&for=county:*&get=B01001_001E" % apikey
blocks   = apibase + "2015/acs5?key=%s&for=block+group:*&in=state:{}+county:{}&get=" % apikey

# Scale back your ambitions on the ACS 5-year: ethnic and age break-down not available at BG.
blocks  += "B01001_001E" # ,B01001B_001E,B01001I_001E,B05003_008E,B05003_019E,B05003B_008E,B05003B_019E,B05003I_008E,B05003I_019E"

counties = pd.DataFrame(data = requests.get(counties).json()[1:],
                        columns = ["p", "s", "c"])[["s", "c"]].astype(int)

with open ("block_pop2015.csv", "w") as out: 
    out.write("s,c,t,bg,population\n") # ,black,hispanic,vap,bvap,hvap\n")

for ni, row in counties.iterrows():
    print(row.s, row.c, blocks.format(row.s, row.c))
    jsre = requests.get(blocks.format(row.s, row.c)).json()
    df = pd.DataFrame(data = jsre[1:], columns = ["population", 
                                                  # "black", "hispanic", "vap_male", "vap_female", "bvap_male", "bvap_female", "hvap_male", "hvap_female",
                                                  "s", "c", "t", "bg"]).fillna(0).astype(int)
  
    # df["vap"]  = df["vap_male"]  + df["vap_female"]
    # df["bvap"] = df["bvap_male"] + df["bvap_female"]
    # df["hvap"] = df["hvap_male"] + df["hvap_female"]

    # df = df[["s", "c", "t", "bg", "population", "black", "hispanic", "vap", "bvap", "hvap"]].fillna(0).astype(int)
    df = df[["s", "c", "t", "bg", "population"]].fillna(0).astype(int)

    with open("block_pop2015.csv", "a") as out:
        df.to_csv(out, header = False, index = False)
  
sys.exit()

sql_con = psycopg2.connect(database = 'census', user = user, password = passwd,
                           host = "saxon.harris.uchicago.edu", port = 5432)

cur = sql_con.cursor()

## create a temp table with an appropriate schema
cur.execute("CREATE TEMP TABLE dpop (p int, s smallint, c smallint, t int, b int);")
print("CREATE TEMP TABLE dpop (p int, s smallint, c smallint, t int, b int);")

## copy the data into the table.
with open("block_pop2015.csv") as f: 
  cur.copy_from(f, "dpop", columns = ["p", "s", "c", "t", "b"], sep = ",")
  sql_con.commit()

print(pd.read_sql("select * from dpop;", sql_con))
cur.execute("""UPDATE census_bg_2015 SET pop = p FROM dpop 
               WHERE state = s AND county = c AND tract = t and bgroup = b;""")
sql_con.commit()

print(pd.read_sql("select * from census_bg_2015;", sql_con).head())

cur.execute("DROP TABLE dpop;")



