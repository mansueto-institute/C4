#!/usr/bin/env python

import sys
import requests
import pandas as pd

import psycopg2
from netrc import netrc
user, acct, passwd = netrc().authenticators("harris")
user, acct, apikey = netrc().authenticators("census")

url = "https://api.census.gov/data/2015/acs5?get=B01001_001E&for=county:*"

jsre = requests.get(url).json()
df = pd.DataFrame(data = jsre[1:], columns = ["pop", "state", "county"]).fillna(0).astype(int)

df[["state", "county", "pop"]].to_csv("county_pop2015.csv", header = False, index = False)


sql_con = psycopg2.connect(database = 'census', user = user, password = passwd,
                           host = "saxon.harris.uchicago.edu", port = 5432)

cur = sql_con.cursor()

## create a temp table with an appropriate schema
cur.execute("CREATE TEMP TABLE dpop (s smallint, c int, N int);")

## copy the data into the table.
with open("county_pop2015.csv") as f: 
  cur.copy_from(f, "dpop", columns = ["s", "c", "N"], sep = ",")
  sql_con.commit()

print(pd.read_sql("SELECT * FROM dpop;", sql_con))
cur.execute("""UPDATE counties_2015 SET pop = N
               FROM dpop 
               WHERE state = s AND county = c;""")

sql_con.commit()

print(pd.read_sql("select * from counties_2016;", sql_con).head())

cur.execute("DROP TABLE dpop;")



