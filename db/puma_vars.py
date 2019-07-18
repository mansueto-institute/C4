#!/usr/bin/env python

import sys
import requests
import pandas as pd

import psycopg2
from netrc import netrc
user, acct, passwd = netrc().authenticators("harris")
user, acct, apikey = netrc().authenticators("census")

url = "https://api.census.gov/data/2015/acs5?get=B01001_001E,B23006_001E,B23006_009E,B23006_016E,B23006_023E&for=public%20use%20microdata%20area:*"

jsre = requests.get(url).json()
df = pd.DataFrame(data = jsre[1:], columns = ["pop", "Adults", "HS", "Coll", "BA", "state", "puma"]).fillna(0).astype(int)

df[["state", "puma", "pop", "Adults", "HS", "Coll", "BA"]].to_csv("puma_pop2015.csv", header = False, index = False)


sql_con = psycopg2.connect(database = 'census', user = user, password = passwd,
                           host = "saxon.harris.uchicago.edu", port = 5432)

cur = sql_con.cursor()

## create a temp table with an appropriate schema
cur.execute("CREATE TEMP TABLE dpop (s smallint, p int, N int, A int, H int, C int, B int);")

## copy the data into the table.
with open("puma_pop2015.csv") as f: 
  cur.copy_from(f, "dpop", columns = ["s", "p", "N", "A", "H", "C", "B"], sep = ",")
  sql_con.commit()

cur.execute("""UPDATE puma SET pop = N, adults = A, hs = H + C + B, ba = B
               FROM dpop 
               WHERE state = s AND puma = p;""")
sql_con.commit()

print(pd.read_sql("select * from puma;", sql_con).head())

cur.execute("DROP TABLE dpop;")



