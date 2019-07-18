#!/usr/bin/env python

# CREATE TABLE acs2015 AS SELECT state, county, tract FROM census_tracts_2015;

import pandas as pd
import numpy as np

import psycopg2

import argparse

from netrc import netrc
user, acct, passwd = netrc().authenticators("harris")
user, acct, apikey = netrc().authenticators("census")

apikey = "a4b2eab7c7050050923fffa485fb81e22be63e68"

import requests, json
import sys

def main(table = "acssf5y2015"):

    year = int(table.lower().split("y")[1])

    cols, cvars = [], []
    readable = {}
    description = {}
    print(table.lower() + ".csv")
    for line in open(table.lower() + ".csv", "r"):

        spline = [x.strip() for x in line.split(",")]

        var = spline[0]
        col = "" if len(spline) < 2 else spline[1].strip()
        
        cvars.append(var)
        cols.append(col)

        readable[var]    = col
        description[col] = "" if len(spline) < 3 else spline[2].strip()

##      dfs = []
##      for st in [1, 2, 4, 5, 6, 8, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 53, 54, 55, 56]:
##  
##          st = int(st)
##          print(st, ",".join(cvars))
##          if st > 57: continue
##  
##          url = "https://api.census.gov/data/{}/acs/acs5/profile?for=tract:*&in=state:{:02d}&key={}&get=".format(year, st, apikey)
##          url += ",".join(cvars)
##          print(url)
##          
##          r = requests.get(url)
##          while r.status_code != 200:
##              r = requests.get(url)
##  
##          j = r.json()
##  
##          dfs.append(pd.DataFrame(data = j[1:], columns = j[0]))
##          # print(dfs[-1])
##  
##      df = pd.concat(dfs).rename(columns = readable)
##  
##      for c in cols:
##          df[c] = pd.to_numeric(df[c], errors = "coerce")
##          if "Frac" in description[c]: df[c] = np.round(df[c] / 100., 3)
##  
##      df[["state", "county", "tract"] + cols].to_csv(table.lower() + "_data.csv", index = False, header = False, na_rep = "")
##  
##      ##  sys.exit()

    sql_con = psycopg2.connect(database = 'census', user = user, password = passwd,
                               host = "saxon.harris.uchicago.edu", port = 5432)
    cur = sql_con.cursor()
    print("Got cursor.")

    # delete it if it already exists.
    # try: cur.execute("DROP TABLE {};".format(table))
    # except psycopg2.ProgrammingError: print("nothing to delete!!")

    print("Creating table.")
    ## create a table with an appropriate schema
    varis = ", ".join(["{} float".format(c) for c in cols])
    print("CREATE TABLE {} (state smallint, county smallint, tract int, {});".format(table, varis))
    cur.execute("CREATE TABLE {} (state smallint, county smallint, tract int, {});".format(table, varis))

    for v in cols:
        print("COMMENT ON COLUMN {}.{} IS '{}';".format(table, v, description[v]))
        cur.execute("COMMENT ON COLUMN {}.{} IS '{}';".format(table, v, description[v]))

    ## copy the data into the table.
    with open(table.lower() + "_data.csv") as f: 
        cur.copy_from(f, table, columns = ["state", "county", "tract"] + cols, sep = ",", null = '')
        sql_con.commit()

    print("ALTER TABLE {} ADD PRIMARY KEY (state, county, tract);".format(table))
    cur.execute("ALTER TABLE {} ADD PRIMARY KEY (state, county, tract);".format(table))

    if table is "acssf5y2015":
        cur.execute("alter table acssf5y2015 add column total_vap integer, add column black_vap integer, add column hispanic_vap integer;")
        cur.execute("update acssf5y2015 set total_vap = b05003_008e + b05003_019e, black_vap = b05003b_008e + b05003b_019e, hispanic_vap = b05003i_008e + b05003i_019e;")

        cur.execute("""UPDATE census_tracts_2015 AS tr SET
		         pop      = acs.b01001_001e,
		         black    = acs.b01001b_001e, 
		         hispanic = acs.b01001i_001e,
		         vap      = acs.total_vap,
		         bvap     = acs.black_vap,
		         hvap     = acs.hispanic_vap 
		       FROM acssf5y2015 AS acs WHERE 
		         tr.state = acs.state AND tr.county = acs.county AND tr.tract = acs.tract;""")

    if "ACSProfile5Y" in table: acsprofile5y(table)

def acsprofile5y(table):

   table = table.lower()

   print("""
PLEASE RUN::

psql census << EOD
ALTER TABLE {} ADD log_mhi      FLOAT;
ALTER TABLE {} ADD recent_immig FLOAT;
ALTER TABLE {} ADD geoid        BIGINT;

UPDATE {} SET geoid        = state::bigint * 1000000000 + county * 1000000 + tract;
UPDATE {} SET log_mhi      = CASE WHEN mhi > 0 THEN ROUND(LN(mhi)::numeric, 3) ELSE NULL END;
UPDATE {} SET recent_immig = CASE WHEN total_pop > 0 THEN ROUND((xrecent_immig / total_pop)::numeric, 3) ELSE NULL END;

ALTER TABLE {} DROP COLUMN xrecent_immig;
ALTER TABLE {} DROP COLUMN pop16_up;
EOD""".format(table, table, table, table, table, table, table, table))

        
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--table", default = "ACSProfile5Y2015")
    args = parser.parse_args()
    main(**vars(args))

