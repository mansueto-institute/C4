#!/usr/bin/env python

# CREATE TABLE acs2015 AS SELECT state, county, tract FROM census_tracts_2015;

import cenpy as cpy
import pandas as pd

import psycopg2

import argparse

from netrc import netrc
user, acct, passwd = netrc().authenticators("harris")
user, acct, apikey = netrc().authenticators("census")

import sys

def main(v, table = "ACSSF5Y2015"):

    con = cpy.base.Connection("ACSSF5Y2015", apikey = apikey)
    
    cols = []
    for vi in v: cols.extend(con.varslike(vi))

    dfs = []
    for st in list(cpy.explorer.fips_table('state')['FIPS Code']):
            
        st = int(st)
        if st > 57: continue

        dfs.append(con.query(cols, geo_unit = 'tract:*', geo_filter = {'state':str(st)},
                                   apikey = apikey))

    df = pd.concat(dfs).fillna(0)

    print(list(df.columns))
    df[["state", "county", "tract"] + cols].to_csv("temp.csv", index = False, header = False)



    sql_con = psycopg2.connect(database = 'census', user = user, password = passwd,
                               host = "saxon.harris.uchicago.edu", port = 5432)

    cur = sql_con.cursor()

    # delete it if it already exists.
    try: cur.execute("DROP TABLE {};".format(table))
    except psycopg2.ProgrammingError: print("nothing to delete!!")

    ## create a table with an appropriate schema
    varis = ", ".join(["{} float".format(c) for c in cols])
    cur.execute("CREATE TABLE {} (state smallint, county smallint, tract int, {});".format(table, varis))
    print("CREATE TABLE {} (state smallint, county smallint, tract int, {});".format(table, varis))

    ## copy the data into the table.
    with open("temp.csv") as f: 
        cur.copy_from(f, table, columns = ["state", "county", "tract"] + cols, sep = ",")
        sql_con.commit()

        
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--vars", nargs = "+")
    parser.add_argument("--table", default = "ACSSF5Y2015")
    args = parser.parse_args()
    main(args.vars, args.table)

