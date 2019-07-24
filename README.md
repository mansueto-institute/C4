# <img src="img/c4_logo.png" width=25px> C4: Contiguity-Constrained Clustering in C++

The goal of C4 is to perform fast, iterative, contiguity-preserving optimization
  of many of the compactness objective functions found in the gerrymandering literature.

The implemented algorithms include:
* Isoperimeter Quotient: e.g., [Polsby & Popper](https://digitalcommons.law.yale.edu/ylpr/vol9/iss2/6)
* Moments of Inertia: [Weaver & Hess](http://dx.doi.org/10.2307/794769)
* Largest-Inscribed Circle: [Ehrenburg](https://babel.hathitrust.org/cgi/pt?id=mdp.39015076673584;view=1up;seq=41)
* Smallest Circumscribing Circle: [Reock](http://dx.doi.org/10.2307/2109043)
* Convex Hull Area or Population, e.g., [Hofeller & Grofman](http://www.socsci.uci.edu/~bgrofman/B48-Comparing-the-Compactness.pdf)
* Various Radii: e.g., [Frolov 1975, for a good historical review](http://dx.doi.org/10.1080/00385417.1975.10640104)
* Exchange: [Angel et al](http://dx.doi.org/10.1111/j.1541-0064.2009.00304.x)
* Power Diagrams, e.g., [Fryer & Holden](http://dx.doi.org/10.1086/661511)
* Path Fraction ("Bizareness"): [Chambers & Miller](http://dx.doi.org/10.1561/100.00009022)
* Distance assignment: [Chen & Rodden](http://dx.doi.org/10.1561/100.00012033)
* Split-line: [Forrest](http://dx.doi.org/10.1177/000276426400800407)

To do this, C4 includes three c++ classes:
  (1) a universe, (2) a region, and (3) cells.
In districting parlance, this means states, legislative districts, and Census tracts (or block groups or blocks).
The Universe class (mainly) is exposed to python through cython.
All plotting and most data management happens in python.

## Running C4 as a Docker Container (Simple!)

To facilitate use and replication, all of the dependences and the built software are included in a docker container,
  hosted on [DockerHub](https://hub.docker.com/r/jamessaxon/c4).
The container was generated using the included [Dockerfile](Dockerfile) and is about 2 GB.
Running from scratch is as simple as:
```
docker run -v $(pwd)/res/:/C4/res/ -e STATE=pa -e SEED=2 -e METHOD=POWER --rm -it  jamessaxon/c4:replication
```

This means:
* Run the c4 software as in this image (`docker run [...] jamessaxon/c4:replication`)
* Mount the local directory called `res/` to `/C4/res/` in the container.  Results written here will be available when the jobs completes. You **must make the `res/` directory!!**
* Environment variables / arguments: simulate districts for the `STATE` of Pennsylvania (USPS code `pa`), with `SEED` of 2 (any number, but < 1000 will format better), using the `POWER` diagram `METHOD`.  You can also specify which maps to draw with the `SHADING` variable (`-e SHADING=all` is _all_ of them).
   * The possible methods are `POWER DIST RADII IPQ CIRCLES HULL_P HULL_A INERTIA AXIS SPLIT PATH_FRAC`.  Several of these run several methods in sequence.  For instance `CIRCLES` runs the `exchange`, `reock` and `ehrenburg` methods, and `RADII` includes `rohrbach`, `harm_radius`, `mean_radius`, and `dyn_radius `.
   * The SHADING options are: `district` (just colors), `target` (ratio to target population), `density` (show population centers), `scores` (show spatial scores), `counties` (overlay county geometries), and `all` or `none`.  The default is `none`.
* Remove (`--rm`) the container when it exits.
* Run interactively and allow input (`-it`).

Of course, you can also run interactively with `/bin/bash` and just `cd` to `C4` to use `run.py` with all of its arguments.
In that case, skip to [Running C4](#running-c4), below.

For large-scale jobs, I run C4 as a Docker container on AWS.
The Dockerfile changes very slightly, with AWS keys and `s3` tools (see [DockerfileAWS](DockerfileAWS)).
The scripts for [building the container](docker_build.sh) and [launching jobs](aws_launch.sh) 
  are also in this directory (though they won't run, without dependencies or passwords!).

## Running C4 with a Local Install

### Building

Some components require Armadillo, which in turn requires OpenBlas.  Note that if you create a python environment with geopandas, OpenBLAS will come for free.  So you may not need to to do this.
* OpenBLAS: First [download](https://github.com/xianyi/OpenBLAS/zipball/master) it, then `make && sudo make install`.  (Yes, it's that easy!)
* Armadillo:
  * Linux: `sudo apt-get install libarmadillo-dev libarmadillo6 libarmadillo6-dbgsym`
  * Mac [download](http://arma.sourceforge.net/download.html) then `cmake . && make && sudo make install`
  
You will also need all of the compiled and python packages listed in the [Dockerfile](Dockerfile):
* Compiled: `libboost-all-dev` `libgeos-dev`, `libgdal-dev`, `python3-gdal`, `gdal-bin`
* Python: `cython`, `matplotlib`, `fiona`, `pysal`, `geopandas`, `psycopg2`
GEOS and GDAL can be finnicky with respect to the OS, so this installation is on the user.  
However, Anaconda has made it much easier.  On relatively modern (few year-old) Macs, `conda install geopandas pysal cython psycopg2 libboost` seems to give you everything you need.
  
Then to build **C4**, it's just `python setup.py build_ext --inplace`.

### Running C4

Setting up all of the shapefiles, topologies, and voting records from scratch is pretty involved, but these are all included in the package.  To run C4, I suggest checking out `run.py` which will show you its options (`-h`), or `run_iter.sh` where you can find the default settings for any of the methods.

For example, to run power diagrams for Pennsylvania, do: 
```
./run.py -s pa -i power:100000 -t 0.01  -x300 -l0 -c100 --power_restart --print_init -w pa/power/s300 -m power
```
This means: run Pennsylvania, for up to 100k iterations initialized through the power diagram method, using a tolerance of 1%, a seed of 300.  Run this for 100 cycles, and restart after completion using the power restart method.  Write to file.  The method is power.

This is a little different from most methods, since power diagrams do not use the standard greedy optimizer.  More typical is:
```
./run.py -s pa -m hull_p -t 0.01 -x300 -n10000 --conv_iter 500 -c20  --destrand_min 3 --destrand_max 50 --allow_trades
```
This means use the hull population method.  Run up to 10000 iterations total, and stop after 500 iterations with no improvement.
Restart the search 20 times.
Do remove "strands" from the regions (larger than 3, but smaller than 50 cells).

Alternatively, you can just accept my defaults, and do
```
export STATE=pa; export SEED=300; export METHOD=IPQ; ./run_iter.sh
```
which will run Pennsylvania for seed 300 with the Isoperimeter quotient method.
This is the default script that the Docker container runs.
The possible methods are: `POWER`, `DIST`, `RADII`, `IPQ`, `CIRCLES`, `HULL_P`, `HULL_A`, `INERTIA`, `AXIS`, `SPLIT`, `PATH_FRAC`.
Several of these options run multiple methods in series.

## Outputs

Four types of files will be written to `res/` (your local directory).
Note again that, running with docker, the directory must be mounted to a local file.
1. JSON files containing a summary of the simulation will be written to `res/json/[state usps]_[method]_s[seed]_c[cycle].json`.  These files contain a summary of the entire run: the tract to district assignment, the spatial parameters of the districts, the partisan voting (if available for that state), race and ethnicity, voter balance (`PopulationDeviation`), the method used, and so forth (run `jq keys file.json` to see this).  These data 
2. CSV files containing simply the tract to district assignment.   This is just a two-column assignment: row number (equivalent to county + tract geoid, though perhaps a poor technical choice) and the district.  This will be written to `res/[state usps]/[method]/s[seed]/c[cycle]/final.csv`. 
3. GeoJSON files of the districting plans, `final.geojson`, in the same directory as (2).  These contain the basic plan and some of the information of (1).  They are in in the EPSG 4326 coordinate reference system, and display nicely on GitHub or Gist or in Leaflet etc.
4. Finally, `final_*.pdf` are static maps of the districting plans.  There will be one map for each shading method used.


## Browsing Maps:

If you're more interested in the results,
  just head over to my webspace to play with the outputs in an interactive map:

http://saxon.harris.uchicago.edu/redistricting_map/

