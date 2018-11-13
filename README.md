# C4: Contiguity-Constrained Clustering in C++

C4 is a collection of three c++ classes exposed to python through cython.
The goal is to perform fast, iterative, contiguity-preserving optimization
  of many of the compactness objective functions found in the gerrymandering literature.

This includes:
* Isoperimeter Quotient: e.g., [Polsby & Popper](https://digitalcommons.law.yale.edu/ylpr/vol9/iss2/6)
* Moments of Inertia: [Weaver & Hess](10.2307/794769)
* Largest-Inscribed Circle: [Ehrenburg](https://babel.hathitrust.org/cgi/pt?id=mdp.39015076673584;view=1up;seq=41)
* Smallest Circumscribing Circle: [Reock](10.2307/2109043)
* Convex Hull Area or Population, e.g., [Hofeller & Grofman](http://www.socsci.uci.edu/~bgrofman/B48-Comparing-the-Compactness.pdf)
* Various Radii: e.g., [Frolov 1975, for a good historical review](http://dx.doi.org/10.1080/00385417.1975.10640104)
* Exchange: [Angel et al](http://dx.doi.org/10.1111/j.1541-0064.2009.00304.x)
* Power Diagrams, e.g., [Fryer & Holden](http://dx.doi.org/10.1086/661511)
* Path Fraction ("Bizareness"): [Chambers & Miller](http://dx.doi.org/10.1561/100.00009022)
* Distance assignment: [Chen & Rodden](http://dx.doi.org/10.1561/100.00012033)
* Split-line: [Forrest](http://dx.doi.org/10.1177/000276426400800407)

## Build Instructions

Some components require Armadillo, which in turn requires OpenBlas:
* OpenBLAS: First [download](https://github.com/xianyi/OpenBLAS/zipball/master) it, then `make && sudo make install`.  (Yes, it's that easy!)
* Armadillo:
  * Linux: `sudo apt-get install libarmadillo-dev libarmadillo6 libarmadillo6-dbgsym`
  * Mac [download](http://arma.sourceforge.net/download.html) then `cmake . && make && sudo make install`
  
You will also need all of the compiled and python packages listed in the [Dockerfile](Dockerfile):
* Compiled: `libgeos-dev`, `libgdal-dev`, `python3-gdal`, `gdal-bin`
* Python `cython`, `matplotlib`, `fiona`, `pysal`, `geopandas`, `psycopg2`
GEOS and GDAL can be very finnicky WRT OS, so this is on you.
  
Then to build **C4**, it's just `python setup.py build_ext --inplace`.

#### Docker
For large-scale jobs, I run C4 as a Docker container on AWS.  You can check out the [Dockerfile](Dockerfile), as well as the [build](docker_build.sh) and [launch](aws_launch.sh) scripts I use for this.

## Running C4

Setting up all of the shapefiles, topologies, and voting records from scratch is pretty involved, but these are all included in the package.  To run C4, I suggest checking out `run.py` which will show you its options (`-h`), or `run_iter.sh` where you can find the default settings for any of the methods.

For example, to run power diagrams for Pennsylvania, do: 
```
./run.py -s pa -i power:100000 -t 0.01  -x300 -l0 -c100 --power_restart --print_init -w pa/power/s300 -m power
```
This means: run Pennsylvania, for up to 100k iterations initialized through the power diagram method, using a tolerance of 1%, a seed of 300.  Run this for 100 cycles, and restart after completion using the power restart method.  Write to file.  The method is power.

This is a little different from most methods, since power diagrams do not use the standard greedy optimizer.  More typical is:
```
./run.py -s pa -m hull_p -t 0.01 -x300 -n10000 -c20 --conv_iter 500  --destrand_min 3 --destrand_max 50 --allow_trades
```
This means use the hull population method.  Run up to 10000 iterations total, and stop after 500 iterations with no improvement.  Do remove "strands" from the regions (larger than 3, but smaller than 50 cells).

Aternatively, you can just accept my defaults, and do
```
export STATE=pa; export SEED=300; export METHOD=IPQ
```
which will run Pennsylvania for seed 300 with the Isoperimeter quotient method.

## Browsing Maps:

If you're more interested in the results, just head over to my webspace to play with the outputs in an interactive map:

http://saxon.harris.uchicago.edu/redistricting_map/

