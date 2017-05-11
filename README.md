# cluscious

Cluscious is a collection of three c++ classes exposed to python through cython.
The goal is to perform fast, iterative, contiguity-preserving optimization
  of many of the compactness objective functions found in the gerrymandering literature.



To build:

Some parts need armadillo (which needs OpenBlas):

OpenBLAS:
download https://github.com/xianyi/OpenBLAS/zipball/master
Then 
`make && sudo make install`

(Yes, it's that easy!?)

Armadillo
Linux -- trivial:
sudo apt-get install libarmadillo-dev libarmadillo6 libarmadillo6-dbgsym

Mac -- Download
http://arma.sourceforge.net/download.html

```
cmake . 
make 
sudo make install

g++ -I/usr/local/include/ example1.cpp -o example1 -O2 -larmadillo
```

Then for Cluscious, it's:

```
python setup.py build_ext --inplace
```

Then to test, try
```
./test.py PA 2
```

I have included one set of test shapefiles;
  normally I pull this from postgres.


Thus far, I have implemented: 
* areal and population distance
* areal and population moment of inertia (normalizeed)
* circumscribing circles with Miniball
* convex hull with boost

I have adapted the largest inscribed circle from mapbox/polylabel,
  to boost (plus some simplifications, using boost functions), in `ehrenburg.hpp`.
Alternatively, in `vor.cxx`, I used the boost Voronoi class to find the nodes
  of the Voronoi diagram of the line segments of the polygon;
  one can then check each of these nodes to find the one furthest from the external ring.
Neither of them are sufficiently fast.
The next step, therefore, is to maintain a single ring
  consisting of the centroids of each of the border districts.

