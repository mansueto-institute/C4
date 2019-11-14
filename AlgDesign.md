## Algorithm Design

The aim of the software is to optimize compact, contiguous partitions of states.
Ultimately, efficient optimization 
  means efficiently calculating changes to compactness scores
  from the addition or removal of a single cell, rather than from scratch, each time.

For example, the IPQ (Polsby-Popper measure)
  relates the area of a shape to its perimeter.
It is trivial to add or remove one cell worth of area, 
  but keeping track of perimeters is a little harder.
Adding one cell,
  some parts of the old perimeter will now fall _within_ the shape,
  and those the same parts of the _new_ cell's own perimeter will too!
So one must track all of the adjacencies and add the right bits of cells'
  perimeters each time.
Every algorithm has considerations like this one.

The software is composed of three c++ classes.
These are exposed to python through cython.  

### Classes

#### Universe

This class is, in the districting context, the state.
The state contains a fixed number of regions (districts), and it owns/instantiates 
  all regions and cells (tracts or block groups).
To avoid copies, everyone else just passes pointers around.

It is initiated with simply the number of regions.
It then adds cells, and the edges and nodes between them.
It transforms these into its own (somewhat involved) interlinked set of pointers.
This can be seen in `run.py`, but I would advise against editing this.

Finally, it checks the contiguity of the state itself,
  and connects it if it is not already connected.
This process is important, for instance, for Hawaii.
(These connections can also be specified, as for [California](https://github.com/JamesSaxon/C4/blob/master/shapes/ca_conx.csv).)

The Universe is the object that runs the optimization,
  and with which any user should interact (except from creating single cells).
The three regionalization functions are 
  `oiterate()`, `iterate_power()`, and `split_line_init()`.

#### Region

Region is used mainly to evaluate objective functions.
It saves at all times its area, its various centroids (population or area-weighted), etc.
This entails keeping track of cells along its border -- both internal and external.
It also has a number of convenience functions for distances,
  and can run contiguity checks.

#### Cells

Cells are where a lot of magic lives.
Cells [know](https://github.com/JamesSaxon/C4/blob/master/c4.h#L203)
  their region, their ID, their population and location and -- 
  most important all of their neighbors (to which they have pointers).
If cells are the cut vertex to a portion of the graph,
  that portion will be removed and the cell will record
  a listing of the cells from "beyond the cut."

### Contiguity Constraints

Contiguity constrains are enforced using integer programming.

#### Connecting the Graph

* Function: `Universe::connect_graph()`

The graph is is initially connected with `Universe::connect_graph()`.
This algorithm starts from a star graph for a single cell.
It then creates another one.  If the two graphs overlap, combined them.  Continue.

If there are multiple graphs at the end of this procedure, 
  find the closest cells between two graphs by Euclidean distance.
Draw an edge between these, and join the graph.

Note that these initially-disconnected pieces are 
  likely to be removed again (subsumed, really), by `Universe::trim_graph()`,
  since they cannot be moved without moving the cell to which they are connected.
But this trimming applies to both disconnected islands
  and to enclaves -- like small county seats in broad Texas counties.

#### Preserving Region Connectivity

* Key Function: `neighbor_sets()` (called by `neighbors_connected()` and `neighbor_strands()`)

The requirement is this:
  moving a cell from one region to another 
  must leave its neighbors from its old regions in connected subgraph.
Others (Cho and Liu) have enforced this quite simply by requiring that the neighbors themselves be connected,
  but this only works because they both preclude holes and use Queen weights.
I perform the broader search, so that the neighbors left behind my be connected by cells that the moving cell did not touch.
This is done by checking if neighbors are on different subgraphs, with `Cell:neighbor_sets()`.
This is possibly the most important and most complicated function in the code.
It is commented thoroughly, but it is NON-TRIVIAL!

### Optimization Routines

#### Greedy Optimizer
* Key Function: `Universe::oiterate()` and `Universe::greedy()`.
`Universe::oiterate()` loops over regions, grabbing the best cell to annex, if such a cell exists.
It forwards most of its work to `Universe::greedy()` (which also includes GRASP and tabu functionality),
  which in turn towards `Universe::greedy_evaluate()`.
This is where the population difference is calculated
  and where `Region::obj()` is called with cells to add or subtract.
If a move would result in the (best yet) change in the objective,
  I check if its neighbors are on separate subgraphs, with `neighbor_strands()`.

#### Split-Restart
* Code: `Universe::split_restart()`
Split-Restart is a means of "whacking" the system to start over, without completely starting from scratch.
The worst-performing district is split in two, and two other random districts are joined together.

#### Power 
* Code: `Universe::iterate_power()` and `Universe::voronoi_classify()`
Iterate power slowly changes the effective lambda parameter (`r2_pow`)
  and the power diagram centers (`x/y_pow`) of the regions.
It then runs `voronoi_classify()` to set region membership.

#### Split-Line
* Code: `Universe::split_line_init()`
This function spins an unit vector in steps of pi / 20,
  and divides the cells of a region in two (`Universe::split_region()`) 
  according to their dot product with that vector.
It chooses the angle that results in the smallest perimeter between the two new regions,
  and "finalizes" that split.
This search is very, very slow, but it only runs once per state, so it's OK.

### Objective Functions

There is one objective for each (greedy-optimized) method.
These are all run through `Region::obj()`, with a switch on an enum.
Each of these 
* `obj_distance`: Theis is à la Chen & Roddn.  Used both for population and area-weighted distances, with different centers.
* `obj_inertia_a/p`: Areal or population moment of inertia.  Changes are calculated using the parallel axis theorem.
* `obj_reock`: AKA smallest circumscribing circle (SCC) method.  This is calculated using Bernd Gärtner's incredible Miniball algorithm.  But there is a copy here, to line things up, so it's still not super fast.  The SCC is only recalculated if a cell to be removed is on the boundary or a cell to be added is beyond its radius.
* `obj_hull/_p`: Hulls and population hulls are recalculated using Boost, but only if the new cell would affect the border.
* `obj_polsby`: Polsby is calculated by tracking exactly which parts of each cell's perimeter are shared with its neighbors.  Every move of a cell moves some parts of a boundary into the shape, and exposes other parts of the moved cell.  This is done efficiently by tracking the internal and external neighbors of a region at all times.  Cells also have pointers to their neighbors and know the perimeter they share with them.
* `obj_polsby_w` Similar algorithm, but with perimeters weighted along regimes.
* `obj_path_frac`: This is the fraction of shortest paths between two residents of a district that are contained in the district itself.  This is calculated by caching the Dijkstra shortest paths for the entire state, before running.  To efficiently save all of the paths, I just record first step of each path, since _the next cell_ can tell you where to go from there.)  This is still very slow.  Note that other, apparently closely related algorithms, like any involving the shortest paths between points on perimeters are entirely infeasible.  Doing so requires calculating the entire visibility graph for each region from scratch for every change, instead of caching it at the beginning.
* `obj_ehrenburg`: AKA Largest Inscribed Circle (LIC).  This is entirely my own algorithm, and I quite like it.  First, create a polygon from the nodes of the external "ring" of the region (`get_node_ring()` -- itself, very complicated). Then run Boost's polygon Voronoi method to get the Voronoi diagram (and its points).  Then the LIC is the circle whose center is the node of the Voronoi diagram that is (a) within the ring of the region and (b) furthest from its source points (it must be equidistant from all of the source points -- that's what a Voronoi diagram _is_).
* `obj_axis_ratio`: This is calculated using armadillo from the PCA.  
* `obj_mean_radius`: Adding or removing a cell affects its center, so move that first -- this is just a weighted average.  Then recalculate all of the distances to the modified center.
* `obj_dyn_radius`: Same as above but with distances squared and a new normalization.
* `obj_harm_radius`: Ditto with 1R.
* `obj_rohrbach`: Distances to the perimeter.  I make the approximation of the centroid of the nearest border cell (instead of the nearest border _line_ or _node_).  Very similar to IPQ, we have to track the border cells.  Adding or removing a cell changes who's on the border.  Then just calculate the distances to the edge.   
* `obj_exchange`: Intersection of equal area circle with the polygon.  Calculate the modified center _C_ and radius _R_ of the equal area circle.  Then check which cells have are more than _R_ from _C_!

### Cython Wrapping
* Code: `pyc4.pyx` and `setup.py`

This wraps mainly the Universe class, which is the one that users should interact with 
  (since Universes own the other classes).
Only the functions to create a cell or region are otherwise exposed.
This simply maps the c++ classes and enums to python.

`setup.py` is the build script.

### Python Scripts

Generating districts requires very specific geometries and topologies,
  and so these data are included with the software.
These data are loaded using python scripts that 
  users should NOT manipulate.
 
There are a few files of note:
* `run.py`: runs the code.  `./run.py -h` will show all options.  
* `pyc4_helper.py`: these are plotting functions and data caching functions.
* `ps_query.py`: contains all of the SQL scripts for retrieving the data.

