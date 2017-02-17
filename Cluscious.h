#include <vector>
#include <map>
#include <unordered_set>
#include <algorithm>
#include <limits>
#include <iostream>
#include <stdlib.h>
#include <cassert>

// #include "boost/geometry/geometry.hpp"
// #include "boost/range/adaptor/reversed.hpp"

#include <math.h>
#include "Miniball.hpp"

#include <iostream>
#include <fstream>

#include <boost/geometry.hpp>
#include <boost/geometry/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>

#include <boost/geometry/algorithms/union.hpp>
#include <boost/geometry/algorithms/difference.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/point.hpp>
// #include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/multi/geometries/multi_polygon.hpp>

#include <boost/polygon/voronoi.hpp>
#include <boost/geometry/geometries/adapted/boost_polygon.hpp>

// #include "ehrenburg.hpp"


#define DEBUG_ME cerr << __FILE__ << "::" << __LINE__ << "\t" << __FUNCTION__ << endl

namespace Cluscious {

  using std::cout;
  using std::endl;
  using std::cerr;
  namespace bgeo  = boost::geometry;
  namespace bpoly = boost::polygon;
  
  typedef bpoly::point_data<double> bp_pt;
  typedef bpoly::segment_data<double> bp_seg;
  
  namespace bgmod = boost::geometry::model;
  // typedef bgmod::d2::point_xy<double>   bg_pt;
  typedef bgmod::multi_point<bp_pt>     bg_mpt; // note BP point_data (not BG point_xy), for Voronoi.
  typedef bgmod::polygon<bp_pt>         bg_poly;
  typedef bgmod::linestring<bp_pt>      bg_lstr;
  typedef bgmod::ring<bp_pt>            bg_ring;
  typedef bgmod::multi_polygon<bg_poly> bg_mpoly;



  // forward declare for typedef, sort
  class Cell;
  class Region;

  typedef std::map<int, double> weight_map;
  typedef std::map<Cell*, double> neighbor_map;

  // Miniball types for Reock weights.
  typedef std::list<std::vector<double> >::const_iterator PointIterator; 
  typedef std::vector<double>::const_iterator CoordIterator;
  typedef Miniball::Miniball <Miniball::CoordAccessor<PointIterator, CoordIterator> > MB;

  bool pop_compare(Region* r1, Region* r2);
  bool id_compare(Region* r1, Region* r2);
  bool cellp_set_len_compare(std::unordered_set<Cell*> s1, std::unordered_set<Cell*> s2);

  auto uo_set_nempty = [](std::unordered_set<Cell*> i) { return !i.empty(); };

  int sign(double x);

  enum ObjectiveMethod {DISTANCE_A, DISTANCE_P, INERTIA_A, INERTIA_P, HULL_A, POLSBY, REOCK, EHRENBURG, POLSBY_W, PATH_FRAC};

  class Cell {

    public:

      Cell();
      Cell(const Cell&);
      Cell(int i, int p, double x, double y,
           double a, weight_map wm, std::string mp_wkt);

      float d2(float xi, float yi) { return (x-xi)*(x-xi) + (y-yi)*(y-yi); }

      void contiguity_to_neighbors(std::vector<Cell*>&);
      void merge(Cell*);

      bool neighbors_connected();
      int neighbor_sets(std::vector<std::unordered_set<Cell*> >& graphs, bool for_merging);

      int region;
      int id, pop;
      double x, y;  
      double area;
      weight_map wm;
      neighbor_map nm;

      bp_pt    pt;
      bg_mpoly poly;

      std::vector<int> linked;

  };

  class Region {

    public:

      Region();
      Region(int rid, Cell*);
      Region(int rid, double x, double y);

      void add_cell(Cell*, bool);
      void remove_cell(Cell*, bool);

      double obj(ObjectiveMethod omethod, Cell* add, Cell* sub, bool verbose);

      double obj_distance (Cell* add, Cell* sub, float xx, float yy);
      double obj_inertia_a(Cell* add, Cell* sub, bool verbose);
      double obj_inertia_p(Cell* add, Cell* sub, bool verbose);
      double obj_reock    (Cell* add, Cell* sub, bool verbose);
      double obj_hull     (Cell* add, Cell* sub, bool verbose);
      double obj_polsby   (Cell* add, Cell* sub, bool verbose);

      int id;
      int pop;
      int ncells;
      double area;
      double xctr, yctr;
      double xpctr, ypctr;

      double sumw_border;

      // The Universe creates the cells.
      std::unordered_set<Cell*> cells;
      std::unordered_set<Cell*> ext_borders;
      std::unordered_set<Cell*> int_borders;

      double x_mb, y_mb, r2_mb, eps_mb;
      float update_miniball(Cell* add, Cell* sub, bool UPDATE);

      bg_mpoly poly;
      bg_mpt   mpt;

      double   ch_area;
      bg_poly  ch_poly;

  };

  class Universe {

    public:

      Universe();
      Universe(int);
      int pop;
      unsigned long rcount, nregions;
      double target;

      void add_cell(Cell);
      int  get_ncells();

      std::map<int, int> cell_region_map();
      std::vector<int> border_cells(int rid);
      std::vector<int> clipped_cells();

      void connect_graph();

      void trim_graph();
      int  merge_strands(Cell* c, int max);


      void contiguity_to_neighbors();

      std::vector<Cell*>   cells;
      std::vector<Region*> regions;

      void rand_districts(int s);

      void grow_kmeans(int popgrow);
      void iterate(int niter, float tol, int r);

      void oiterate(ObjectiveMethod omethod, int niter, float tol, float alpha, int r, int verbose);


  };

}



