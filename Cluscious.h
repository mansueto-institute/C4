#include <execinfo.h>
#include <assert.h>

#include <vector>
#include <map>
#include <unordered_set>
#include <algorithm>
#include <limits>
#include <iostream>
#include <stdlib.h>

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
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/multi/geometries/multi_polygon.hpp>

#include <boost/polygon/voronoi.hpp>
#include <boost/geometry/geometries/adapted/boost_polygon.hpp>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/astar_search.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/random.hpp>

#define DEBUG_ME std::cerr << __FILE__ << "  " << __FUNCTION__ << "::" << __LINE__ << "\t" << __FUNCTION__ << std::endl

#include "cluscious/topology.h"

// #include "ehrenburg.hpp"



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

  // Types for the Dijkstra graph.
  typedef boost::property<boost::edge_weight_t, float> edge_weight_prop_t;
  typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS, boost::no_property, edge_weight_prop_t> bgraph;
  typedef bgraph::vertex_descriptor bvertex;


  // forward declare for typedef, sort
  class Cell;
  class Region;

  typedef std::map<int, double> weight_map;
  typedef std::vector< std::pair<Cell*, double> > neighbor_map;

  // Miniball types for Reock weights.
  typedef std::list<std::vector<double> >::const_iterator PointIterator; 
  typedef std::vector<double>::const_iterator CoordIterator;
  typedef Miniball::Miniball <Miniball::CoordAccessor<PointIterator, CoordIterator> > MB;

  bool pop_compare(Region* r1, Region* r2);
  bool id_compare(Region* r1, Region* r2);
  bool cellp_set_len_compare(std::unordered_set<Cell*> s1, std::unordered_set<Cell*> s2);

  auto uo_set_nempty = [](std::unordered_set<Cell*> i) { return !i.empty(); };

  bool maxed_or_empty(unsigned int max, std::vector<std::unordered_set<Cell*> >& graphs, 
                                        std::vector<std::unordered_set<Cell*> >& to_add);

  int sign(double x);

  enum ObjectiveMethod {DISTANCE_A, DISTANCE_P, INERTIA_A, INERTIA_P, HULL_A, POLSBY, REOCK, EHRENBURG, POLSBY_W, PATH_FRAC};

  class Cell {

    public:

      Cell();
      Cell(const Cell&);
      Cell(int i, int p, double x, double y,
           double a, weight_map wm, bool is_univ_edge); // , std::string mp_wkt);

      float d2(float xi, float yi) { return (x-xi)*(x-xi) + (y-yi)*(y-yi); }
      float dist(float xi, float yi) { return sqrt(d2(xi, yi)); }

      void add_edge(int edge_id, int nodea, int nodeb);
      void adjacency_to_pointers(std::vector<Cell*>&);
      void node_ids_to_pointers(std::vector<Node*>&);
      void merge(Cell*);

      bool next_edge_in_region(Node* node, int start_edge_id, Cell*& next_cell, Edge*& next_edge, bool CW);
      bool neighbors_connected();
      int  neighbor_sets(std::vector<std::unordered_set<Cell*> >& graphs, unsigned int max_merge, bool QUEEN);
      int  neighbor_strands(std::unordered_set<Cell*>& strand, int dreg, unsigned int max_merge, bool QUEEN);

      int  get_ext_edge_idx();
      int  get_edge_idx(int id);
      bool has_edge(int edge_id);
      Edge* get_edge(int edge_id);

      std::pair<Node*, Node*> get_cw_node(int reg);

      int region;
      int id, pop;
      double x, y;  
      double area;
      weight_map wm;
      neighbor_map nm;

      bool is_univ_edge;

      bp_pt    pt;

      std::map<int, Cell*> dijkstra_step;

      std::vector<int> enclaves_and_islands;

      std::vector<Edge> edges;
      std::set<Node*> nodes;

  };

  class Region {

    public:

      Region();
      Region(int rid);
      Region(int rid, Cell*);

      void add_cell(Cell* c, bool b);
      void remove_cell(Cell* c, bool b);

      void add_cell_int_ext_neighbors(Cell* c);
      void remove_cell_int_ext_neighbors(Cell* c);

      double obj(ObjectiveMethod omethod, Cell* add, Cell* sub, bool verbose);

      double obj_distance (Cell* add, Cell* sub, float xx, float yy);
      double obj_inertia_a(Cell* add, Cell* sub, bool verbose);
      double obj_inertia_p(Cell* add, Cell* sub, bool verbose);
      double obj_reock    (Cell* add, Cell* sub, bool verbose);
      double obj_hull     (Cell* add, Cell* sub, bool verbose);
      double obj_polsby   (Cell* add, Cell* sub, bool verbose);
      double obj_path_frac(Cell* add, Cell* sub, bool verbose);

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

      void make_ring();
      void add_cell_to_ring(Cell* c);
      std::vector<std::pair<float, float> > get_point_ring();
      void get_node_ring(std::vector<Node*>& nr, Cell* cell = 0, int i = -1);
      // void sub_cell_to_ring(Cell* c);

      bg_mpoly poly;
      bg_mpt   mpt;

      double   ch_area;
      bg_poly  ch_poly;
      std::vector<Node*> node_ring;

  };

  class Universe {

    public:

      Universe();
      Universe(int);
      int pop;
      unsigned long rcount, nregions;
      double target;

      void add_cell(Cell);
      void add_edge(int cell_id, int edge_id, int nodea, int nodeb);
      void add_node(int node_id, float x, float y);
      void add_node_edge(int node_id, int edge_id);
      int  get_ncells();

      std::vector<std::pair<float, float> > get_point_ring(size_t rid);
      void add_cell_to_region(size_t rid, int cid);

      std::map<int, int> cell_region_map();
      std::vector<int> border_cells(bool EXT, int rid);
      std::vector<int> clipped_cells();

      void connect_graph();

      void trim_graph();
      int  merge_strands(Cell* c, int max);


      bgraph dijkstra_graph;
      void adjacency_to_pointers();
      void node_ids_to_pointers();
      void build_dijkstra_graph();

      std::vector<int> do_dijkstra(int start_id, int end_id);

      std::vector<Cell*>   cells;
      std::vector<Region*> regions;
      std::vector<Node*>   nodes;


      void rand_init(int s);
      void grow_kmeans(int popgrow);

      void load_partition(std::map<int, int> reg_map);
      void iterate(int niter, float tol, int r);

      void oiterate(ObjectiveMethod omethod, int niter, float tol, float alpha, int r, int verbose);

      void transfer_strand(std::unordered_set<Cell*>& strand, Region* source, Region* dest);


  };


}



