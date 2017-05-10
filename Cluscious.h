#include <execinfo.h>
#include <assert.h>

#include <deque>
#include <vector>
#include <iterator>
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

#include <armadillo>

#define DEBUG_ME std::cerr << __FILE__ << "  " << __FUNCTION__ << "::" << __LINE__ << std::endl

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

  const double VOR_SCALE = 1e4;

  bool pop_compare(Region* r1, Region* r2);
  bool id_compare(Region* r1, Region* r2);
  bool cellp_set_len_compare(std::unordered_set<Cell*> s1, std::unordered_set<Cell*> s2);

  auto uo_set_nempty = [](std::unordered_set<Cell*> i) { return !i.empty(); };
  auto uo_set_empty  = [](std::unordered_set<Cell*> i) { return  i.empty(); };

  bool maxed_or_empty(unsigned int max, std::vector<std::unordered_set<Cell*> >& graphs, 
                                        std::vector<std::unordered_set<Cell*> >& to_add);

  template <typename T> int sign(T x);
  template <typename T> T clip(const T& n, T clipval);

  enum ObjectiveMethod {DISTANCE_A, DISTANCE_P, INERTIA_A, INERTIA_P, HULL_A, HULL_P, POLSBY, REOCK, EHRENBURG,
                        PATH_FRAC, AXIS_RATIO, MEAN_RADIUS, DYN_RADIUS, HARM_RADIUS, ROHRBACH, EXCHANGE, POLSBY_W};

  enum RadiusType      {EQUAL_AREA, EQUAL_AREA_POP, EQUAL_CIRCUMFERENCE, SCC, LIC, HULL, POWER};

  class Cell {

    public:

      Cell();
      Cell(const Cell&);
      Cell(int i, int p, double x, double y,
           double a, weight_map wm, float edge_perim, bool is_split); // , std::string mp_wkt);

      float d2(float xi, float yi) { return (x-xi)*(x-xi) + (y-yi)*(y-yi); }
      float d2(Cell* ci) { return (x-ci->x)*(x-ci->x) + (y-ci->y)*(y-ci->y); }
      float dist(float xi, float yi) { return sqrt(d2(xi, yi)); }
      float dist(Cell* ci) { return sqrt(d2(ci->x, ci->y)); }

      void add_edge(int edge_id, int nodea, int nodeb);
      void adjacency_to_pointers(std::vector<Cell*>&);
      void node_ids_to_pointers(std::vector<Node*>&);
      void check_internal_connectedness();
      void merge(Cell*);

      bool next_edge_in_region(Node* node, int start_edge_id, Cell*& next_cell, Edge*& next_edge, bool CW, int region_i);
      bool neighbors_connected(bool QUEEN, bool FAST);
      int  neighbor_sets(std::vector<std::unordered_set<Cell*> >& graphs, unsigned int max_merge = 0, bool QUEEN = false, bool FAST = false);
      int  neighbor_strands(std::unordered_set<Cell*>& strand, unsigned int max_merge, bool QUEEN);

      std::unordered_set<int> neighbor_regions(bool QUEEN = false);
      bool is_neighbor(Cell* c) { for (auto& n : nm) if (n.first == c) return true; return false; }
      bool touches_region(int r) { for (auto& n : nm) if (n.first->region == r) return true; return false; }

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

      float edge_perim;
      bool  is_univ_edge;
      bool  is_split;
      bool  split_neighbor;

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

      double obj(ObjectiveMethod omethod, Cell* add = 0, Cell* sub = 0, bool verbose = false);

      double obj_distance   (Cell* add, Cell* sub, float xx, float yy);
      double obj_inertia_a  (Cell* add = 0, Cell* sub = 0, bool verbose = false);
      double obj_inertia_p  (Cell* add = 0, Cell* sub = 0, bool verbose = false);
      double obj_reock      (Cell* add = 0, Cell* sub = 0, bool verbose = false);
      double obj_hull       (Cell* add = 0, Cell* sub = 0, bool verbose = false);
      double obj_hull_p     (Cell* add = 0, Cell* sub = 0, bool verbose = false);
      double obj_polsby     (Cell* add = 0, Cell* sub = 0, bool verbose = false);
      double obj_path_frac  (Cell* add = 0, Cell* sub = 0, bool verbose = false);
      double obj_ehrenburg  (Cell* add = 0, Cell* sub = 0, bool verbose = false);
      double obj_axis_ratio (Cell* add = 0, Cell* sub = 0, bool verbose = false);
      double obj_mean_radius(Cell* add = 0, Cell* sub = 0, bool verbose = false);
      double obj_dyn_radius (Cell* add = 0, Cell* sub = 0, bool verbose = false);
      double obj_harm_radius(Cell* add = 0, Cell* sub = 0, bool verbose = false);
      double obj_rohrbach   (Cell* add = 0, Cell* sub = 0, bool verbose = false);
      double obj_exchange   (Cell* add = 0, Cell* sub = 0, bool verbose = false);

      float d2(float x, float y, RadiusType rt = RadiusType::EQUAL_AREA);
      float d2(Region* r, RadiusType rt = RadiusType::EQUAL_AREA);
      float d2(Cell* c, RadiusType rt = RadiusType::EQUAL_AREA) { return d2(c->x, c->y, rt); }
      float dist(float x, float y, RadiusType rt = RadiusType::EQUAL_AREA) { return sqrt(d2(x, y, rt)); }
      float dist(Cell*c, RadiusType rt = RadiusType::EQUAL_AREA) { return dist(c, rt); }

      double inertia_parallel_axis(double I0, double x0, double y0, double w0, double xc, double yc, double wc);


      int id;
      int pop;
      size_t ncells;
      double area;
      double xctr, yctr;
      double xpctr, ypctr;

      double sumw_border;

      double Ip, Ia;
      
      bool has_topo;

      // The Universe creates the cells.
      std::unordered_set<Cell*> cells;
      std::unordered_set<Cell*> ext_borders;
      std::unordered_set<Cell*> int_borders;

      double eps_scc, x_scc, y_scc, r2_scc;
      double eps_lic, x_lic, y_lic, r2_lic;
      float update_scc(Cell* add, Cell* sub, bool UPDATE);
      float update_lic(Cell* add, Cell* sub, bool UPDATE);

      double x_pow, y_pow, r2_pow;

      void make_ring();
      void divert_ring_at_cell(Cell* c, bool CW);
      std::pair<std::pair<float, float>, float> get_circle_coords(RadiusType rt);
      std::vector<std::pair<float, float> > get_point_ring();
      void get_node_ring(std::vector<Node*>& nr, Cell* cell = 0, int i = -1);

      std::pair<float, float> update_pca(Cell* add = 0, Cell* sub = 0, bool vec = false, bool UPDATE = false);
      float pca0, pca1;
      std::pair<float, float> pca_vec0, pca_vec1;

      bg_mpoly poly;
      bg_mpt   mpt;

      double   ch_area;
      bg_poly  ch_poly;
      bp_pt    ch_pt;
      float    ch_x, ch_y;
      std::vector<Node*> node_ring;

  };

  class Universe {

    public:

      Universe();
      Universe(int);
      int pop;
      unsigned long rcount, ncells, nregions;
      double target;

      int total_iterations;

      void add_cell(Cell);
      void add_edge(int cell_id, int edge_id, int nodea, int nodeb);
      void add_node(int node_id, float x, float y);
      void add_node_edge(int node_id, int edge_id);

      std::pair<std::pair<float, float>, float> get_circle_coords(size_t rid, RadiusType rt);
      std::vector<std::pair<float, float> > get_point_ring(size_t rid);
      void add_cell_to_region(int cid, size_t rid);

      std::map<int, int> cell_region_map();
      std::vector<int> border_cells(bool EXT, int rid);
      std::vector<int> clipped_cells();

      void connect_graph();

      void trim_graph(float max_frac = 0.9);
      int  merge_strands(Cell* c, float max_frac = 0.9);


      bool loaded_topo;
      bgraph dijkstra_graph;
      void adjacency_to_pointers();
      void node_ids_to_pointers();
      void build_dijkstra_graph();

      std::vector<int> do_dijkstra(int start_id, int end_id);

      std::vector<Cell*>   cells;
      std::vector<Region*> regions;
      std::vector<Node*>   nodes;

      void assign_to_zero();
      bool split_region(int r = 0, float a = -1, bool connect = true);
      bool merge_regions(int rA, int rB);
      void split_line_init();
      void rand_init(int s);
      void grow_kmeans(int popgrow);
      void grow_random(int seed = 0);

      void iterate_power(float tol, int niter, int reset = false); // , int popgrow = false);
      void voronoi_classify();
      // void center_power_cells();

      void load_partition(std::map<int, int> reg_map);
      void iterate(int niter, float tol, int r);

      float ALPHA;
      void oiterate(ObjectiveMethod omethod, int niter, float tol, int seed, int r, int verbose);
      bool greedy(Region* rit, ObjectiveMethod omethod, float tol, float best_move = 0, bool random = false, int r = -1, bool verbose = false);
      bool greedy_evaluate(Region* r, Cell* b, float tol, ObjectiveMethod omethod, float& best_move, Cell*& b_opt_c, std::unordered_set<Cell*>& opt_strands, bool verbose = false);

      int  RANDOM;
      std::mt19937 mersenne;

      int  TRADE;
      bool trade(Region* a, Region* b, ObjectiveMethod om);

      size_t TABU_LENGTH;
      std::deque<Cell*>    tabu;
      void set_tabu(Cell* c) { tabu.push_front(c); while (tabu.size() > TABU_LENGTH) tabu.pop_back(); }
      bool is_tabu(Cell* c)  { return TABU_LENGTH && find(tabu.begin(), tabu.end(), c) != tabu.end(); }
      bool is_tabu_strand(std::unordered_set<Cell*> s) { for (auto c : s) if (is_tabu(c)) return true; return false;}

      size_t DESTRAND_MIN, DESTRAND_MAX;
      void transfer_strand(std::unordered_set<Cell*>& strand);
      int  destrand(int mini, size_t maxi);


  };


}



