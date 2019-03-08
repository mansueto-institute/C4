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
#include <iomanip>
#include <stdlib.h>

// #include "boost/geometry/geometry.hpp"
// #include "boost/range/adaptor/reversed.hpp"

#include <math.h>
#include "c4/Miniball.hpp"

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

#include "c4/topology.h"

// #include "ehrenburg.hpp"



namespace c4 {

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

  const double VOR_SCALE = 1e1;

  bool pop_compare(Region* r1, Region* r2);
  bool id_compare(Region* r1, Region* r2);

  template <typename T1, typename T2>
  bool compare_first(std::pair<T1, T2> a, std::pair<T1, T2> b);

  template <typename T1, typename T2>
  bool compare_second(std::pair<T1, T2> a, std::pair<T1, T2> b);

  bool cellp_set_len_compare(std::unordered_set<Cell*> s1, std::unordered_set<Cell*> s2);

  // Helper functions for std::any_of calls.
  auto uo_set_nempty = [](std::unordered_set<Cell*> i) { return !i.empty(); };
  auto uo_set_empty  = [](std::unordered_set<Cell*> i) { return  i.empty(); };

  // Helper function for std::none_of calls.
  bool maxed_or_empty(unsigned int max, std::vector<std::unordered_set<Cell*> >& graphs, 
                                        std::vector<std::unordered_set<Cell*> >& to_add);

  template <typename T> int sign(T x);
  template <typename T> T clip(const T& n, T clipval);

  // Possible objective functions.
  enum ObjectiveMethod {DISTANCE_A, DISTANCE_P, INERTIA_A, INERTIA_P, HULL_A, HULL_P, POLSBY, REOCK, EHRENBURG,
                        PATH_FRAC, AXIS_RATIO, MEAN_RADIUS, DYN_RADIUS, HARM_RADIUS, ROHRBACH, EXCHANGE, POLSBY_W};

  // Potential radii for d2 or dist, and get_circle_coords
  enum RadiusType      {EQUAL_AREA, EQUAL_AREA_POP, EQUAL_CIRCUMFERENCE, SCC, LIC, HULL, POWER};

  // Base class for a single node.
  class Cell {

    public:

      Cell();
      Cell(const Cell&);
      Cell(int i, int p, double x, double y,
           double a, weight_map wm, float edge_perim, bool is_split); // , std::string mp_wkt);

      // Distance squared to a point.
      float d2(float xi, float yi) { return (x-xi)*(x-xi) + (y-yi)*(y-yi); }

      // Distance squared to another cell.
      float d2(Cell* ci) { return (x-ci->x)*(x-ci->x) + (y-ci->y)*(y-ci->y); }

      // Distance to a point
      float dist(float xi, float yi) { return sqrt(d2(xi, yi)); }
      
      // Distance to another cell.
      float dist(Cell* ci) { return sqrt(d2(ci->x, ci->y)); }
      
      // Convert Cell adjacency map to pointers to other cells.
      void adjacency_to_pointers(std::vector<Cell*>&);

      // Add an edge (side of a polygon), between two node ids.
      void add_edge(int edge_id, int nodea, int nodeb);

      // Convert node IDs to pointers to nodes, for edges.
      void node_ids_to_pointers(std::vector<Node*>&);

      // void check_internal_connectedness();

      // Disconnect two nodes.  This must be done when merging them.
      void disconnect();

      // Merge two nodes -- for example, if there is a cut vertex.
      void merge(Cell*);


      // At each node, cycle around the incoming edges, and report if the next edge is in the same region.
      // This is used by Region::get_node_ring to determine the boundary of a region.
      bool next_edge_in_region(Node* node, int start_edge_id, Cell*& next_cell, Edge*& next_edge, bool CW, int region_i);

      // Neighbors connected is just a convenience function for neighbor_sets().
      bool neighbors_connected(bool QUEEN, bool FAST);

      // This is _the key_ function.  It determines (quickly if possible, then with a full search)
      //  the connected components that the neigbors are members of of.  
      // If there are multiple such connected components, then the node is a cut vertex.
      int  neighbor_sets(std::vector<std::unordered_set<Cell*> >& graphs, unsigned int max_merge = 0, bool QUEEN = false, bool FAST = false);

      // Find components, less than a certain size,
      //   that would be separated from this graph if this vertex were removed.
      int  neighbor_strands(std::unordered_set<Cell*>& strand, unsigned int max_merge, bool QUEEN);

      // Region IDs of this cell's neighbors.
      std::unordered_set<int> neighbor_regions(bool QUEEN = false);

      // Are cell's neighbors of eachother?
      bool is_neighbor(Cell* c) { for (auto& n : nm) if (n.first == c) return true; return false; }

      // Does this cell have a neighbor in region r?
      bool touches_region(int r) { for (auto& n : nm) if (n.first->region == r) return true; return false; }

      // Get an edge to either a foreign region or the edge of the "universe."
      int  get_ext_edge_idx();
      int  get_edge_idx(int id);

      // Does this cell have this edge?
      bool has_edge(int edge_id);

      // Get an edge pointer for an edge id.
      Edge* get_edge(int edge_id);

      // std::pair<Node*, Node*> get_cw_node(int reg);

      int region;
      int id, pop;
      double x, y;  
      double area;
  
      // IDs or pointers to neighbors.
      weight_map wm;
      neighbor_map nm;

      // Fixed quantities.
      float edge_perim;
      bool  is_univ_edge;

      // Census tracts are occasionally split -- 
      //   this can really complicate the topologies.
      bool  is_split;
      bool  split_neighbor;

      bp_pt    pt;

      std::map<int, Cell*> dijkstra_step;

      // Cells can be removed from the graph
      std::vector<int> enclaves_and_islands;

      // This is the topology of the PERIMETERS, 
      //   rather than the connectivity between cells (tracts).
      std::vector<Edge> edges;
      std::set<Node*> nodes;

  };

  class Region {

    public:

      Region();
      Region(int rid);
      Region(int rid, Cell*);

      void reindex(int rid);

      void add_cell(Cell* c, bool b);
      void remove_cell(Cell* c, bool b);

      void add_cell_int_ext_neighbors(Cell* c);
      void remove_cell_int_ext_neighbors(Cell* c);

      // Check if the region is contiguous, optionally returning 
      // a set of disconnected components.
      bool contiguous();
      bool contiguous(std::vector<std::unordered_set<Cell*> > &graphs);

      // Evaluate an objective function possibly with 
      //  the addition or removal of another cell.
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

      // Distances and distances squared.  Can specify how to calculate the region position.
      float d2(float x, float y, RadiusType rt = RadiusType::EQUAL_AREA);
      float d2(Region* r, RadiusType rt = RadiusType::EQUAL_AREA);
      float d2(Cell* c, RadiusType rt = RadiusType::EQUAL_AREA) { return d2(c->x, c->y, rt); }
      float dist(float x, float y, RadiusType rt = RadiusType::EQUAL_AREA) { return sqrt(d2(x, y, rt)); }
      float dist(Cell*c, RadiusType rt = RadiusType::EQUAL_AREA) { return dist(c, rt); }

      // For moment of inertia calculations, use the parallel axis theorem 
      //  to do this efficiently.
      double inertia_parallel_axis(double I0, double x0, double y0, double w0, double xc, double yc, double wc);


      // Fundamental properties of the region must be constantly maintained.
      int id;
      int pop;
      size_t ncells;
      double area;

      // Areal center
      double xctr, yctr;

      // Population center
      double xpctr, ypctr;

      // Border length
      double sumw_border;

      // Moments of inertia.
      double Ip, Ia;
      
      // Topology (nodes and edges) have been loaded.
      bool has_topo;

      // The Universe creates the cells; 
      //   the region tracks its contents,
      //   the internal borders, and the external borders.
      std::unordered_set<Cell*> cells;
      std::unordered_set<Cell*> ext_borders;
      std::unordered_set<Cell*> int_borders;

      // Tolerance, centers and radii for LIC and SCC methods.
      double eps_scc, x_scc, y_scc, r2_scc;
      double eps_lic, x_lic, y_lic, r2_lic;

      void  reset_borders();
      float update_scc(Cell* add, Cell* sub, bool UPDATE);
      float update_lic(Cell* add, Cell* sub, bool UPDATE);


      // Center and power for power diagram method.
      double x_pow, y_pow, r2_pow;

      // Very complicated functions to extract the border
      //   from the nodes/edges of the cells.
      bool make_ring();
      void divert_ring_at_cell(Cell* c, bool CW);

      // Get a circle center and radius, for a given radius type.
      std::pair<std::pair<float, float>, float> get_circle_coords(RadiusType rt);

      // Floating point coordinates of the node ring, if more than 2 cells.
      std::vector<std::pair<float, float> > get_point_ring();
      std::vector<std::pair<float, float> > hull(bool IB);
      void get_node_ring(std::vector<Node*>& nr, Cell* cell = 0, int i = -1);

      // Principle components for length to width ratio.
      std::pair<float, float> update_pca(Cell* add = 0, Cell* sub = 0, bool vec = false, bool UPDATE = false);
      float pca0, pca1;
      std::pair<float, float> pca_vec0, pca_vec1;

      bg_mpoly poly;

      // Boost geometry multipoint of points in region.
      bg_mpt   mpt;

      // Convex Hull stuff...
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

      // Total population
      int pop;

      // Number of cells and region.
      unsigned long ncells, nregions;

      // Target population
      double target;

      // Total iterations run.
      int total_iterations;

      // New cells, edges and nodes.
      void add_cell(Cell);
      void add_edge(int cell_id, int edge_id, int nodea, int nodeb);
      void add_node(int node_id, float x, float y);
      void add_node_edge(int node_id, int edge_id);

      // Coordinates of all circles.
      std::pair<std::pair<float, float>, float> get_circle_coords(size_t rid, RadiusType rt);

      // Coordinates of the point runt or the hull.
      std::vector<std::pair<float, float> > get_point_ring(size_t rid);
      std::vector<std::pair<float, float> > hull(size_t rid, int IB = 0);

      // For testing: explicitly move a cell to a region.
      void add_cell_to_region(int cid, size_t rid);

      // Return the cell to region mapping (python dict)
      std::map<int, int> cell_region_map();

      // Return cells on the border of a region.
      // Can be internal or external --
      //   not equal, since internal includes edges of the universe.
      std::vector<int> border_cells(bool EXT, int rid);

      // Get list of cells that are "wards" of other cells.
      std::vector<int> clipped_cells();

      void connect_graph();

      // Consolidate parts of the graph that are less
      //  than a fraction of a single region.
      void trim_graph(float max_frac = 0.9);

      // At the single cell level, 
      int  merge_strands(Cell* c, float max_frac = 0.9);

      // When moving strands, contiguity can get broken.
      //   If this happens, force it.
      void force_contiguity(int rid, bool verbose = false);

      bool loaded_topo;
      bgraph dijkstra_graph;

      // Convert adjacency matrix and edge node networks to pointers.
      void adjacency_to_pointers();
      void node_ids_to_pointers();

      // This is used for the path fraction.
      void build_dijkstra_graph();

      std::vector<int> do_dijkstra(int start_id, int end_id);

      std::vector<Cell*>   cells;
      std::vector<Region*> regions;
      std::vector<Node*>   nodes;

      // Reset everything.
      void assign_to_zero();

      // Create new regions with power diagrams.
      void power_restart(int seed, int niter, float tol, int = false);

      // Split one region, merge two others, and restart.
      void split_restart(int seed, ObjectiveMethod om);
      bool split_region(int r = 0, float a = -1, bool connect = true, float margin = 0.1);
      bool merge_regions(int rA, int rB);

      // Initialize from the split line algo.
      void split_line_init();

      // Initialize with random draws.
      void rand_init(int s);

      // Grow from random draws by kmeans or random.
      void grow_kmeans(int popgrow);
      void grow_random(int seed = 0);

      // Run power diagrams.
      void iterate_power(float tol, int niter, int reset = false, int verbose = false); // , int popgrow = false);
      void voronoi_classify();
      // void center_power_cells();
      
      float best_solution_val;
      float best_tolerance_val;
      int   iterations_since_improvment;
      std::map<int, int> best_solution;
      // std::vector<std::pair<float, std::map<int, int> > > best_solutions;
      
      // Keep track of the best contender seen so far.
      void update_best_solutions(ObjectiveMethod omethod, float tol, bool verbose = false);
      std::map<int, int> get_best_solution() { return best_solution; }
      void load_best_solution() { load_partition(best_solution); }

      // Be able to load an arbitrary starting configuration.
      void load_partition(std::map<int, int> reg_map);

      // Single type iterator, using Chen and Rodden method.
      void iterate(int niter, float tol, int r);

      // Actually run the optimization.
      float ALPHA;
      bool oiterate(ObjectiveMethod omethod, int niter, float llh_tol, float cut_tol, int conv_iter, int seed, int r, int verbose);
      bool greedy(Region* rit, ObjectiveMethod omethod, float tol, float best_move = 0, bool random = false, int r = -1, bool verbose = false);
      bool greedy_evaluate(Region* r, Cell* b, float tol, ObjectiveMethod omethod, float& best_move, Cell*& b_opt_c, std::unordered_set<Cell*>& opt_strands, bool verbose = false);

      int  RANDOM;
      std::mt19937 mersenne;

      int  TRADE;
      bool trade(Region* a, Region* b, ObjectiveMethod om);

      // Tabu list settings.
      size_t TABU_LENGTH;
      std::deque<Cell*>    tabu;
      void set_tabu(Cell* c) { tabu.push_front(c); while (tabu.size() > TABU_LENGTH) tabu.pop_back(); }
      bool is_tabu(Cell* c)  { return TABU_LENGTH && find(tabu.begin(), tabu.end(), c) != tabu.end(); }
      bool is_tabu_strand(std::unordered_set<Cell*> &s) { for (auto c : s) if (is_tabu(c)) return true; return false;}
      bool is_uninit_strand(std::unordered_set<Cell*> &s);

      size_t DESTRAND_MIN, DESTRAND_MAX;
      bool transfer_strand(std::unordered_set<Cell*>& strand);
      int  destrand(int mini = 1, size_t maxi = 1e9, float ctol = 1);


  };


}



