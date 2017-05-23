// g++ sip.cxx -I ~/anaconda3/include/ -std=c++11 -o sip

const bool debug = true;

#include <boost/graph/astar_search.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/random.hpp>
#include <boost/random.hpp>
#include <sys/time.h>
#include <vector>
#include <list>
#include <iostream>
#include <fstream>
#include <math.h>    // for sqrt, atan2
#include <utility>
#include <functional>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/segment.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/ring.hpp>
#include <boost/geometry/multi/geometries/multi_linestring.hpp>
#include <boost/geometry/multi/geometries/multi_polygon.hpp>
#include <boost/geometry/algorithms/union.hpp>

#include <boost/geometry/geometries/adapted/boost_polygon.hpp>

#include "astar_test.h"

#define DEBUG_ME cerr << __FILE__ << "::" << __LINE__ << "\t" << __FUNCTION__ << endl


using std::cout;
using std::cerr;
using std::endl;

// using namespace boost;
// using namespace std;

namespace bgeo  = boost::geometry;
namespace bpoly = boost::polygon;

typedef bpoly::point_data<double> bp_pt;
// typedef bpoly::segment_data<double> bp_seg;

namespace bgmod = boost::geometry::model;
typedef bgmod::d2::point_xy<double> bg_pt;
typedef bgmod::multi_point<bp_pt>   bg_mpt; // note BP point (not BG), for Voronoi.
typedef bgmod::polygon<bp_pt>       bg_poly;
typedef bgmod::segment<bp_pt>       bg_seg;
typedef bgmod::linestring<bp_pt>    bg_lstr;
typedef bgmod::ring<bp_pt>          bg_ring;
typedef bgmod::box<bp_pt>           bg_box;

typedef bgmod::multi_linestring<bg_lstr> bg_mlstr;
typedef bgmod::multi_polygon<bg_poly>    bg_mpoly;

// euclidean distance heuristic
template <class Graph, class CostType, class LocMap>
class distance_heuristic : public boost::astar_heuristic<Graph, CostType> {

  public:
    distance_heuristic(LocMap& l, int goal) : m_location(l), m_goal(goal) {}

    CostType operator()(int u) { return bgeo::distance(m_location[m_goal], m_location[u]); }

  private:
    LocMap m_location;
    int m_goal;
};


struct found_goal {}; // exception for termination

// visitor that terminates when we find the goal
template <class Vertex>
class astar_goal_visitor : public boost::default_astar_visitor {
  public:
    astar_goal_visitor(Vertex goal) : m_goal(goal) {}
    template <class Graph>
    void examine_vertex(Vertex u, Graph& g) {
      // return;
      if(u == m_goal) throw found_goal();
    }
  private:
    Vertex m_goal;
};

bool within_angles(float dx, float dy,
                   float thlA_lo, float thlA_hi,
                   float thlB_lo, float thlB_hi) {

  float thAB_A = atan2(dy, dx);
  if (thlA_hi < thlA_lo) thlA_hi += 2 * M_PI;
  if (thAB_A  < thlA_lo) thAB_A += 2 * M_PI;
  if (thlA_hi < thAB_A) return false;

  float thAB_B = thAB_A - M_PI;
  if (thlB_hi < thlB_lo) thlB_hi += 2 * M_PI;
  if (thAB_B  < thlB_lo) thAB_B += 2 * M_PI;
  if (thlB_hi < thAB_B) return false;

  return true;

}

bool segments_cross(bg_seg &A, bg_seg &B) {

  if (A.first  == B.first)  return false;
  if (A.first  == B.second) return false;
  if (A.second == B.first)  return false;
  if (A.second == B.second) return false;

  if (bgeo::intersects(A, B)) return true;
}


int main(int argc, char **argv)
{

  // specify some types
  typedef boost::property<boost::edge_weight_t, float> edge_weight_prop_t;
  typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS, boost::no_property, edge_weight_prop_t> mygraph_t;
  typedef mygraph_t::vertex_descriptor vertex;

  bg_mpt mpt;

  bool concave = false; int nmline = 0; bg_mlstr mline;

  bg_poly poly;
  boost::geometry::read_wkt(pa7, poly);

  float dxn, dyn, dxp, dyp;
  std::vector<std::pair<float, float> > angles;
  std::vector<int> concave_borders;

  mygraph_t g;

  bg_ring exr = bgeo::exterior_ring(poly);
  int npts = exr.size();

  // Check if the last point is concave?
  // Not necesary if everyone closes their polygons...
  dxn = exr[npts-1].x() - exr[npts-2].x();
  dyn = exr[npts-1].y() - exr[npts-2].y();
  dxp = exr[npts-1].x() - exr[0].x();
  dyp = exr[npts-1].y() - exr[0].y();
  if (dxp * dyn < dyp * dxn) concave = true;

  for (int pti = 0; pti < npts; pti++) {

    int ptp = (npts+pti-1)%npts;
    dxn = exr[(pti+1)%npts].x() - exr[pti].x();
    dyn = exr[(pti+1)%npts].y() - exr[pti].y();

    if (dxp * dyn < dyp * dxn) {

      bgeo::append(mpt, exr[pti]);

      angles.push_back(std::make_pair(atan2(dyp, dxp),
                                      atan2(dyn, dxn)));

      concave_borders.push_back(concave);
      concave = true; // pay it forward.

    } else concave = false;

    dxp = -dxn; dyp = -dyn;
  }

  std::vector<bg_seg> segments;
  auto ipt  = boost::begin(bgeo::exterior_ring(poly));
  auto ipte = boost::end(bgeo::exterior_ring(poly));
  segments.push_back(bg_seg(*ipt, *(ipte-1)));
  for (; ipt != ipte; ++ipt) {
    segments.push_back(bg_seg(*ipt, *(ipt+1)));
  }

  // For the short cache set.
  int hcount = 0;
  const int nhopeful = 10;
  int hopefuls[nhopeful];
  for (int ii = 0; ii < nhopeful; ii++) hopefuls[ii] = ii;

  int ncpts = mpt.size();
  for (int pti = 0; pti < ncpts;  pti++) {
    
    cout << "point i=" << pti << " (ncpts=" << ncpts << ")" << endl;
    for (int ptj = pti+1; ptj < ncpts; ptj++) {

      float dx = mpt[ptj].x() - mpt[pti].x();
      float dy = mpt[ptj].y() - mpt[pti].y();
      float th = atan2(dy, dx);

      // Check if the path would *immediately* 
      // leave the polygon, on either side.
      if (!within_angles(dx, dy, angles[pti].first, angles[pti].second, 
                                 angles[ptj].first, angles[ptj].second)) continue;

      // Otherwise look for explicit *crossings*
      // (not intersections, since it's ok to hit 
      // the target point itself!
      bg_seg segAB(mpt[pti], mpt[ptj]);
      bool intersected = false;

      // Maintain a short list of the most-recent
      // segments crossed: 2.5x improvement.
      for (int hi = 0; hi < nhopeful; hi++) {
        if (segments_cross(segAB, segments[hopefuls[hi]])) {
          intersected = true;
        }
      }
      if (intersected) continue;

      int nseg = 0;
      for (auto seg : segments) {
        if (segments_cross(segAB, seg)) {
          hopefuls[hcount] = nseg;
          hcount = (hcount+1)%nhopeful;

          intersected = true;
          break;
        }

        nseg++;
      }
      if (intersected) continue;

      if (debug) {
        bg_lstr line;
        line.push_back(mpt[pti]);
        line.push_back(mpt[ptj]);
        mline.resize(++nmline);
        bgeo::append(mline[nmline-1], line);
      }

      float dist = ::sqrt(dx*dx + dy*dy);
      add_edge(pti, ptj, edge_weight_prop_t(dist), g);
    }

    if (concave_borders[pti]) {

      int ptj = (ncpts+pti-1)%ncpts;

      // Keep the line here, for the sanity plot.
      if (debug) {
        bg_lstr line;
        line.push_back(mpt[pti]);
        line.push_back(mpt[ptj]);
        mline.resize(++nmline);
        bgeo::append(mline[nmline-1], line);
      }

      // This is the actual work -- making the graph.
      float dx = mpt[ptj].x() - mpt[pti].x();
      float dy = mpt[ptj].y() - mpt[pti].y();
      float dist = ::sqrt(dx*dx + dy*dy);
      add_edge(pti, ptj, edge_weight_prop_t(dist), g);
    }

  }

  int found = 0;
  bg_box env;
  bgeo::envelope(poly, env);
  float min_x = env.min_corner().x();
  float min_y = env.min_corner().y();
  float max_x = env.max_corner().x();
  float max_y = env.max_corner().y();
  
  srand (time(0));
  while (found < 2) {

    float r1 = static_cast <float> (rand()) / RAND_MAX;
    float r2 = static_cast <float> (rand()) / RAND_MAX;

    bp_pt pt(min_x + (max_x - min_x) * r1,
             min_y + (max_y - min_y) * r2);

    cout << bgeo::within(pt, poly) << "   (x,y)=(" << pt.x() << "," << pt.y() << ")" << endl;
    if (!bgeo::within(pt, poly)) continue;

    for (int pti = 0; pti < ncpts+found;  pti++) {

      bg_lstr line;
      line.push_back(mpt[pti]);
      line.push_back(pt);

      if (bgeo::within(line, poly)) {
        mline.resize(++nmline);
        bgeo::append(mline[nmline-1], line);

        float dx = mpt[pti].x() - pt.x();
        float dy = mpt[pti].y() - pt.y();
        float dist = ::sqrt(dx*dx + dy*dy);
        
        add_edge(ncpts+found, pti, edge_weight_prop_t(dist), g);
      }

    }

    bgeo::append(mpt, pt);

    found++;
  }

  // pick random start/goal
  vertex start = ncpts;
  vertex goal  = ncpts+1;
  bg_lstr path;
  
  std::vector<mygraph_t::vertex_descriptor> p(num_vertices(g));
  std::vector<float> d(num_vertices(g));
  try {
    // call astar named parameter interface
    astar_search(g, start, distance_heuristic<mygraph_t, float, bg_mpt>(mpt, goal),
                           boost::predecessor_map(&p[0]).distance_map(&d[0]).visitor(astar_goal_visitor<vertex>(goal)));
  
  } catch(found_goal fg) { // found a path to the goal
    std::list<vertex> shortest_path;
    for(vertex v = goal;; v = p[v]) {
      shortest_path.push_front(v);
      if(p[v] == v) break;
    }
    cout << "Shortest path from " << start << " to "
         << goal << ": ";
    std::list<vertex>::iterator spi = shortest_path.begin();
    cout << start;
    bgeo::append(path, mpt[*spi]);
    for(++spi; spi != shortest_path.end(); ++spi) {
      cout << " -> " << *spi;
      bgeo::append(path, mpt[*spi]);
    }
    cout << endl << "Total travel time: " << d[goal] << endl;
  }


  // Declare a stream and an SVG mapper
  std::ofstream svg("path.svg");

  float aspect = sqrt((max_x - min_x)/(max_y - min_y));
  boost::geometry::svg_mapper<bp_pt> mapper(svg, 800*aspect, 800/aspect);

  // Add geometries such that all these geometries fit on the map
  mapper.add(poly);
  mapper.add(mpt);
  mapper.add(mline);
  mapper.add(path);

  // Draw the geometries on the SVG map, using a specific SVG style
  // Chicago grey=(128,122,117) and red=(139,0,33)
  mapper.map(poly,  "fill-opacity:0.6;fill:rgb(128,122,117);");
  mapper.map(mpt,   "fill-opacity:1.0;fill:rgb(139,0,33);stroke:rgb(139,0,33);stroke-width:0.10");
  mapper.map(mline, "fill-opacity:1.0;fill:rgb(139,0,33);stroke:rgb(139,0,33);stroke-width:0.10");
  mapper.map(path,  "fill-opacity:1.0;fill:rgb(0,0,255);stroke:rgb(0,0,255);stroke-width:3");

  return 0;
  
}
