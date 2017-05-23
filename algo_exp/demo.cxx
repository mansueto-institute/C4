// g++ demo.cxx -I ~/anaconda3/include/ -std=c++11 -o demo

#include <iostream>
#include <fstream>

#include <unordered_set>
#include <algorithm>
#include <random>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/multi/geometries/multi_polygon.hpp>
#include <boost/geometry/algorithms/union.hpp>

#include <boost/polygon/voronoi.hpp>
#include <boost/geometry/geometries/adapted/boost_polygon.hpp>

#include "astar_test.h"



using std::cout;
using std::cerr;
using std::endl;

namespace bgeo  = boost::geometry;
namespace bpoly = boost::polygon;

typedef bpoly::point_data<double> bp_pt;
typedef bpoly::segment_data<double> bp_seg;

namespace bgmod = boost::geometry::model;
typedef bgmod::d2::point_xy<double> bg_pt;
typedef bgmod::multi_point<bp_pt>   bg_mpt; // note BP point (not BG), for Voronoi.
typedef bgmod::polygon<bp_pt>       bg_poly;
typedef bgmod::linestring<bp_pt>    bg_lstr;
typedef bgmod::ring<bp_pt>          bg_ring;
typedef bgmod::multi_polygon<bg_poly> bg_mpoly;


int main() {

    bg_poly poly;
    boost::geometry::read_wkt(pa7, poly);

    bp_pt ctr;
    bgeo::centroid(poly, ctr);

    bg_poly hull;
    bgeo::convex_hull(poly, hull);

    bg_lstr line;
    std::vector<bp_seg> sdata2;
    for (auto ipt  = boost::begin(bgeo::exterior_ring(poly));
              ipt != boost::end(bgeo::exterior_ring(poly)); ++ipt) {
      line.push_back(*ipt);
      if (ipt != boost::end(bgeo::exterior_ring(poly))) sdata2.push_back(bp_seg(*ipt, *(ipt+1)));
    }

    bpoly::voronoi_diagram<double> vd;
    construct_voronoi(sdata2.begin(), sdata2.end(), &vd);

    double d2;
    bp_pt  vd_ctr;
    double max_dist = boost::numeric::bounds<double>::lowest();

    for (auto it : vd.vertices()) {

      bp_pt pt(it.x(), it.y());
      d2 = bpoly::euclidean_distance(sdata2[it.incident_edge()->cell()->source_index()], pt);

      if (max_dist < d2 && bgeo::within(pt, poly)) {
        max_dist = d2;
        bgeo::set<0>(vd_ctr, it.x());
        bgeo::set<1>(vd_ctr, it.y());
      }
    } 

    // Declare the point_circle strategy
    int pp_circ = 360;
    boost::geometry::strategy::buffer::point_circle point_strategy(pp_circ);

    // Declare other strategies
    boost::geometry::strategy::buffer::distance_symmetric<double> vor_ds(max_dist);
    boost::geometry::strategy::buffer::join_round join_strategy(pp_circ);
    boost::geometry::strategy::buffer::end_round end_strategy(pp_circ);
    boost::geometry::strategy::buffer::point_circle circle_strategy(pp_circ);
    boost::geometry::strategy::buffer::side_straight side_strategy;

    // Create the buffer of a multi point
    bgeo::model::multi_polygon<bg_poly> vd_buff;
    bgeo::buffer(vd_ctr, vd_buff, vor_ds, side_strategy,
                                  join_strategy, end_strategy, circle_strategy);

    // Create the buffer of a multi point
    float epR = bgeo::perimeter(poly) / (2 * M_PI);
    boost::geometry::strategy::buffer::distance_symmetric<double> epR_ds(epR);
    bgeo::model::multi_polygon<bg_poly> ep_buff;
    bgeo::buffer(ctr, ep_buff, epR_ds, side_strategy,
                               join_strategy, end_strategy, circle_strategy);

    
    // Create the buffer of a multi point
    float eaR = sqrt(bgeo::area(poly) / M_PI);
    boost::geometry::strategy::buffer::distance_symmetric<double> eaR_ds(eaR);
    bgeo::model::multi_polygon<bg_poly> ea_buff;
    bgeo::buffer(ctr, ea_buff, eaR_ds, side_strategy,
                              join_strategy, end_strategy, circle_strategy);


    // Declare a stream and an SVG mapper
    std::ofstream svg("demo.svg");
    bgeo::svg_mapper<bp_pt> mapper(svg, 1200, 1200, "width=\"120%\" height=\"120%\"");

    // Add geometries such that all these geometries fit on the map
    mapper.add(poly);
    mapper.add(hull);
    mapper.add(vd_ctr);
    mapper.add(vd_buff);
    mapper.add(ea_buff);
    mapper.add(ep_buff);

    // Draw the geometries on the SVG map, using a specific SVG style
    mapper.map(hull,    "fill-opacity:0;stroke:rgb(255,0,0);stroke-width:3,2", 5);
    mapper.map(poly,    "fill-opacity:0.5;fill:rgb(128, 122, 117);");
    mapper.map(vd_buff, "fill-opacity:0;stroke:rgb(255,0,0);stroke-width:3,2", 5);
    mapper.map(ea_buff, "fill-opacity:0;stroke:rgb(0,255,0);stroke-width:3,2", 5);
    mapper.map(ep_buff, "fill-opacity:0;stroke:rgb(0,0,255);stroke-width:3,2", 5);

    return 0;

}

