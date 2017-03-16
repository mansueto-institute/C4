// g++ sip.cxx -I ~/anaconda3/include/ -std=c++11 -o sip

#include <iostream>
#include <fstream>

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

#include <boost/polygon/voronoi.hpp>
#include <boost/geometry/geometries/adapted/boost_polygon.hpp>


#include "ehrenburg.hpp"


#define DEBUG_ME cerr << __FILE__ << "::" << __LINE__ << "\t" << __FUNCTION__ << endl


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
typedef bgmod::segment<bp_pt>       bg_seg;
typedef bgmod::linestring<bp_pt>    bg_lstr;
typedef bgmod::ring<bp_pt>          bg_ring;

typedef bgmod::multi_linestring<bg_lstr> bg_mlstr;
typedef bgmod::multi_polygon<bg_poly>    bg_mpoly;


int main()
{
    bg_pt pt1(2, 1.3);

    bg_mpoly mpoly;
    bg_mpt mpt;

    // boost::geometry::read_wkt("MULTIPOLYGON(())", mpoly); 
    boost::geometry::read_wkt("MULTIPOLYGON(((2 0,2 8,3 8,3 2,4 8,6 8,11 4,8 10,9 11,15 7,17 12,11 13,15 17,18 17,20 10,21 17,23 17,24 10,26 17,29 17,33 13,27 12,29 7,35 11,36 10,33 4,38 8,40 8,41 2,41 8,42 8,42 0,39 0,38 5,37 0,34 0,31 3,34 8,29 5,26 8,23 8,22 13,21 8,18 8,15 5,10 8,13 3,10 0,7 0,6 6,5 0,2 0)))", mpoly);

    int pi(0);
    for (auto p : mpoly) {
      cout << "poly" << pi << ", area: " << bgeo::area(p) << endl;

      bg_ring exr = bgeo::exterior_ring(p);
      int npts = exr.size();

      float dxn, dyn;
      float dxp = exr[npts-1].x() - exr[0].x();
      float dyp = exr[npts-1].y() - exr[0].y();

      for (int pti = 0; pti < npts; pti++) {

        int ptn = (pti+1)%npts;
        int ptp = (npts+pti-1)%npts;
        dxn = exr[(pti+1)%npts].x() - exr[pti].x();
        dyn = exr[(pti+1)%npts].y() - exr[pti].y();

        cout << "   point " << pti << ", (x,y)=(" << exr[pti].x() << "," << exr[pti].y() << ")    "
             << "   next "  << ptn << ", (x,y)=(" << exr[ptn].x() << "," << exr[ptn].y() << ")    "
             << "   prev "  << ptp << ", (x,y)=(" << exr[ptp].x() << "," << exr[ptp].y() << ")    "
             << "  (dxp,dyp)=" << dxp << "," << dyp
             << "  (dxn,dyn)=" << dxn << "," << dyn << endl;

        if (dxp * dyn < dyp * dxn) {
          bgeo::append(mpt, exr[pti]);
          cout << "     >>>> CONCAVE" << endl;
        }

        dxp = -dxn; dyp = -dyn;

      }
    }

    int ncpts = mpt.size();
    bg_mlstr mline; int nmline = 0;
    for (int pti = 0; pti < ncpts;  pti++) {
      for (int ptj = 0; ptj < ncpts; ptj++) {
        bg_lstr line;
        line.push_back(mpt[pti]);
        line.push_back(mpt[ptj]);

        if (bgeo::within(line, mpoly)) {
          cout << "line length=" << bgeo::length(line) << "   within:" << 1 << endl;
          cout << "   line: " << bgeo::dsv(line) << endl;
          mline.resize(++nmline);
          bgeo::append(mline[nmline-1], line);
        }
      }
    }
    
    cout << "mline: " << bgeo::dsv(mline) << endl;


    // const int points_per_circle = 20;
    // boost::geometry::strategy::buffer::distance_symmetric<double> distance_strategy(max_dist);
    // boost::geometry::strategy::buffer::join_round join_strategy(points_per_circle);
    // boost::geometry::strategy::buffer::end_round end_strategy(points_per_circle);
    // boost::geometry::strategy::buffer::point_circle circle_strategy(points_per_circle);
    // boost::geometry::strategy::buffer::side_straight side_strategy;

    // Create the buffer of a multi point
    // boost::geometry::model::multi_polygon<bg_poly> vd_buff;
    // boost::geometry::buffer(vd_ctr, vd_buff, distance_strategy, side_strategy,
    //                                          join_strategy, end_strategy, circle_strategy);

    // Declare a stream and an SVG mapper
    std::ofstream svg("my_map.svg");
    boost::geometry::svg_mapper<bp_pt> mapper(svg, 400, 400, "width=\"120%\" height=\"120%\"");

    // Add geometries such that all these geometries fit on the map
    mapper.add(mpoly);
    mapper.add(mpt);
    mapper.add(mline);

    // Draw the geometries on the SVG map, using a specific SVG style
    // Chicago grey=(128,122,117) and red=(139,0,33)
    mapper.map(mpoly, "fill-opacity:0.8;fill:rgb(128,122,117);");
    mapper.map(mpt,   "fill-opacity:1.0;fill:rgb(0,0,255);stroke:rgb(0,0,255);stroke-width:0,2", 3);
    mapper.map(mline, "fill-opacity:1.0;fill:rgb(0,0,255);stroke:rgb(0,0,255);stroke-width:0,2", 1);

    return 0;

}

