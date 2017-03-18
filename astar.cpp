// g++ sip.cxx -I ~/anaconda3/include/ -std=c++11 -o sip

#include "shortest_path.hpp"
#include "astar_test.h"

typedef boost::polygon::point_data<double> bp_pt;
typedef boost::geometry::model::polygon<bp_pt> bg_poly;

int main(int argc, char **argv) {

  bg_poly poly;
  boost::geometry::read_wkt(pa7, poly);
  cl_shortest_path::shortest_path(poly);

  return 0;
  
}
