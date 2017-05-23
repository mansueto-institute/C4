#pragma once

// Adapted to boost from --
// https://github.com/mapbox/polylabel

#include <algorithm>
#include <cmath>
#include <iostream>
#include <queue>

#include <iostream>
#include <fstream>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/linestring.hpp>

#include <boost/polygon/voronoi.hpp>
#include <boost/geometry/geometries/adapted/boost_polygon.hpp>


namespace ehrenburg {

  using std::cout;
  using std::endl;
  namespace bgeo  = boost::geometry;
  namespace bpoly = boost::polygon;
  
  typedef bpoly::point_data<double> bp_pt;
  typedef bpoly::segment_data<double> bp_seg;
  
  namespace bgmod   = boost::geometry::model;
  namespace bgtrans = boost::geometry::strategy::transform;


  const double sqrt2 = 1.414213562373095;

  // signed distance from point to polygon outline (negative if point is outside)
  template <class T, class pT>
  T pointToPolygonDist(const pT& pt, 
                       const bgmod::polygon<pT>& poly,
                       const bgmod::linestring<pT>& line) {

    bool inside  = false;

    if (bgeo::within(pt, poly)) { inside = true; }

    // would be nice to switch to comparable distance, 
    // but we add on the h * sqrt2, to get the max distance,
    // which anyway requires the square root...
    T d = bgeo::distance(pt, line); 

    return (inside ? 1 : -1) * d;

  }

  template <class T, class pT>
  struct Cell {
    Cell(const pT& c_, T h_, 
         const bgmod::polygon<pT>& polygon,
         const bgmod::linestring<pT>& line)
      : c(c_),
        h(h_),
        d(pointToPolygonDist<T, pT>(c, polygon, line)),
        max(d + h * sqrt2)
    {}

    pT c; // cell center
    T h; // half the cell size
    T d; // distance from cell center to polygon
    T max; // max distance to polygon within a cell
  };

  template <class T, class pT>
  std::pair<pT, T>
  pole_inacc(const pT& guess_point,
             const bgmod::polygon<pT >& polygon, 
             double precision = 1, bool debug = false) {

    // find the bounding box of the outer ring
    bgmod::box<pT> env;
    bgeo::envelope(polygon, env);

    const pT size { (env.max_corner().x() - env.min_corner().x()) / 2.0,
                    (env.max_corner().y() - env.min_corner().y()) / 2.0 };

    const T cellSize = std::min(size.x(), size.y());
    T h = cellSize;

    // a line string for calculating distances.
    bgmod::linestring<pT> line;
    for (auto ipt : bgeo::exterior_ring(polygon)) line.push_back(ipt);


    // Perhaps, since this is iterative, we can benefit 
    // from starting at the last location....
    auto bestCell = Cell<T, pT>(guess_point, 0, polygon, line);


    // a priority queue of cells in order of their "potential" (max distance to polygon)
    auto compareMax = [] (const Cell<T, pT>& a, const Cell<T, pT>& b) { return a.max < b.max; };
    std::priority_queue<Cell<T, pT>, std::vector<Cell<T, pT>>, decltype(compareMax)> cellQueue(compareMax);

    if (cellSize == 0) return std::make_pair(env.min_corner(), 0);

    // cover polygon with initial cells
    for (T x = env.min_corner().x(); x < env.max_corner().x(); x += cellSize) {
      for (T y = env.min_corner().y(); y < env.max_corner().y(); y += cellSize) {
        cellQueue.push(Cell<T, pT>({x + h, y + h}, h, polygon, line));
      }
    }

    // special case for rectangular polygons
    pT cell_pt;
    bgtrans::translate_transformer<T, 2, 2> translate(size.x(), size.y());
    bgeo::transform(env.min_corner(), cell_pt, translate);
    Cell<T, pT> bboxCell(cell_pt, 0, polygon, line);
    if (bboxCell.d > bestCell.d) bestCell = bboxCell;

    auto numProbes = cellQueue.size();
    while (!cellQueue.empty()) {
      // pick the most promising cell from the queue
      auto cell = cellQueue.top();
      cellQueue.pop();

      // update the best cell if we found a better one
      if (cell.d > bestCell.d) {
        bestCell = cell;
        if (debug) std::cout << "found best " << std::round(1e4 * cell.d) / 1e4 << " after " << numProbes << " probes" << std::endl;
      }

      // do not drill down further if there's no chance of a better solution
      if (cell.max - bestCell.d <= precision) continue;

      // split the cell into four cells
      h = cell.h / 2;
      cellQueue.push(Cell<T, pT>({cell.c.x() - h, cell.c.y() - h}, h, polygon, line));
      cellQueue.push(Cell<T, pT>({cell.c.x() + h, cell.c.y() - h}, h, polygon, line));
      cellQueue.push(Cell<T, pT>({cell.c.x() - h, cell.c.y() + h}, h, polygon, line));
      cellQueue.push(Cell<T, pT>({cell.c.x() + h, cell.c.y() + h}, h, polygon, line));
      numProbes += 4;
    }

    if (debug) {
      std::cout << "num probes: " << numProbes << std::endl;
      std::cout << "best distance: " << bestCell.d << std::endl;
    }

    return std::make_pair(bestCell.c, bestCell.d);
  }


  template <class T, class pT>
  std::pair<pT, T>
  pole_inacc(const bgmod::polygon<pT>& polygon, 
             double precision = 1, bool debug = false) {
  
    // take centroid as the first best guess
    pT ctrd;
    bgeo::centroid(polygon, ctrd);

    return pole_inacc<T, pT>(ctrd, polygon, precision, debug);

  }

}


