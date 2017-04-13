#include "Cluscious.h" 


namespace Cluscious {

  static std::map<ObjectiveMethod, RadiusType> obj_radius = {{DISTANCE_A, EQUAL_AREA}, {DISTANCE_P, EQUAL_AREA_POP},
                                                             {INERTIA_A, EQUAL_AREA}, {INERTIA_P, EQUAL_AREA_POP}, 
                                                             {HULL_A, HULL}, {POLSBY, EQUAL_CIRCUMFERENCE},{REOCK, SCC}, 
                                                             {EHRENBURG, LIC}, {POLSBY_W, EQUAL_CIRCUMFERENCE}, 
                                                             {PATH_FRAC, EQUAL_AREA}};

  Cell::Cell()
    : region(-1), id(0), pop(0), x(0), y(0), area(0) {}
  Cell::Cell(const Cell &c) 
    : region(-1), id(c.id), pop(c.pop), x(c.x), y(c.y), area(c.area), wm(c.wm),
      is_univ_edge(c.is_univ_edge), is_split(c.is_split), split_neighbor(c.split_neighbor),
      pt(c.pt) {} // , poly(c.poly){}

  Cell::Cell(int i, int p, double x, double y,
             double a, weight_map wm,
             bool ue = false, bool split = false)
    : region(-1), id(i), pop(p), x(x), y(y), area(a), wm(wm), is_univ_edge(ue), is_split(split), split_neighbor(false), pt(x, y) {

    if (!p) pop = 1; // N.B. that we're setting the minimum population per cell to 1!!

    // boost::geometry::read_wkt(mp_wkt, poly);  
  }


  void Cell::adjacency_to_pointers(std::vector<Cell*>& univ_cells) {

    weight_map::iterator wit;
    weight_map::iterator wmend = wm.end();

    // iterate over the full universe
    for (auto c : univ_cells) {

      // if this cell's id is in our weight map,
      // add the cell and its weight to the 
      // neighbor map.
      wit = wm.find(c->id);
      if (wit != wmend) nm.push_back(std::make_pair(c, wit->second));

    }
  }

  void Cell::add_edge(int edge_id, int nodea, int nodeb) {

    edges.push_back(Edge(edge_id, nodea, nodeb));

  }


  void Cell::node_ids_to_pointers(std::vector<Node*>& univ_nodes) {

    for (auto& n : univ_nodes) {
      for (auto& e : edges) {
        if (e.na_id == n->id) { e.set_na(n); nodes.insert(n); }
        if (e.nb_id == n->id) { e.set_nb(n); nodes.insert(n); }
      }
    }

    if (is_split) {
      for (auto& n : nm) {
        n.first->split_neighbor = true;
      }
    }
  }

  int Cell::get_edge_idx(int id) {
    
    int eidx = -1;
    for (auto& e : edges) { eidx++; 
      if (e.id == id) return eidx;
    }

    return eidx; // shit got fucked.

  }

  int Cell::get_ext_edge_idx() {

    int eidx = -1;

    // If this cell is not on the border of the universe,
    // we're looking for an edge shared with a foreign neighbor.
    if (!is_univ_edge) {

      for (auto e : edges) { eidx++;                 // For each edge,
        for (auto nei : nm) {                        // loop over the 
          if (nei.first->region == region) continue; // foreign neighbors'
          for (auto nei_e : nei.first->edges) {      // edges;
            if (nei_e.id == e.id) {                  // if the edge is in both,
              return eidx;                         // return the edge index. 
            }
          }
        }
      }

    } else { // otherwise, looking for an edge shared with no one.

      for (auto e : edges) { eidx++;                 // For each edge,
        bool found_edge = false;
        for (auto nei : nm) {                        // loop over the neighbors'
          for (auto nei_e : nei.first->edges) {      // edges;
            if (nei_e.id == e.id) {                  // if it's shared, 
              found_edge = true;                     // it's not what we're looking for.
              break;
            }
          }
          if (found_edge) break;
        }
        if (found_edge) continue;

        return eidx;
      }
    }

    return -1; // otherwise -1, to signify failure.

  }

  std::pair<Node*, Node*> Cell::get_cw_node(int reg) {

    std::pair<Node*, Node*> rval(0, 0);
    cout << "On get_cw_node for cell " << id << " (x,y)=(" << x << "," << y << ") in region=" << reg << ")" << endl;

    // Do the last one first, to deal with 
    // the first transition.
    bool na_dom_nei, nb_dom_nei;
    for (auto& e : edges) {
      na_dom_nei = nb_dom_nei = false;
      for (auto na_e : e.na->edges) {
        for (auto n : nm) {
          if (n.first->region == reg &&
              n.first->has_edge(na_e)) {

            na_dom_nei = true;
            break;
          }
        }
        if (na_dom_nei) break;
      }

      for (auto nb_e : e.nb->edges) {
        for (auto n : nm) {
          if (n.first->region == reg &&
              n.first->has_edge(nb_e)) {
            nb_dom_nei = true;
            break;
          }
        }
        if (nb_dom_nei) break;
      }

      // Is THIS edge in the region?
      // Must do false -> true, because univ. edges are not in the region
      // but won't give any contradiction -- there is no cell OUT of region on the edge.
      bool edge_in_region = false;
      for (auto n : nm) {
        if (n.first->region == reg &&
            n.first->has_edge(e.id)) {
          edge_in_region = true;
        }
      }

      cout << "  Edge id=" << e.id << "  inreg=" << edge_in_region << "  na_dom_nei=" << na_dom_nei << "  nb_dom_nei=" << nb_dom_nei << endl;
      if ( na_dom_nei && !nb_dom_nei) rval.first  = e.na;
      if (!na_dom_nei &&  nb_dom_nei) rval.second = e.nb;

      if (!edge_in_region && na_dom_nei && nb_dom_nei)
        rval.first = e.na; rval.second = e.nb;

      // if (rval.first && rval.second) break;
    }

    // if (!rval.first || !rval.second) 
    //   throw std::runtime_error(std::string("get_cw_node() got an empty start or end!"));

    return rval;

  }

  bool Cell::has_edge(int edge_id) {

    for (auto e : edges) {
      if (edge_id == e.id) return true;
    }

    return false;
  }

  Edge* Cell::get_edge(int edge_id) {

    for (auto e = edges.begin(); e != edges.end(); e++) {
      if (edge_id == e->id) return &(*e);
    }

    return 0;
  }

  bool Cell::next_edge_in_region(Node* node, int start_edge_id,
                                 Cell*& next_cell, Edge*& next_edge,
                                 bool CW = true, int region_id = -1) {

    if (region_id < 0) region_id = region;

    int edge_id;
    node->set_edge(start_edge_id);
    // cout << "cell " << id << "(R=" << region << "), " << (CW ? "adding to" : "removing from") << " region_id=" << region_id << endl;
    // cout << "We're starting at edge " << start_edge_id << " :: follows :: ";
    for (int ni = 0; ni < node->size(); ni++) {

      if (CW) edge_id = node->next();
      else    edge_id = node->prev();
      // cout << edge_id << " ";

      int border = (region == region_id) && has_edge(edge_id);
      if (border) {
        next_cell = this;
        next_edge = get_edge(edge_id);
        // cout << "C";
      }

      for (auto& nei : nm) {
        if (nei.first->region == region_id && 
            nei.first->has_edge(edge_id)) {

          next_cell = nei.first;
          next_edge = nei.first->get_edge(edge_id);
          border++;
          // cout << "N";
        }
      }
      // cout << "/ ";

      if (border == 1) {
        // cout << endl;
        return true;
      }
    }
    // cout << endl;

    // Should never arrive here.
    assert(0 && "Never found the next border -- made a full loop!!\nIs the cell's region set correctly?");

    return false;

  }

  std::unordered_set<int> Cell::neighbor_regions(bool QUEEN) {

    std::unordered_set<int> rval;
    for (auto& n : nm)
      if (n.first->region != region && 
          (QUEEN || n.second > 0))
        rval.insert(n.first->region);

    return rval;

  }


  // Overloading...
  bool Cell::neighbors_connected(bool QUEEN = false, bool FAST = false) {

    std::vector<std::unordered_set<Cell*> > graphs;
    return (neighbor_sets(graphs, false, QUEEN, FAST) == 1);

  }

  int Cell::neighbor_strands(std::unordered_set<Cell*>& strand, // int dest_reg,
                             unsigned int max_merge = 0, bool QUEEN = false) {

    std::vector<std::unordered_set<Cell*> > graphs;
    int nsets = neighbor_sets(graphs, max_merge, false);

    if (nsets == 1) return 1;
    if (nsets >  2) return nsets;

    for (auto& g : graphs) {
      if (g.size() < max_merge) {
        strand = g;
        return 2;
      }
    }

    return nsets;
  }

  // This is quite similar to "connect graph."
  int Cell::neighbor_sets(std::vector<std::unordered_set<Cell*> >& graphs,
                          unsigned int max_merge, bool QUEEN, bool FAST) {

    for (auto& n : nm) {

      if (n.first->region != region) continue;

      // First make the local star graph, 
      // for now, just with cells connected to the original,
      // and _in the same region._
      std::unordered_set<Cell*> lone_star;
      lone_star.insert(n.first);
      for (auto& nn : n.first->nm) 
        if (nn.first->region == region &&
            // best this way?  removing this would allow corners...
            nn.first->wm.find(id) != nn.first->wm.end() &&
            (QUEEN || nn.second > 0)) // ROOK
          lone_star.insert(nn.first);

      // Then get a vector of all the existing
      // graphs that it's connected to:
      // if any element in lone star is in an 
      // existing graph, just mark it.
      std::vector<int> gidx;
      for (size_t gi = 0; gi < graphs.size(); gi++) {
        for (auto& c : lone_star) {
          if (graphs[gi].find(c) != graphs[gi].end()) {
            gidx.push_back(gi);
            break;
          }
        }
      } // gid should now hold a list of overlapping graphs

      // If it has (yet) no neighbors, add it on its own.
      if (!gidx.size()) { 
        graphs.push_back(lone_star);

      // Otherwise, add it to the first,
      // and add _others_ on to that one.
      } else {
        graphs[gidx[0]].insert(lone_star.begin(), lone_star.end());

        // backwards, to keep removals consistent.
        for (auto gi : boost::adaptors::reverse(gidx)) {
          if (gi == gidx[0]) continue;
          graphs[gidx[0]].insert(graphs[gi].begin(), graphs[gi].end());
          graphs.erase(graphs.begin() + gi);
        }
      }
    }

    // If the cell's neighbors are connected, we're done.
    if (graphs.size() <= 1) return 1;
    else if (FAST) return 2;

    // Otherwise, expand search.

    // Create a set of cells "to add", for each set.
    std::vector<std::unordered_set<Cell*> > to_add;

    // Add the neighbors of the current sets to each list.
    for (auto g : graphs) {
      to_add.push_back(std::unordered_set<Cell*>());
      for (auto c : g) for (auto n : c->nm) if (QUEEN || n.second > 0) { // ROOK
        if (n.first->region == region &&   // in the region
            g.find(n.first) == g.end() &&  // not yet in graph, and 
            n.first != this) {             // not the cell in question.
          to_add.back().insert(n.first);
        }
      }
    }

    // cout << " >> We've now expanded the search: c" << id << " >> " << endl;

    // for (int gidx = 0; gidx < int(graphs.size()); gidx++) {
    //   cout << "    Graph contains [[ ";
    //   for (auto c : graphs[gidx]) cout << c->id << " ";
    //   cout << "]]  +  [[ ";
    //   for (auto c : to_add[gidx]) cout << c->id << " ";
    //   cout << "]]" << endl;
    // }

    std::unordered_set<Cell*> next;

    // Now go about building these sets up.
    // If they are merged or we can't go further, we're done.
    while (graphs.size() > 1 && 

              // may not be empty, but can still flesh out the sets.
           (( max_merge && !maxed_or_empty(max_merge, graphs, to_add)) ||   
           // (( max_merge && std::any_of(to_add.begin(), to_add.end(), uo_set_nempty)) ||  

              // all of them have neighbors, so we've gotta keep searching.
            (!max_merge && std::all_of(to_add.begin(), to_add.end(), uo_set_nempty)))) {


      for (int gidx = 0; gidx < int(graphs.size()); gidx++) {

        for (int giidx = int(graphs.size())-1; giidx > gidx; giidx--) {

          for (auto c : to_add[giidx]) {
            if (to_add[gidx].find(c) != to_add[gidx].end() ||
                graphs[gidx].find(c) != graphs[gidx].end()) {

              if (graphs.size() == 2) {
                // cout << "    cid=" << c->id << "  ... could return at this point..." << endl << endl;
                return 1;
              }

              to_add[gidx].insert(to_add[giidx].begin(), to_add[giidx].end());
              graphs[gidx].insert(graphs[giidx].begin(), to_add[giidx].end());
              
              to_add.erase(to_add.begin() + giidx);
              graphs.erase(graphs.begin() + giidx);

              giidx--;
              break;

            }
          }
        }

      }

      for (int gidx = 0; gidx < int(graphs.size()); gidx++) {

        graphs[gidx].insert(to_add[gidx].begin(), to_add[gidx].end());

        next.clear();
        for (auto c : to_add[gidx]) {
          for (auto n : c->nm) if (QUEEN || n.second > 0) { // ROOK
            // cout << "nid=" << n.first->id << "  r" << n.first->region << " v. " << region << " :: ";
            if (n.first->region == region && c != this &&
                graphs[gidx].find(n.first) == graphs[gidx].end()) {
              // cout << " (added)";
              next.insert(n.first);
            }
            // cout << endl;
          }
        }
        to_add[gidx].clear();
        to_add[gidx].insert(next.begin(), next.end());

      }

      // cout << "    >>>>>>>>>>>>>>>>    " << endl;
      // for (int gidx = 0; gidx < int(graphs.size()); gidx++) {
      //   cout << "    Graph contains [[ ";
      //   for (auto c : graphs[gidx]) cout << c->id << " ";
      //   cout << "]]  +  [[ ";
      //   for (auto c : to_add[gidx]) cout << c->id << " ";
      //   cout << "]]" << endl;
      // }

    }

    // cout << "    *****************" << endl << "    ***  RETURNING " << graphs.size() << endl << endl;

    // cout << endl;

    return graphs.size();

  }

  void Cell::merge(Cell* m) {

    // Do we want to move the coords?  I think not.
    // x = (x * area + m->x * m->area)/(area + m->area);
    // y = (y * area + m->y * m->area)/(area + m->area);

    area += m->area;
    pop  += m->pop;

    wm.erase(m->id);

    int nm_idx = 0;
    for (auto nei : nm) {
      if (nei.first == m) {
        nm.erase(nm.begin()+nm_idx);
        break;
      }
      nm_idx++;
    }

    enclaves_and_islands.push_back(m->id);

  }

  Region::Region() 
    : id(0), pop(0), ncells(0), area(0), xctr(0), yctr(0), xpctr(0), ypctr(0), sumw_border(0),
      has_topo(false), eps_scc(10), x_scc(0), y_scc(0), r2_scc(0),
      eps_lic(2), x_lic(0), y_lic(0), r2_lic(std::numeric_limits<double>::infinity()) {}

  Region::Region(int rid) 
    : id(rid), pop(0), ncells(0), area(0), xctr(0), yctr(0), xpctr(0), ypctr(0), sumw_border(0),
      has_topo(false), eps_scc(10), x_scc(0), y_scc(0), r2_scc(0), 
      eps_lic(2), x_lic(0), y_lic(0), r2_lic(std::numeric_limits<double>::infinity()) {}

  Region::Region(int rid, Cell* c) 
    : id(rid), pop(0), ncells(0), area(0), xctr(c->x), yctr(c->y), xpctr(c->x), ypctr(c->y), sumw_border(0),
      has_topo(false), eps_scc(10), x_scc(0), y_scc(0), r2_scc(0),
      eps_lic(2), x_lic(0), y_lic(0), r2_lic(std::numeric_limits<double>::infinity())
  { add_cell(c, true); }


  void Region::add_cell(Cell* c, bool UPDATE_CTR = true) {

    ncells++;
    c->region = id;
    cells.insert(c);

    add_cell_int_ext_neighbors(c);

    // if (has_topo && node_ring.size()) divert_ring_at_cell(c, true);
    if (has_topo && node_ring.size()) make_ring();

    if (UPDATE_CTR) {
      xctr  = (xctr  * area + c->x * c->area)/(area + c->area);
      yctr  = (yctr  * area + c->y * c->area)/(area + c->area);
      xpctr = (xpctr * pop  + c->x * c->pop )/(pop  + c->pop );
      ypctr = (ypctr * pop  + c->y * c->pop )/(pop  + c->pop );

      update_scc(c, 0, true); // add cell, do update.
      if (has_topo && node_ring.size()) update_lic(c, 0, true);
    }

    area += c->area;
    pop  += c->pop;

    // bgeo::union_(poly, c->poly, poly);
    bgeo::append(mpt, c->pt);

    if (!bgeo::within(c->pt, ch_poly)) {
      bgeo::clear(ch_poly);
      bgeo::convex_hull(mpt, ch_poly);
      ch_area = bgeo::area(ch_poly);
    }

  }

  void Region::add_cell_int_ext_neighbors(Cell *c) {

    ext_borders.erase(c);

    // require check, because otherwise cells can
    // get assigned after all ext. neighbors are,
    // and then never re-visited.
    if (c->is_univ_edge) int_borders.insert(c);
    for (auto n : c->nm) if (n.first->region != id) {
      int_borders.insert(c); break;
    }

    // Iterate over the MOVING cell's neighbors.
    for (auto const& n : c->nm) {

      // If the neighbor is not a member of the
      // destination region, then this is a new border.
      // We try to add it to the list of border cells
      // (it may already be there), but at any rate
      // add its newly-exposed border to sumw_border.
      if (n.first->region != id) {
        if (n.second > 0) ext_borders.insert(n.first); // ROOK
        sumw_border += n.second;

      // If the cell's neighbor IS a member of
      // the destination region, then this border was 
      // previously exposed.  Remove its contribution 
      // to sumw_border.
      } else {
        sumw_border -= n.second;

        // if, further, that neighbor has no other connections
        // to the outside world, then remove it from int_borders.
        bool indep_connections = false;
        for (auto nn : n.first->nm) {

          if (nn.first->region != id) {
            indep_connections = true;
            break;
          }
        }
        if (!indep_connections && !n.first->is_univ_edge) {
          int_borders.erase(n.first);
        }
      }
    }
    
    // Sanity check.
    // for (auto b : ext_borders) {
    //   float total_border = 0;
    //   for (auto bn : b->nm) {
    //     if (bn.first->region == id)
    //       total_border += bn.second;
    //   }
    //   if (total_border == 0) {
    //     DEBUG_ME;
    //     cout << "cell " << b->id << " is in the external border of Region " << id << " but has no non-zero (ROOK) connection to it." << endl;
    //     exit(1);
    //   }
    // }

  }

  void Region::remove_cell(Cell *c, bool UPDATE_CTR = true) {

    if (ncells == 1) return; // FIXME how to treat this.

    if (c->region == id) c->region = -1;

    area -= c->area;
    pop  -= c->pop;

    ncells--;
    cells.erase(c);

    remove_cell_int_ext_neighbors(c);

    // if (has_topo && node_ring.size()) divert_ring_at_cell(c, false);
    if (has_topo && node_ring.size()) make_ring();

    if (UPDATE_CTR) {
      xctr  = (xctr  * area - c->x * c->area)/(area - c->area);
      yctr  = (yctr  * area - c->y * c->area)/(area - c->area);
      xpctr = (xpctr * pop  - c->x * c->pop )/(pop  - c->pop );
      ypctr = (ypctr * pop  - c->y * c->pop )/(pop  - c->pop );

      update_scc(0, c, true); // sub cell, do update.
      if (has_topo && node_ring.size()) update_lic(0, c, true);
    }

    mpt.clear(); // faster to just re-copy.
    for (auto c : cells) bgeo::append(mpt, c->pt);

    // check "within" because corners count as outside.
    if (!bgeo::within(c->pt, ch_poly)) {
      bgeo::clear(ch_poly);
      bgeo::convex_hull(mpt, ch_poly);
      ch_area = bgeo::area(ch_poly);

      bgeo::centroid(ch_poly, ch_pt);
      ch_x = bgeo::get<0>(ch_pt);
      ch_y = bgeo::get<1>(ch_pt);
    }

  }

  void Region::remove_cell_int_ext_neighbors(Cell* c) {

    int_borders.erase(c);
    ext_borders.insert(c);

    // Iterate over the MOVING cell's neighbors.
    for (auto const& n : c->nm) {

      // If the neighbor cell is not a member of the 
      // source (= id, this) region, remove this boundary's 
      // contribution to sumw_border.  If, furthermore,
      // that cell does not have any *other* connections
      // to the source, erase it from the border.
      if (n.first->region != id) {
        sumw_border -= n.second;

        bool indep_connections = false;
        for (auto nn : n.first->nm) {

          if (nn.first->region == id && 
              nn.second > 0) { // ROOK
            indep_connections = true;
            break;
          }
        }
        if (!indep_connections) {
          ext_borders.erase(n.first);
        } 

      // If the moving cell's neighbor IS a member of
      // source region, add its contribution to 
      // the sumw_border.  (It is an internal border.)
      } else {

        int_borders.insert(n.first);

        sumw_border += n.second;
      }
    }

    // Sanity check.
    for (auto b : ext_borders) {
      float total_border = 0;
      for (auto bn : b->nm) {
        if (bn.first->region == id)
          total_border += bn.second;
      }
      if (total_border == 0) {
        DEBUG_ME;
        cout << "cell " << b->id << " is in the external border of Region " << id << " but has no non-zero (ROOK) connection to it." << endl;
      }
    }

  }

  void Region::make_ring() { 

    node_ring.clear();
    get_node_ring(node_ring);
  }

  std::pair<std::pair<float, float>, float> Region::get_circle_coords(RadiusType rt) {

    switch (rt) {

      case RadiusType::EQUAL_AREA_POP:       return std::make_pair(std::make_pair(xpctr, ypctr), sqrt(area/M_PI));   
      case RadiusType::EQUAL_CIRCUMFERENCE:  return std::make_pair(std::make_pair(xctr,  yctr),  sumw_border/2/M_PI );  
      case RadiusType::SCC:                  return std::make_pair(std::make_pair(x_scc, y_scc), sqrt(r2_scc)); 
      case RadiusType::LIC:                  return std::make_pair(std::make_pair(x_lic, y_lic), sqrt(r2_lic));  
      case RadiusType::HULL:                 return std::make_pair(std::make_pair(ch_x, ch_y)  , sqrt(ch_area/M_PI));
      case RadiusType::EQUAL_AREA:
      default:                               return std::make_pair(std::make_pair(xctr,  yctr),  sqrt(area/M_PI));   
    }

  }

  std::vector<std::pair<float, float> > Region::get_point_ring() {

    if (!node_ring.size()) make_ring();
    // make_ring();

    std::vector<std::pair<float, float> > point_ring;

    for (auto n : node_ring) {
      point_ring.push_back(std::make_pair(n->x, n->y));
    }
    if (node_ring.size()) point_ring.push_back(point_ring.front()); // close the loop.
    else {
      cerr << "Found a ring with size 0 in region " << id << endl;
      throw std::runtime_error(std::string("Found a ring with size 0"));
    }

    return point_ring;

  }

  // Not directly returning the float pairs makes it easier
  // to add in the queen/corner loops.
  void Region::get_node_ring(std::vector<Node*> &ring,
                             Cell* curr_cell, int this_edge_idx) {

    bool recursed = curr_cell; // safety -- one layer deep.

    // If a starting cell/edge is not provided, find a random one.
    // if we happen to choose a cell with only a corner connection
    // to the outside, the external ege will fail; hence we loop.
    for (auto ib = int_borders.begin();
         ib != int_borders.end() && this_edge_idx < 0; ++ib) {
      curr_cell = *ib;
      if (curr_cell->region != id) continue;
      if (curr_cell->is_split) continue; // just not worth starting here.
      if (curr_cell->split_neighbor) continue;
      this_edge_idx = curr_cell->get_ext_edge_idx();
    }
    if (this_edge_idx < 0) {
      throw std::runtime_error(std::string("No good starting edge cell found!!"));
    }

    Edge* this_edge = &curr_cell->edges[this_edge_idx];
    // cout << "Starting cell for node ring " << curr_cell->id << endl;
    Node* node = this_edge->nb;

    // cout << "Region " << id << ", cell " << curr_cell->id << " at node " << node->id << endl; 

    Node* start_node = this_edge->na;
    if (start_node->size() > 2) ring.push_back(start_node);

    // This will store any corner cases we miss.
    std::map<Node*, std::pair<Cell*, int> > diagonals;

    // Looping around nodes of of the polygon.
    int loop_safe = 0; // because I lack confidence, for the time being....
    // while ((node != start_node || ring.size() <= 1) && loop_safe < 1000) {
    while (node != start_node && loop_safe < 1000) {

      loop_safe++; 

      // Loop around CW, until the next edge on the border.
      Cell* next_cell = 0; Edge* next_edge = 0;
      curr_cell->next_edge_in_region(node, this_edge->id, next_cell, next_edge);

      // If it's a new cell, we mark the node in our ring.
      if (node->size() > 2) ring.push_back(node);
      // if (next_cell != curr_cell) ring.push_back(node);

      // If there's the possibility of catty-corner cells, look in reverse too.
      if (node->size() > 3) {

        Cell* corner_cell; Edge* corner_edge;
        curr_cell->next_edge_in_region(node, this_edge->id, corner_cell, corner_edge, false); // (false -> CCW)

        // If the two directions give different edges, we've got a corner.
        if (corner_edge != next_edge) {

          // If it's the first time we've seen this cell, 
          if (diagonals.find(node) == diagonals.end()) {
            // add it to the list of potential diagonals.
            diagonals[node] = std::make_pair(corner_cell, corner_cell->get_edge_idx(corner_edge->id));
          } else {
            diagonals.erase(node); // otherwise, remove it from that list. 
          }
        }
      }

      node = next_edge->nb;  // that will be our next node
      this_edge = next_edge; // edge, and
      curr_cell = next_cell; // cell...

    }

    if (recursed) return;
      
    while (diagonals.size()) {
      // cout << "Diagonal nodes :: " << diagonals.size() << endl;

      auto nd = diagonals.begin();
      // cout << "  > " << nd->first->id << " (" << nd->second.first << ", " << nd->second.second << ")" << endl;

      std::vector<Node*> insert_ring;
      get_node_ring(insert_ring, nd->second.first, nd->second.second);

      std::vector<Node*>::iterator ring_iter = find(ring.begin(), ring.end(), nd->first);
      ring.insert(ring_iter, insert_ring.begin(), insert_ring.end());

      for (auto ni : insert_ring) {
        diagonals.erase(ni);
      }

      diagonals.erase(nd->first);
    }

    if (loop_safe == 1000) {
      cerr << "WARNING :: Region " << id << " :: get_node_ring reached max loop value. "
           << "Output size is " << ring.size() << "." << endl;
    }

    if (!ring.size()) {

      assert(ring.size());
    }

  }

  void Region::divert_ring_at_cell(Cell* c, bool CW = true) {

    DEBUG_ME; cout << "\n >>>>>>>>>>>>>>>>> \nREGION=" << id << ", " << (CW ? "adding" : "removing") << " cell" << c->id << endl;

    int number_on_cell = 0;
    int insertion_point = 0;
    Node *snip_a = 0, *snip_b = 0;

    if (!c->is_split && !c->split_neighbor) {

      int ni = 0; int node_ring_size = node_ring.size();
      if (c->nodes.find(node_ring.back()) != c->nodes.end()) {
        ni = node_ring_size - c->nodes.size() - 1;
      }

      // cout << "I'm on cell " << c->id << ", which is not split." << endl;

      // not wasteful, since the second expression evaluates, 
      int loop_safety = 0;
      while ((!snip_b || c->nodes.find(node_ring[ni]) != c->nodes.end()) && loop_safety < 2) {

        // If snip_b is defined, we already tested, so don't retest...
        if (snip_b || c->nodes.find(node_ring[ni]) != c->nodes.end()) { 

          // Otherwise, this node IS on the cell.
          number_on_cell++;
          cout << "node=" << node_ring[ni]->id << endl;
          if (!snip_a) {
            snip_a = node_ring[ni];
            insertion_point = ni;
          }

          snip_b = node_ring[ni];
        }

        ni = (ni+1) % node_ring_size;
        if (!ni) loop_safety++;
      }


    } else {

      // cout << "We've got cell " << c->id << ", which is unconnected!" << endl;

      bool last, curr;
      int nr_size = node_ring.size();

      last = (c->nodes.find(node_ring.back()) != c->nodes.end());
      // cout << "last at fwd start=" << last << endl;
      for (int ni = 0; ni < nr_size; ni++) {
        curr = (c->nodes.find(node_ring[ni]) != c->nodes.end());
        if (curr && !last) {
          snip_a = node_ring[ni];
          insertion_point = ni;
          break;
        } 
        last = curr;
      }
      
      last = (c->nodes.find(node_ring.front()) != c->nodes.end());
      // cout << "last at bkwd start=" << last << endl;
      for (int ni = nr_size-1; ni >= 0; ni--) {
        curr = (c->nodes.find(node_ring[ni]) != c->nodes.end());
        if (curr && !last) {
          snip_b = node_ring[ni];
          number_on_cell = ni - insertion_point + 1;
          if (number_on_cell < 0) number_on_cell += nr_size;
          break;
        } 
        last = curr;
      }

      // cout << "region=" << id << "  cell id=" << c->id << "  number=" << number_on_cell << "  insertion_point=" << insertion_point << "  snip_a=" << snip_a->id << "  snip_b=" << snip_b->id << endl;
    }

    assert(snip_a && snip_b);


    Edge* start_edge = &c->edges[0];
    for (auto& e : c->edges) {
      if (( CW && e.na == snip_a) ||
          (!CW && e.nb == snip_a)) {
        start_edge = &e;
        break;
      }
    }

    Cell* this_cell = c;
    Edge* this_edge = start_edge;
    Node* node = CW ? start_edge->nb : start_edge->na; // if removing do node a.

    std::vector<Node*> cell_path;

    // Looping around nodes of of the polygon.
    size_t loop_safe = 0; // because I lack confidence, for the time being....
    Cell* next_cell = 0; Edge* next_edge = 0;
    cout << "cell id=" << c->id << "  snipA=" << snip_a->id << "  snipB=" << snip_b->id
         << "  number_on_cell=" << number_on_cell
         << "  edge=" << this_edge->id 
         << "  (a=" << this_edge->na_id << ", b=" << this_edge->nb_id << ")" << endl;
    while (node != snip_b && loop_safe < c->nodes.size()) {

      cout << "this_cell=" << this_cell->id 
           << "  this node=" << node->id << "(x,y)=(" << node->x << "," << node->y << ")"
           << "  this_edge=" << this_edge->id 
           << "  (a=" << this_edge->na_id << ", b=" << this_edge->nb_id << ")" << endl;

      loop_safe++; 
      DEBUG_ME;

      // If it's a trijunction, add it.
      if (node->size() > 2) cell_path.push_back(node);

      // Loop around CW, until the next edge on the border.
      // If removing, go CCW (false).  But also be careful about the region ID.
      this_cell->next_edge_in_region(node, this_edge->id, next_cell, next_edge, CW, id);

      // node = CW ? next_edge->nb : next_edge->na; // that will be our next node
      node = next_edge->nb;  // that will be our next node

      this_edge = next_edge; // edge, and
      this_cell = next_cell; // We have to do cells, because the tracts are not necessarily contiguous.
    }

    cout << "node ring BEFORE :: ";
    for (auto n : node_ring) {
      if (n == snip_a || n == snip_b) cout << "**";
      cout << n->id;
      if (n == snip_a || n == snip_b) cout << "**";
      cout << " ";
    }
    cout << endl;

    cout << "snip_a=" << snip_a->id << "(" << snip_a << ")" << " snip_b=" << snip_b->id << "(" << snip_b << ")" << "  number_on_cell=" << number_on_cell << " cell_path=" << cell_path.size() << endl;
    if (snip_a == snip_b && number_on_cell > 1) {
      cell_path.clear();
      number_on_cell++;
    }

    int removed = 0;
    if (number_on_cell > 2) {

      if (insertion_point+number_on_cell-1 <= (int) node_ring.size()) {

        node_ring.erase(node_ring.begin()+insertion_point+1,
                        node_ring.begin()+insertion_point+number_on_cell-1);
        removed = number_on_cell - 2;

      } else {

        int unremoved = (insertion_point+number_on_cell-1) - node_ring.size();

        DEBUG_ME; cout << "insertion_point=" << insertion_point << "  number_on_cell=" << number_on_cell << "  unremoved=" << unremoved << "  and node_ring.size()=" << node_ring.size() << endl;
        if (insertion_point+1 < (int) node_ring.size()) 
          node_ring.erase(node_ring.begin()+insertion_point+1, node_ring.end());
        node_ring.erase(node_ring.begin(), node_ring.begin()+unremoved);
        insertion_point = node_ring.size()-1; // insert at the end.
        cout << "new insertion point=" << insertion_point << endl;
      }
    }

    cout << "node ring AFTER removing " << removed << " :: ";
    for (auto n : node_ring) {
      if (n == snip_a || n == snip_b) cout << "**";
      cout << n->id;
      if (n == snip_a || n == snip_b) cout << "**";
      cout << " ";
    }
    cout << endl;

    // If there was originally just a corner on this cell,
    // we need to add another copy of this node to return to.
    DEBUG_ME;  cout << "node_ring.size()=" << node_ring.size() << ", cell_path.size()=" << cell_path.size() << endl;
    if (number_on_cell == 1 && cell_path.size()) cell_path.push_back(node_ring[insertion_point]);

    DEBUG_ME; cout << "about to insert at " << insertion_point+1 << " on node_ring size=" << node_ring.size()
                   << " adding " << cell_path.size() << " nodes." << endl;
    node_ring.insert(node_ring.begin()+insertion_point+1, cell_path.begin(), cell_path.end());
    DEBUG_ME;

    cout << "Adding to region " << id << " cell=" << c->id << " with snipA=" << snip_a->id << " with snipB=" << snip_b->id
         << "\n  > cell_path size=" << cell_path.size() << ", elements :: ";
    for (auto n : cell_path) cout << n->id << " ";
    cout << endl;

    cout << "node ring AFTER insertion :: ";
    for (auto n : node_ring) {
      if (n == snip_a || n == snip_b) cout << "**";
      cout << n->id;
      if (n == snip_a || n == snip_b) cout << "**";
      cout << " ";
    }
    cout << endl;

    Node* last = node_ring.back();
    for (auto n : node_ring) {
      assert(last != n);
      last = n;
    }

    return;

  }

  float Region::update_scc(Cell* add = 0, Cell* sub = 0, bool UPDATE = false) { 

    // Avoid if at all possible
    if ((!add || add->d2(x_scc, y_scc) + eps_scc <= r2_scc) &&
        (!sub || sub->d2(x_scc, y_scc) + eps_scc <= r2_scc)) return r2_scc;

    // For now, we're updating as needed with Bernd Gaertner's Miniball.
    // https://people.inf.ethz.ch/gaertner/subdir/software/miniball.html
    // If the above condition fails, we know (xi, yi) is one of the 
    // support points, which dramatically changes the optimal algorithm.
    // On the other hand, if you remove a point, you don't know that one
    // of the remaining points is on the border:
    //  .
    // :             x
    //  .
    
    // move data structure as a list of vectors
    // -->> faster to expose the data directly??
    std::list<std::vector<double> > lp;
    for (auto c : cells) {
      if (sub && sub->id == c->id) continue;
      lp.push_back(std::vector<double>{c->x, c->y});
    }
    if (add) lp.push_back(std::vector<double>{add->x, add->y});

    MB mb(2, lp.begin(), lp.end()); 
    
    if (UPDATE) { 
      r2_scc = mb.squared_radius();
      x_scc  = mb.center()[0];
      y_scc  = mb.center()[1];
    }

    return mb.squared_radius();

  }

  float Region::update_lic(Cell* add = 0, Cell* sub = 0, bool UPDATE = false) { 

    if (!node_ring.size()) {
      cerr << "Found a ring with size 0 in region " << id << endl;
      throw std::runtime_error(std::string("Found a ring with size 0"));
    }

    if (!add && !sub && !UPDATE) return r2_lic;

    bool in_circle = false;

    if (add) for (auto n : add->nodes) {
      if (n->d2(x_lic, y_lic) < r2_lic + eps_lic) in_circle = true;
      if (in_circle) break;
    }

    if (sub) for (auto n : sub->nodes) {
      if (n->d2(x_lic, y_lic) < r2_lic + eps_lic) in_circle = true;
      if (in_circle) break;
    }

    if (!in_circle && !UPDATE) return r2_lic;

    int add_reg = -1, sub_reg = -1;
    if (add) { add_reg = add->region; add->region = id; }
    if (sub) { sub_reg = sub->region; sub->region = -1; }

    std::vector<Node*> ring;
    get_node_ring(ring);

    if (add) { add->region = add_reg; }
    if (sub) { sub->region = sub_reg; }


    bg_poly poly;
    std::vector<bp_pt> poly_pts;
    bpoly::voronoi_diagram<double> vd;

    for (auto n : ring) {
      bgeo::append(poly, n->pt);
      poly_pts.push_back(bp_pt(n->x*VOR_SCALE, n->y*VOR_SCALE));
    }
    bgeo::append(poly, ring.front()->pt);

    construct_voronoi(poly_pts.begin(), poly_pts.end(), &vd);

    double dist, x_max = xctr, y_max = yctr;
    double dist_max= boost::numeric::bounds<double>::lowest();

    for (auto it : vd.vertices()) {

      bp_pt pt_vor(it.x(), it.y());
      bp_pt pt_nom(it.x()/VOR_SCALE, it.y()/VOR_SCALE);

      dist = bpoly::euclidean_distance(poly_pts[it.incident_edge()->cell()->source_index()], pt_vor)/VOR_SCALE;

      if (dist_max < dist && bgeo::within(pt_nom, poly)) {
        dist_max = dist;
        x_max    = bgeo::get<0>(pt_nom);
        y_max    = bgeo::get<1>(pt_nom);
      }
    }

    if (UPDATE) { 
      r2_lic = dist_max*dist_max;
      x_lic  = x_max;
      y_lic  = y_max;
    }

    return dist_max*dist_max;

  }



  Universe::Universe() {}
  Universe::Universe(int n) 
    : pop(0), rcount(0), nregions(n), total_iterations(0), ALPHA(4), RANDOM(false), TRADE(true), TABU_LENGTH(0),
    DESTRAND_MIN(2), DESTRAND_MAX(15)
  {

    assert(nregions > 0);

    for (size_t ri = 0; ri < nregions; ri++) 
      regions.push_back(new Region(ri));
  }

  void Universe::add_cell(Cell c) {
  
    cells.push_back(new Cell(c));

    pop += c.pop;
    target = pop/nregions;

  }

  void Universe::add_edge(int cell_id, int edge_id, int nodea, int nodeb) {

    for (auto& c : cells) {
      if (c->id == cell_id) {
        c->add_edge(edge_id, nodea, nodeb);
        return;
      }
    }

    throw std::runtime_error(std::string("Failed to find cell for edge!"));
  }

  void Universe::add_node(int node_id, float x, float y) {

    nodes.push_back(new Node(node_id, x, y));

  }

  void Universe::add_node_edge(int node_id, int edge_id) {

    for (auto n : nodes) {
      if (n->id == node_id) {
        n->add_edge(edge_id);
        return;
      }
    }
  }

  void Universe::adjacency_to_pointers() {
    for (auto c : cells) c->adjacency_to_pointers(cells);
  }

  void Universe::node_ids_to_pointers() {
    for (auto r : regions) r->has_topo = true;
    for (auto c : cells)   c->node_ids_to_pointers(nodes);
  }

  void Universe::build_dijkstra_graph() {

    for (auto& c : cells) {
      for (auto& n : c->nm) {
        if (c->id < n.first->id) { // only make nodes one way...
          float w_edge = c->dist(n.first->x, n.first->y);
          boost::add_edge(c->id, n.first->id,
                          edge_weight_prop_t(w_edge),
                          dijkstra_graph);
        }
      }
    }

    for (auto& c : cells) {

      unsigned int start_id = c->id;
      bvertex start = vertex(start_id, dijkstra_graph);
      std::vector<bvertex> p(num_vertices(dijkstra_graph));

      dijkstra_shortest_paths(dijkstra_graph, start, boost::predecessor_map(&p[0]));

      for (auto& endc : cells) {
        int end_id = endc->id;
        bvertex v = end_id;
        while (p[v] != start_id) v = p[v];
        for (auto& n : c->nm) {
          if (v == ((unsigned int) n.first->id)) {
            c->dijkstra_step[end_id] = n.first;
            break;
          }
        }
      }
    }

  }


  std::vector<int> Universe::do_dijkstra(int start_id, int end_id) {

    std::vector<int> path;

    Cell *crawler = NULL;
    for (auto& c : cells) 
      if (c->id == start_id)
        crawler = c;

    if (!crawler) return path;

    while (crawler->id != end_id) {
      path.push_back(crawler->id);
      crawler = crawler->dijkstra_step[end_id];
    }
    path.push_back(end_id);
 
    return path;

  }


  int Universe::get_ncells() { return cells.size(); }

  void Universe::rand_init(int seed = 0) {

    if (cells.size() < nregions) {
      throw std::runtime_error(std::string("Cannot make districts: fewer cells than districts!"));
    }

    while (nregions < regions.size())
      regions.push_back(new Region(regions.size()-1));

    mersenne.seed(seed);
    int ncells = get_ncells();

    // Draw without replacement.
    int c = -1;
    std::unordered_set<int> cell_idx({-1});
    for (auto r : regions) {

      while (cell_idx.find(c) != cell_idx.end())
        c = mersenne() % ncells;

      cell_idx.insert(c);
      r->add_cell(cells[c]);
    }

  }

  std::pair<std::pair<float, float>, float> Universe::get_circle_coords(size_t rid, RadiusType rt) {

    if (rid >= regions.size()) return std::pair<std::pair<float, float>, float>();
    return regions[rid]->get_circle_coords(rt);

  }

  std::vector<std::pair<float, float> > Universe::get_point_ring(size_t rid) { 
    
    if (rid >= regions.size()) return std::vector<std::pair<float, float> >();
    return regions[rid]->get_point_ring();

  }

  // Function for testing.
  void Universe::add_cell_to_region(int cid, size_t rid) {

    if (rid >= regions.size()) {
      cout << "Warning :: rid > nregions.  Returning.\n";
      return;
    }

    for (auto c : cells) {
      if (c->id == cid) {

        if (!c->neighbors_connected()) {
          // cout << "Moving this cell would break continuity.  Returning.\n";
          return;
        }

        int rem = c->region;

        std::vector<Node*> rem_ring, add_ring;
        c->region = rid;
        regions[rem]->get_node_ring(rem_ring);
        regions[rid]->get_node_ring(add_ring);

        // regions[rem]->remove_cell(c);
        // regions[rid]->add_cell(c);

        regions[rem]->make_ring();
        regions[rid]->make_ring();


        auto iter = find(add_ring.begin(), add_ring.end(), regions[rid]->node_ring.front());
        if (iter != add_ring.end())
          std::rotate(add_ring.begin(), iter, add_ring.end());


        if (add_ring != regions[rid]->node_ring) {

          cout << "Mismatch!" << endl;

          cout << "Default :: ";
          for (auto n : regions[rem]->node_ring) {
            bool on_cell = (c->nodes.find(n) != c->nodes.end());
            if (on_cell) cout << "**";
            cout << n->id;
            if (on_cell) cout << "**";
            cout << " ";
          }
          cout << endl;

          cout << "Add ring :: ";
          for (auto n : add_ring) {
            bool on_cell = (c->nodes.find(n) != c->nodes.end());
            if (on_cell) cout << "**";
            cout << n->id;
            if (on_cell) cout << "**";
            cout << " ";
          }
          cout << endl;
      
        } else cout << "Huzzah!!!" << endl;

        c->region = rem;

        // cout << "Running divert for regions " << rid << " and " << rem << endl;
        // regions[rem]->divert_ring_at_cell(c, false);
        // regions[rid]->divert_ring_at_cell(c, true);

        // cout << "running make_ring for regions " << rid << " and " << rem << endl;
        // regions[rid]->make_ring();
        // regions[rem]->make_ring();

        regions[rid]->remove_cell(c);
        regions[rem]->add_cell(c);
      }
    }

  }

  std::map<int, int> Universe::cell_region_map() {

    std::map<int, int> crm;
    for (auto c : cells) {
      crm[c->id] = c->region;
      for (auto l : c->enclaves_and_islands) 
        crm[l] = c->region;
    }

    return crm;

  }

  std::vector<int> Universe::border_cells(bool EXT = false, int rid = -1) {

    std::vector<int> bc;
    for (auto r : regions) {
      if (rid >= 0 && r->id != rid) continue;

      if (EXT) for (auto const& c : r->ext_borders) bc.push_back(c->id);
      else     for (auto const& c : r->int_borders) bc.push_back(c->id);
    }

    return bc;

  }

  void Universe::connect_graph() {

    std::vector<std::unordered_set<Cell*> > graphs;
    for (auto c : cells) {

      // First make the local star graph.
      std::unordered_set<Cell*> lone_star;
      lone_star.insert(c);
      for (auto n : c->nm) lone_star.insert(n.first);

      // Then get a vector of all the existing
      // graphs that it's connected to:
      // if any element in lone star is in an 
      // existing graph, just mark it.
      std::vector<int> gidx;
      for (size_t gi = 0; gi < graphs.size(); gi++) {
        for (auto c : lone_star) {
          if (graphs[gi].find(c) != graphs[gi].end()) {
              gidx.push_back(gi);
              break;
          }
        }
      } // gid should now hold a list of overlapping graphs

      // If it has (yet) no neighbors, add it on its own.
      if (!gidx.size()) { 
        graphs.push_back(lone_star);

      // Otherwise, add it to the first,
      // and add _others_ on to that one.
      } else {
        graphs[gidx[0]].insert(lone_star.begin(), lone_star.end());

        // backwards, to keep removals consistent.
        for (auto gi : boost::adaptors::reverse(gidx)) {
          if (gi == gidx[0]) continue;
          graphs[gidx[0]].insert(graphs[gi].begin(), graphs[gi].end());
          graphs.erase(graphs.begin() + gi);
        }
      }
    }

    if (graphs.size() == 1) {
      cout << "Congratulations, it was already connected!" << endl;
      return;
    }


    while (graphs.size() > 1) {

      std::sort(graphs.rbegin(), graphs.rend(), cellp_set_len_compare); 

      // We'll work on the back element (which we'll pop, afterwards).

      // One up from the reverse end, till one before the first.
      // Search for the closest "foreign" cell and its graph.
      int for_min_g = -1;
      Cell* loc_min_c = 0;
      Cell* for_min_c = 0;
      double min_dist2 = std::numeric_limits<double>::infinity();
      for (int fg = int(graphs.size())-2; fg >= 0; fg--) {
        cout << "Looking at a foreign graph!" << endl;

        for (auto fc : graphs[fg]) {
          for (auto lc : graphs.back()) {

            double d2 = lc->d2(fc->x, fc->y);
            if (d2 < min_dist2) {
              min_dist2 = d2;
              loc_min_c = lc;
              for_min_c = fc;
              for_min_g = fg;
            }
          }
        }
      }

      // link the closest foreign and local cells;
      // add the contents of the smaller graph
      // to the larger one, then pop off the graph.
      loc_min_c->nm.push_back(std::make_pair(for_min_c, 1.));
      for_min_c->nm.push_back(std::make_pair(loc_min_c, 1.));
      for (auto lc : graphs.back()) graphs[for_min_g].insert(lc);
      cout << "lc" << loc_min_c->id << " to "
           << "fc" << for_min_c->id << "   "
           << "d=" << sqrt(min_dist2) << endl;
      cout << "Added a graph of size " << graphs.back().size() 
                << " to a graph of size " << graphs[for_min_g].size() << endl;

      graphs.pop_back();

    }

    cout << "ending with : ";
    for (auto g : graphs) cout << g.size() << " ";
    cout << endl;

  }

  std::vector<int> Universe::clipped_cells() {

    std::vector<int> cc;
    for (auto c : cells) {
      for (auto l : c->enclaves_and_islands) {
        cc.push_back(l);
      }
    }

    return cc;
  }


  void Universe::grow_kmeans(int popgrow = 0) {

    bool growth = true;
    while (growth) {

      growth = false;

      double d2;

      std::sort(regions.begin(), regions.end(), pop_compare);
      for (auto r : regions) {

        Cell* b_min_c = 0; double b_min_d2 = std::numeric_limits<double>::infinity();
        for (auto const& b : r->ext_borders) {
          if (b->region >= 0) continue;

          if (popgrow) d2 = b->d2(r->xpctr, r->ypctr);
          else         d2 = b->d2(r->xctr,  r->yctr);

          if (d2 < b_min_d2) {
            b_min_d2 = d2;
            b_min_c = b;
          }
        }

        if (b_min_c) {
          r->add_cell(b_min_c);
          growth = true;
          break;
        }
      }

    }

    std::sort(regions.begin(), regions.end(), id_compare);

    for (auto r : regions) if (r->has_topo) {
      r->make_ring();
      r->update_lic(0, 0, true);
    }

  }

  void Universe::load_partition(std::map<int, int> reg_map) {

    for (auto& c : cells) 
      regions[reg_map[c->id]]->add_cell(c);

    for (auto r : regions) if (r->has_topo) {
      r->make_ring(); 
      r->update_lic(0, 0, true);
    }

  }

  void Universe::iterate(int niter = 1, float tol = 0.05, int r = -1) {


    if (r >= int(regions.size())) return;

    for (int i = 0; i < niter; i++) { // The number of iterations.

      for (auto rit : regions) { // over all regions...

        if (r >= 0 && rit->id != r) continue;

        int rem_reg = 0;
        Cell* b_min_c = 0; double b_min_d2 = std::numeric_limits<double>::infinity();
        for (auto b : rit->ext_borders) {

          double d2 = b->d2(rit->xctr, rit->yctr);

          if (d2 < b_min_d2 &&
              rit->pop < regions[b->region]->pop && 
              b->neighbors_connected()) {
            b_min_d2 = d2;
            b_min_c = b;
            rem_reg = b->region;
          }
        }

        if (b_min_c) {
          regions[rem_reg]->remove_cell(b_min_c);
          rit->add_cell(b_min_c);
        }


        Cell* b_max_c = 0; int b_dest_reg = -1;
        double b_max_d2 = -std::numeric_limits<double>::infinity();
        for (auto b : rit->int_borders) {

          double d2 = b->d2(rit->xctr, rit->yctr);

          int dreg = -1;
          for (auto n : b->nm) if (n.first->region != rit->id) {
            dreg = n.first->region;
            break;
          }

          if (dreg >= 0 && d2 > b_max_d2 &&
              rit->pop > regions[dreg]->pop && 
              b->neighbors_connected()) {
            b_dest_reg = dreg;
            b_max_d2 = d2;
            b_max_c = b;
          }
        }

        if (b_max_c) {
          rit->remove_cell(b_max_c);
          regions[b_dest_reg]->add_cell(b_max_c);
        }

      }

    }
      
    return;

  }


  double Region::obj(ObjectiveMethod omethod, Cell* add = 0, Cell* sub = 0, bool verbose = false) {

    if (verbose) cout << "In objective fn for region " << id << ", to ";
    if (verbose && add) cout << "add cid" << add->id << "(r" << add->region << ")";
    if (verbose && sub) cout << "sub cid" << sub->id << "(r" << sub->region << ")";
    if (verbose) cout << endl;

    switch (omethod) {

      case ObjectiveMethod::DISTANCE_A: return obj_distance (add, sub, xctr,  yctr ); 
      case ObjectiveMethod::DISTANCE_P: return obj_distance (add, sub, xpctr, ypctr); 
      case ObjectiveMethod::INERTIA_A:  return obj_inertia_a(add, sub, verbose); 
      case ObjectiveMethod::INERTIA_P:  return obj_inertia_p(add, sub, verbose); 
      case ObjectiveMethod::REOCK:      return obj_reock    (add, sub, verbose);
      case ObjectiveMethod::HULL_A:     return obj_hull     (add, sub, verbose);
      case ObjectiveMethod::POLSBY:     return obj_polsby   (add, sub, verbose);
      case ObjectiveMethod::PATH_FRAC:  return obj_path_frac(add, sub, verbose);
      case ObjectiveMethod::EHRENBURG:  return obj_ehrenburg(add, sub, verbose);
      case ObjectiveMethod::POLSBY_W:   // weight by population -- units... ???

      default:
        cout << "Not yet implemented." << endl;
        return 0.;
    }

  }

  double Region::obj_distance(Cell* add, Cell* sub, float xx, float yy) {

    // This is a borderline normalization for area-based,
    // but not appropriate for the population-weighted version.
    float r = area / M_PI;

    if (add) return   r / add->d2(xx, yy);
    if (sub) return - r / sub->d2(xx, yy);

    return 0;

  }

  double Region::obj_inertia_a(Cell* add = 0, Cell* sub = 0, bool verbose = false) {
    
    // Keep these variables for comparison.
    float xmod(xctr * area), ymod(yctr * area), amod(area);

    if (add) {
      xmod += add->area * add->x;
      ymod += add->area * add->y;
      amod += add->area;
    }

    if (sub) {
      xmod -= sub->area * sub->x;
      ymod -= sub->area * sub->y;
      amod -= sub->area;
    }

    xmod /= amod; ymod /= amod;


    float onom = 0, omod = 0;
    for (auto c : cells) {
      onom += c->area * c->d2(xctr, yctr);
      omod += c->area * c->d2(xmod, ymod);
    }

    if (add) omod += add->area * add->d2(xmod, ymod);
    if (sub) omod -= sub->area * sub->d2(xmod, ymod);

    if (verbose) {
      cout << " >> NOM  (x, y)=(" << xctr << ", " << yctr << ")    onom=" << onom << ",  normalized: " << onom /((area * area) / (2 * M_PI)) << "\n"  
           << " >> MOD  (x, y)=(" << xmod << ", " << ymod << ")    omod=" << omod << ",  normalized: " << omod /((amod * amod) / (2 * M_PI)) 
           << endl;
    }

    onom /= (area * area) / (2 * M_PI);
    omod /= (amod * amod) / (2 * M_PI);
    // omod /= onom;

    if (verbose) {
      cout << " >> onom - omod = " << onom - omod << " (positive means improvement.)  \n";
      if (add) cout << " >> add" << add->id << " coords at (x, y) = (" << add->x << ", " << add->y << ")" << endl;
      if (sub) cout << " >> sub" << sub->id << " coords at (x, y) = (" << sub->x << ", " << sub->y << ")" << endl;
    }

    onom -= omod;

    return onom;

  }


  double Region::obj_inertia_p(Cell* add = 0, Cell* sub = 0, bool verbose = false) {
    
    // Keep these variables for comparison.
    int pmod(pop); 
    float xmod(xpctr * pop), ymod(ypctr * pop), amod(area);

    if (add) {
      xmod += add->pop * add->x;
      ymod += add->pop * add->y;
      pmod += add->pop;
      amod += add->area;
    }

    if (sub) {
      xmod -= sub->pop * sub->x;
      ymod -= sub->pop * sub->y;
      pmod -= sub->pop;
      amod -= sub->area;
    }

    xmod /= pmod; ymod /= pmod;


    float onom = 0, omod = 0;
    for (auto c : cells) {
      onom += c->pop * c->d2(xpctr, ypctr);
      omod += c->pop * c->d2(xmod,  ymod );
    }

    if (add) omod += add->pop * add->d2(xmod, ymod);
    if (sub) omod -= sub->pop * sub->d2(xmod, ymod);

    if (verbose) {
      cout << " >> NOM  a=" << area << "  pop=" << pop  << " (x, y)=(" << xpctr << ", " << ypctr << ")    onom=" << onom << ",  normalized: " << onom /((area * pop ) / (2 * M_PI)) << "\n"  
           << " >> MOD  a=" << amod << "  pop=" << pmod << " (x, y)=(" << xmod  << ", " << ymod  << ")    omod=" << omod << ",  normalized: " << omod /((amod * pmod) / (2 * M_PI)) 
           << endl;
    }

    onom /= (area * pop ) / (2 * M_PI);
    omod /= (amod * pmod) / (2 * M_PI);

    if (verbose) {
      cout << " >> onom - omod = " << onom - omod << " (positive means improvement.)  \n";
      if (add) cout << " >> add" << add->id << " coords at (x, y) = (" << add->x << ", " << add->y << ")" << endl;
      if (sub) cout << " >> sub" << sub->id << " coords at (x, y) = (" << sub->x << ", " << sub->y << ")" << endl;
    }

    onom -= omod;

    return onom;

  }


  double Region::obj_reock(Cell* add = 0, Cell* sub = 0, bool verbose = false) {

    float r2_mod = update_scc(add, sub);
    float amod   = area;

    if (add) amod += add->area;
    if (sub) amod -= sub->area;

    float reock_nom = area / (M_PI * r2_scc);
    float reock_mod = amod / (M_PI * r2_mod);

    if (verbose) {
      cout << "r_scc: " << sqrt(r2_scc) << "  area: " << area << "  reock: " << reock_nom << endl;
      if (add) cout << "r_mod add: " << sqrt(r2_mod) << "  area: " << amod << "  reock: " << reock_mod << endl;
      if (sub) cout << "r_mod sub: " << sqrt(r2_mod) << "  area: " << amod << "  reock: " << reock_mod << endl;
    }

    return reock_mod - reock_nom; // Higher number is improving.

  }

  double Region::obj_ehrenburg(Cell* add = 0, Cell* sub = 0, bool verbose = false) {

    if ((add && !add->neighbors_connected()) ||
        (sub && !sub->neighbors_connected())) {
      return boost::numeric::bounds<float>::lowest();
    }

    float r2_mod = update_lic(add, sub);
    float amod   = area;

    if (add) amod += add->area;
    if (sub) amod -= sub->area;

    float ehrenburg_nom = M_PI * r2_lic / area;
    float ehrenburg_mod = M_PI * r2_mod / amod;

    if (verbose) {
      cout << "r_lic nom: " << sqrt(r2_lic) << "  area: " << area << "  ehrenburg: " << ehrenburg_nom << endl;
      if (add) cout << "r_mod add: " << sqrt(r2_mod) << "  area: " << amod << "  ehrenburg: " << ehrenburg_mod << endl;
      if (sub) cout << "r_mod sub: " << sqrt(r2_mod) << "  area: " << amod << "  ehrenburg: " << ehrenburg_mod << endl;
    }

    return ehrenburg_mod - ehrenburg_nom; // Higher number is improving.

  }


  double Region::obj_hull(Cell* add = 0, Cell* sub = 0, bool verbose = false) {

    double hull_nom = area / ch_area;

    float     amod = area;
    float   chamod = ch_area;
    bg_mpt mpt_mod = mpt;

    if (add) amod += add->area;
    if (sub) amod -= sub->area;

    // could be faster to just rebuild.
    if (sub && !bgeo::within(sub->pt, ch_poly)) {
      bg_mpt mpt_diff;
      bgeo::difference(mpt, sub->pt, mpt_diff);
      mpt_mod = mpt_diff;
    }

    if (add && !bgeo::within(add->pt, ch_poly)) {
      bgeo::append(mpt_mod, add->pt);
    }

    if (add || sub) {
      bg_poly ch_poly_mod;
      bgeo::convex_hull(mpt_mod, ch_poly_mod);
      chamod = bgeo::area(ch_poly_mod);
    }

    double hull_mod = amod / chamod;

    if (verbose) {
      if (add) cout << "ADDED   :: nom=" << std::setprecision(5) << hull_nom << "=" << area << "/" << ch_area << "  mod=" << hull_mod << "=" << amod << "/" << chamod << endl;
      if (sub) cout << "SUBB'ED :: nom=" << std::setprecision(5) << hull_nom << "=" << area << "/" << ch_area << "  mod=" << hull_mod << "=" << amod << "/" << chamod << endl;
    }

    return hull_mod - hull_nom;

  }

  double Region::obj_path_frac(Cell* add, Cell* sub, bool verbose) {

    if (verbose) cout << "In obj_path_frac..." << endl;

    // Yes, this is necessary....
    long long num_nom = 0, den_nom= 0;
    long long num_mod = 0, den_mod= 0;

    // If we're adding, add every other pair....
    if (add) {
      int end_id = add->id;
      for (auto crawler : cells) {

        int pop_prod = (add->pop) * (crawler->pop);
        // cout << __LINE__ << " :: " << pop_prod << endl;
        den_mod += pop_prod;

        while (crawler->id != end_id) {
          if (crawler->region != id) break;
          if (crawler == sub)        break;

          crawler = crawler->dijkstra_step[end_id];
        }

        if (crawler->id == end_id) num_mod += pop_prod;
      }
    }


    // Now visit all of the other pairs of cells.
    bool mod_contained, nom_contained;
    for (auto& dest : cells) {

      int end_id = dest->id;
      for (auto crawler : cells) {

        if (dest->id < crawler->id) continue;

        nom_contained = true;
        mod_contained = (sub != crawler);

        int pop_prod = (dest->pop) * (crawler->pop);
        // cout << __LINE__ << " :: dest_pop=" << dest->pop << " crawler_pop=" << crawler->pop << endl;
        // cout << __LINE__ << " :: " << pop_prod << endl;

        den_nom += pop_prod;
        // cout << __LINE__ << " " << den_nom << " + " << pop_prod << " = " << den_mod + pop_prod << endl;
        if (mod_contained) den_mod += pop_prod;

        while (crawler->id != end_id && 
               (mod_contained || nom_contained)) {

          if (crawler->region != id) nom_contained = false;
          if ((crawler != add && crawler->region != id) || 
               crawler == sub) mod_contained = false;

          crawler = crawler->dijkstra_step[end_id];
        }

        if (nom_contained) num_nom += pop_prod;
        if (mod_contained) num_mod += pop_prod;
      }
    }

    if (verbose) cout << "Nom fraction " << 1.*num_nom/den_nom << endl;
    if (verbose) cout << "Mod fraction " << 1.*num_mod/den_mod << endl;

    double frac_mod = 1. * num_mod / den_mod;
    double frac_nom = 1. * num_nom / den_nom;

    return frac_mod - frac_nom;

  }

  double Region::obj_polsby(Cell* add = 0, Cell* sub = 0, bool verbose = false) {

    double polsby_nom = 4 * M_PI * area / sumw_border / sumw_border; // 4*pi*A/P^2

    double area_mod(area), perim_mod(sumw_border);

    // Iterate over the MOVING cell's neighbors.
    // If we're adding this cell, outside neighbors 
    // make the border longer
    if (add) { 
      area_mod += add->area;
      for (auto const& n : add->nm) {
        if (n.first->region != id) perim_mod += n.second;
        else                       perim_mod -= n.second;
      }
    }

    // while it's the inverse for subtracting.
    if (sub) {
      area_mod -= sub->area;
      for (auto const& n : sub->nm) {
        if (n.first->region != id) perim_mod -= n.second;
        else                       perim_mod += n.second;
      }
    }


    double polsby_mod = 4 * M_PI * area_mod / perim_mod / perim_mod;

    if (verbose) {
      if (add) cout << "ADDED   :: nom=" << std::setprecision(5) << polsby_nom << " :: A=" << area << " P=" << sumw_border << "     mod=" << polsby_mod << "=" << " :: A=" << area_mod << " P=" << perim_mod << endl;
      if (sub) cout << "SUB'ED  :: nom=" << std::setprecision(5) << polsby_nom << " :: A=" << area << " P=" << sumw_border << "     mod=" << polsby_mod << "=" << " :: A=" << area_mod << " P=" << perim_mod << endl;
    }

    return polsby_mod - polsby_nom;

  }

  float Region::d2(float x, float y, RadiusType rt) {

    std::pair< std::pair<float, float>, float > circ = get_circle_coords(rt);

    float xR = circ.first.first;
    float yR = circ.first.second;

    return (x-xR)*(x-xR) + (y-yR)*(y-yR);
  }

  void Universe::transfer_strand(std::unordered_set<Cell*>& strand) {

    if (strand.size()) cout << "Removing strand of size " << strand.size() << endl;

    while (strand.size()) {

      bool removed = false;
      for (auto c : strand) {

        std::unordered_set<int> nei_reg = c->neighbor_regions();
        // Move in an order that preserves contiguity.
        if (nei_reg.size() && c->neighbors_connected()) { 
          regions[c->region]->remove_cell(c); 
          regions[*nei_reg.begin()]->add_cell(c);
          strand.erase(c);
          set_tabu(c);
          removed = true;
          break;
        }
      }

      if (!removed) {
        cout << "WARNING :: strand not correctly removed!" << endl;
        break;
      }
    }
  }

  int Universe::destrand(int min = 1, int max = 1e9) {

    std::unordered_set<Cell*> opt_strand;

    for (auto r : regions) {
      for (auto b : r->ext_borders) {

        std::unordered_set<Cell*> strand;
        int strand_max = max;
        if (2*max > regions[b->region]->ncells) max = regions[b->region]->ncells/2-1;
        b->neighbor_strands(strand, strand_max);

        bool strand_is_tabu = false;
        for (auto c : strand) 
          if (is_tabu(c)) strand_is_tabu = true;
        if (strand_is_tabu) continue;

        if (strand.size() > opt_strand.size()) {
          opt_strand = strand;
        }
      }
    }

    int nremoved = opt_strand.size();
    if (nremoved > min) {
      transfer_strand(opt_strand);
      return nremoved;
    } else return 0;

  }


  bool Universe::trade(Region* ir, RadiusType rt) {

    std::unordered_map<int, bool> ib_connected;
    for (auto& ib : ir->int_borders)
      ib_connected[ib->id] = ib->neighbors_connected();

    float best = 0;
    Cell *add = 0, *sub = 0;
    for (auto& eb : ir->ext_borders) {

      if (!eb->neighbors_connected(false, true)) continue;

      for (auto& ib : ir->int_borders) {

        Region* er = regions[eb->region];
        if (er->ext_borders.find(ib) == er->ext_borders.end()) continue;

        if (!ib->neighbors_connected(false, true)) continue; // ROOK, FAST
        if (eb->is_neighbor(ib)) continue; // Don't try this.

        float nom = er->dist(eb->x, eb->y, rt) + ir->dist(ib->x, ib->y, rt);
        float mod = er->dist(ib->x, ib->y, rt) + ir->dist(eb->x, eb->y, rt);

        if (mod - nom < best) {
          best = mod - nom;
          add = eb; sub = ib;
        }
      }
    }

    if (best < 0) {
      Region* er = regions[add->region];
      er->remove_cell(add); ir->add_cell(add);
      ir->remove_cell(sub); er->add_cell(sub);
      return true;
    }

    return false;

  }
  

  bool Universe::greedy(Region* rit, ObjectiveMethod omethod, float tol, float best_move, bool random, int r, bool verbose) {
    
    Cell* b_opt_c = 0;
    std::unordered_set<Cell*> opt_strands;

    // only iterate over adding external; 
    // int. cells will be removed as someone else's external.
    if (random) { // Random greedy is grasp.

      std::vector<std::unordered_set<Cell*>::iterator> ebv(rit->ext_borders.size());
      std::iota(ebv.begin(), ebv.end(), rit->ext_borders.begin());
      std::shuffle(ebv.begin(), ebv.end(), mersenne);

      for (auto ebit : ebv) {

        Cell* b = *ebit;
        if (r >= 0 && b->region != r && rit->id != r) continue;
        if (greedy_evaluate(rit, b, tol, omethod, best_move, b_opt_c, opt_strands, verbose)) break;
      }

    } else {

      for (auto b : rit->ext_borders) {

        if (r >= 0 && b->region != r && rit->id != r) continue;
        greedy_evaluate(rit, b, tol, omethod, best_move, b_opt_c, opt_strands, verbose);

      }
    }


    if (b_opt_c) {
      int rem_reg = b_opt_c->region;

      // up to four cells....
      transfer_strand(opt_strands);

      regions[rem_reg]->remove_cell(b_opt_c);
      rit->add_cell(b_opt_c);

      set_tabu(b_opt_c);

      return true;

    }
        
    return false;

  }


  bool Universe::greedy_evaluate(Region* r, Cell* b, float tol, ObjectiveMethod omethod,
                                 float& best_move, Cell*& b_opt_c, std::unordered_set<Cell*>& opt_strands,
                                 bool verbose) {

     if (is_tabu(b)) return false; 

     float dpop  = r->pop / target - 1;
     float dbpop = regions[b->region]->pop / target - 1;

     float dp_ij = sign(dbpop - dpop) * pow(fabs(dpop - dbpop)/tol, ALPHA);

     float dF_ij = r                 ->obj(omethod, b, 0, verbose) +
                   regions[b->region]->obj(omethod, 0, b, verbose);

     float dO_ij = dp_ij + dF_ij;

     if (dO_ij < best_move) return false;
     if (regions[b->region]->ncells == 1) return false;

     std::unordered_set<Cell*> strands;
     int nsets = b->neighbor_strands(strands, DESTRAND_MAX, false);
     if (nsets != 1 && strands.size() < DESTRAND_MIN) return false;
     if (is_tabu_strand(strands)) return false;

     best_move = dO_ij;
     b_opt_c = b;
     opt_strands = strands;

     return true;

  }


  void Universe::oiterate(ObjectiveMethod omethod, int niter = 1, float tol = 0.05, int seed = 0, int r = -1, int verbose = false) {

    if (r >= int(regions.size()))
      throw std::runtime_error(std::string("Cannot iterate over region -- there aren't that many!"));

    for (int i = 0; i < niter; i++) { // The number of iterations.

      mersenne.seed(total_iterations + seed); // Different splits of loops etc. should be reproducible.
      total_iterations++;

      if (!(i%10)) cout << "iteration " << i << endl;

      for (auto rit : regions) { // over all regions...

        if (TRADE) trade(rit, obj_radius[omethod]);

        greedy(rit, omethod, tol, 0, RANDOM, r, verbose);

      }
    }

    return;
  }



  int Universe::merge_strands(Cell* c, int max = 10) {

    std::vector<std::unordered_set<Cell*> > graphs;
    if (c->neighbor_sets(graphs, max, true) == 1) return 0;

    int removed = 0;
    for (auto s : graphs) {
      if (int(s.size()) > max) continue;

      for (auto nc : s) {

        c->merge(nc);

        cells.erase(std::remove(cells.begin(), cells.end(), nc), cells.end());
        delete nc;

        removed++;

      }
    }

    return removed;

  }


  void Universe::trim_graph() {

    int clipped = 0 ;

    // Trivial case: single strands.  
    // Do this one the easy way.
    std::vector<Cell*>::iterator cit  = cells.begin();
    while (cit != cells.end()) {
      if ((*cit)->nm.size() == 1) {

        (*cit)->nm.begin()->first->merge(*cit);
        delete *cit;

        cit = cells.erase(cit);
        clipped++;

      } else ++cit;
    }

    // int breaker = 0;
    for (int ci = 0; ci < int(cells.size()); ci++) {

      int merged = merge_strands(cells[ci], 10);

      // backtrack by the max cells we could have lost...
      ci -= merged;

      // but increment...
      clipped += merged;

    }
 
  }


  bool cellp_set_len_compare(std::unordered_set<Cell*> s1, std::unordered_set<Cell*> s2) {
    return s1.size() < s2.size();
  }

  bool pop_compare(Region* r1, Region* r2) { return r1->pop < r2->pop; }
  bool id_compare(Region* r1, Region* r2) { return r1->id < r2->id; }

  int sign(double x) { return (x > 0) - (x < 0); }


  bool maxed_or_empty(unsigned int max, std::vector<std::unordered_set<Cell*> >& graphs, 
                                        std::vector<std::unordered_set<Cell*> >& to_add) {

    if (graphs.size() != to_add.size()) 
      throw std::runtime_error(std::string("In maxed_or_empty, graphs and to_add must have the same length"));

    int ngraphs = graphs.size();
    for (int ni = 0; ni < ngraphs; ni++) {
      if (graphs[ni].size() < max && to_add[ni].size()) return false;
    }

    return true;
  }

}


