#include <vector>


#include <boost/geometry/geometries/adapted/boost_polygon.hpp>


namespace c4 {
  
  namespace bpoly = boost::polygon;
  typedef bpoly::point_data<double> bp_pt;

  using std::cout;
  using std::endl;

  struct Node {

    Node(const int id_, const float x_, const float y_)
      : id(id_), pos(0), x(x_), y(y_), pt(x_, y_) {}

    bool set_edge(int edge) {
      auto it = find(edges.begin(), edges.end(), edge);
      if (it == edges.end()) return false;
      pos = distance(edges.begin(), it);
      return true;
    }

    int add_edge(int e) { 
      edges.push_back(e);
      return edges.size();
    }

    int size() { 
      return edges.size();
    }

    int next() {
      pos = (pos+1) % edges.size();
      return edges[pos];
    }

    int prev() {
      pos = (pos + edges.size() - 1) % edges.size();
      return edges[pos];
    }

    float d2(float xi, float yi) { return (x-xi)*(x-xi) + (y-yi)*(y-yi); }

    int id;
    int pos; // use pos instead of iter for easier cyclic.
    float x, y;
    bp_pt pt;
    std::vector<int> edges;

  };


  struct Edge {

    int id, na_id, nb_id;
    Node *na;
    Node *nb;
    Edge(const int id_, const int na_, const int nb_)
      : id(id_), na_id(na_), nb_id(nb_), na(0), nb(0) {}

    void set_na(Node* na_) { na = na_; }
    void set_nb(Node* nb_) { nb = nb_; }

  };

}


