#include <string>
#include <vector>
#include <iostream> 


class Rectangle {

  public:
    int x0, y0, x1, y1;
    std::vector<int> v;
    
    Rectangle(int x0, int y0, int x1, int y1);
    int getArea();

    void incr(int i);
    std::vector<int> myVec() { return v; }

};

