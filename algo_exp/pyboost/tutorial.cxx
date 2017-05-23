#include "tutorial.h"


Rectangle::Rectangle(int X0, int Y0, int X1, int Y1) {

	x0 = X0; x1 = X1; y0 = Y0; y1 = Y1;

}

int Rectangle::getArea() {

	return (x1 - x0) * (y1 - y0);

}

void Rectangle::incr(int i) {

	std::cout << "i=" << i << std::endl;
  v.push_back(i);
  for (int a = 0; a < v.size(); a++)
    std::cout << v[a] << " ";
  std::cout << std::endl;

}


#include <boost/python.hpp>

BOOST_PYTHON_MODULE(tutorial)
{

	  using namespace boost::python;

    class_<Rectangle>("Rectangle", init<int, int, int, int>())
        .def("getArea", &Rectangle::getArea)
        .def("myVec",   &Rectangle::myVec)
        .def("incr",    &Rectangle::incr)
    ;
}
