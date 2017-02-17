# distutils: language = c++
# distutils: sources = Cluscious.cpp

from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.map cimport map
# from libcpp.unordered_map cimport unordered_map
from libcpp.unordered_set cimport unordered_set

cdef extern from "Cluscious.h" namespace "Cluscious" :

    cpdef enum ObjectiveMethod:
        DISTANCE_A, DISTANCE_P, INERTIA_A, INERTIA_P, HULL_A, POLSBY, REOCK, EHRENBURG, POLSBY_W, PATH_FRAC

    cdef cppclass Cell:
        Cell() except + 
        Cell(Cell) except + 
        Cell(int, int, double, double, double, map[int, double], string) except +
        int id, pop
        double x, y, area
        map[int, double] wm

    cdef cppclass Region:
        Region() except + 
        Region(int, Cell*) except + 
        Region(int, double, double) except +
        double xctr, yctr

        # double Region(Objective Method, obj) 

    cdef cppclass Universe:
        Universe() except + 
        Universe(int) except + 
        int pop, nregions
        double target
        vector[Cell*] cells

        void add_cell(Cell)
        int  get_ncells()
        map[int, int] cell_region_map()
        vector[int] border_cells(int)
        vector[int] clipped_cells()

        void rand_districts(int) except +
        void contiguity_to_neighbors()
        void connect_graph()
        void trim_graph()
        void grow_kmeans(int popgrow)
        void iterate(int, float, int)

        void oiterate(ObjectiveMethod, int, float, float, int, int)




cdef class cell:
    cdef Cell c_cell # hold a C++ instance which we're wrapping
    def __cinit__(self, int i, int p, double x, double y, double a, map[int, double] wm, str wkt = "MULTIPOLYGON()"):
        self.c_cell = Cell(i, p, x, y, a, wm, wkt.encode('utf-8'))

    def __str__(self):
        return "<pycl::cell id={} (x,y)=({:.02f},{:.02f})>".format(self.c_cell.id, self.c_cell.x, self.c_cell.y)

    property id:
        def __get__(self): return self.c_cell.id
        def __set__(self, id): self.c_cell.id = id

    property pop:
        def __get__(self): return self.c_cell.pop
        def __set__(self, pop): self.c_cell.pop = pop

    property x:
        def __get__(self): return self.c_cell.x
        def __set__(self, x): self.c_cell.x = x

    property y:
        def __get__(self): return self.c_cell.y
        def __set__(self, y): self.c_cell.y = y

    property area:
        def __get__(self): return self.c_cell.area
        def __set__(self, a): self.c_cell.area = a

    property wm:
        def __get__(self): return self.c_cell.wm
        def __set__(self, wm): self.c_cell.wm = wm




cdef class region:

    cdef Region c_region # hold a C++ instance which we're wrapping
    def __cinit__(self, id, x, y = None):
        if not y: self.c_region = Region(id, x.x, x.y)
        else:     self.c_region = Region(id, x, y)

    property xctr:
        def __get__(self): return self.c_region.xctr
        # def __set__(self, x): self.c_region.x = x

    property yctr:
        def __get__(self): return self.c_region.yctr
        # def __set__(self, y): self.c_region.y = y


                           
cdef class universe:

    cdef Universe c_univ # hold a C++ instance which we're wrapping

    def __cinit__(self, int n):
        self.c_univ = Universe(n)

    def add_cell(self, cell c):
        self.c_univ.add_cell(c.c_cell)

    def get_ncells(self):
        return self.c_univ.get_ncells()

    def rand_districts(self, seed = 0):
        return self.c_univ.rand_districts(seed)

    def connect_graph(self):
        self.c_univ.connect_graph()

    def trim_graph(self):
        self.c_univ.trim_graph()

    def contiguity_to_neighbors(self):
        return self.c_univ.contiguity_to_neighbors()

    def get_cell(self, int c):

        for ci in self.c_univ.cells:
          if ci.id == c:
            return cell(ci.id, ci.pop, ci.x, ci.y, ci.area, ci.wm)

        return None

    def cell_region_map(self):
        return self.c_univ.cell_region_map()

    def border_cells(self, int rid = -1):
        return self.c_univ.border_cells(rid)

    def clipped_cells(self):
        return self.c_univ.clipped_cells()

    def grow_kmeans(self, int popgrow = False):
        self.c_univ.grow_kmeans(popgrow)

    def iterate(self, int niter = 1, float tol = 0.05, int r = -1):
        self.c_univ.iterate(niter, tol, r)

    def oiterate(self, int om_i = 0, int niter = 1, float tol = 0.05, float alpha = 4, int reg = -1, int verbose = 0):
        self.c_univ.oiterate(ObjectiveMethod(om_i), niter, tol, alpha, reg, verbose)


    property pop:
        def __get__(self): return self.c_univ.pop

    property target:
        def __get__(self): return self.c_univ.target

    property nregions:
        def __get__(self): return self.c_univ.nregions

    property cells:
        def __get__(self):
          return [c.id for c in self.c_univ.cells]

