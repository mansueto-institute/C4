# cython: profile=True
# cython: linetrace=True
# distutils: language=c++
# distutils: sources=Cluscious.cpp

from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.pair cimport pair
# from libcpp.unordered_map cimport unordered_map
from libcpp.unordered_set cimport unordered_set

cdef extern from "Cluscious.h" namespace "Cluscious" :

    cpdef enum ObjectiveMethod:
        DISTANCE_A, DISTANCE_P, INERTIA_A, INERTIA_P, HULL_A, POLSBY, REOCK, EHRENBURG, POLSBY_W, PATH_FRAC

    cpdef enum RadiusType:
        EQUAL_AREA, EQUAL_AREA_POP, EQUAL_CIRCUMFERENCE, SCC, LIC

    cdef cppclass Cell:
        Cell() except + 
        Cell(Cell) except + 
        Cell(int, int, double, double, double, map[int, double], int, int) except +
        void add_edge(int, int, int)
        int id, pop
        double x, y, area
        map[int, double] wm

    cdef cppclass Region:
        Region() except + 
        Region(int) except +
        Region(int, Cell*) except + 
        double xctr, yctr

    cdef cppclass Universe:
        Universe() except + 
        Universe(int) except + 
        int pop, nregions
        double target

        int ALPHA
        int RANDOM
        int TRADE
        int TABU_LENGTH
        int DESTRAND_MAX
        int DESTRAND_MIN

        vector[Cell*] cells
        vector[Region*] regions;

        void add_cell(Cell)
        void add_edge(int, int, int, int)
        void add_node(int, float, float)
        void add_node_edge(int, int)
        int  get_ncells()
        map[int, int] cell_region_map()
        vector[int] border_cells(int, int) 
        vector[int] clipped_cells()

        vector[pair[float, float]] get_point_ring(int) 
        pair[pair[float, float], float] get_circle_coords(int, RadiusType)
        void add_cell_to_region(int cid, int rid)

        void rand_init(int) except +
        vector[int] do_dijkstra(int, int)
        void build_dijkstra_graph()
        void adjacency_to_pointers()
        void node_ids_to_pointers()
        void connect_graph()
        void trim_graph()
        void grow_kmeans(int popgrow)
        void load_partition(map[int, int] rmap)
        void iterate(int, float, int)
        int  destrand(int, int)

        void oiterate(ObjectiveMethod, int, float, int, int, int);




cdef class cell:
    cdef Cell c_cell # hold a C++ instance which we're wrapping
    def __cinit__(self, int i, int p, double x, double y, double a, map[int, double] wm, int edge, int split):
        self.c_cell = Cell(i, p, x, y, a, wm, edge, split)

    def add_edge(self, int eid, int nodea, int nodeb):
        self.c_cell.add_edge(eid, nodea, nodeb)

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
        self.c_region = Region(id)

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

    def add_edge(self, int cell_id, int eid, int nodea, int nodeb):
        self.c_univ.add_edge(cell_id, eid, nodea, nodeb)

    def add_node(self, int nid, float x, float y):
        self.c_univ.add_node(nid, x, y)

    def add_node_edge(self, int nid, int eid):
        self.c_univ.add_node_edge(nid, eid)

    def get_ncells(self):
        return self.c_univ.get_ncells()

    def rand_init(self, seed = 0):
        return self.c_univ.rand_init(seed)

    def connect_graph(self):
        self.c_univ.connect_graph()

    def trim_graph(self):
        self.c_univ.trim_graph()

    def build_dijkstra_graph(self):
        self.c_univ.build_dijkstra_graph()

    def do_dijkstra(self, a, b):
        return self.c_univ.do_dijkstra(a, b)

    def adjacency_to_pointers(self):
        self.c_univ.adjacency_to_pointers()

    def node_ids_to_pointers(self):
        self.c_univ.node_ids_to_pointers()

    def get_cell(self, int c):

        for ci in self.c_univ.cells:
          if ci.id == c:
            return cell(ci.id, ci.pop, ci.x, ci.y, ci.area, ci.wm)

        return None

    def get_circle_coords(self, int rid, int rt_i):
        return self.c_univ.get_circle_coords(rid, RadiusType(rt_i))

    def get_point_ring(self, rid):
        return self.c_univ.get_point_ring(rid)

    def add_cell_to_region(self, cid, rid):
        self.c_univ.add_cell_to_region(cid, rid)

    def cell_region_map(self):
        return self.c_univ.cell_region_map()

    def border_cells(self, int ext = 0, int rid = -1):
        return self.c_univ.border_cells(ext, rid)

    def clipped_cells(self):
        return self.c_univ.clipped_cells()

    def grow_kmeans(self, int popgrow = False):
        self.c_univ.grow_kmeans(popgrow)

    def load_partition(self, rmap):
        if type(rmap) is dict:
          self.c_univ.load_partition(rmap)
        elif type(rmap) is str:
          rmapd = {}
          for line in open(rmap):
            sp = [int(x) for x in line.split(",")]
            rmapd[sp[0]] = sp[1]
          self.c_univ.load_partition(rmapd)

    def destrand(self, mini = None, maxi = None):
        if not mini: mini = 1
        if not maxi: maxi = 1e9
        return self.c_univ.destrand(mini, maxi)

    def iterate(self, int niter = 1, float tol = 0.05, int r = -1):
        self.c_univ.iterate(niter, tol, r)

    def oiterate(self, int om_i = 0, int niter = 1, float tol = 0.05, int seed = 0, int reg = -1, int verbose = 0):
        self.c_univ.oiterate(ObjectiveMethod(om_i), niter, tol, seed, reg, verbose)


    property pop:
        def __get__(self): return self.c_univ.pop

    property target:
        def __get__(self): return self.c_univ.target

    property nregions:
        def __get__(self): return self.c_univ.nregions

    property cells:
        def __get__(self):
          return [c.id for c in self.c_univ.cells]

    property ALPHA:
        def __get__(self): return self.c_univ.ALPHA
        def __set__(self, a): self.c_univ.ALPHA = a

    property RANDOM:
        def __get__(self): return self.c_univ.RANDOM
        def __set__(self, l): self.c_univ.RANDOM = l

    property TRADE:
        def __get__(self): return self.c_univ.TRADE
        def __set__(self, l): self.c_univ.TRADE = l

    property TABU_LENGTH:
        def __get__(self): return self.c_univ.TABU_LENGTH
        def __set__(self, l): self.c_univ.TABU_LENGTH = l

    property DESTRAND_MIN:
        def __get__(self): return self.c_univ.DESTRAND_MIN
        def __set__(self, l): self.c_univ.DESTRAND_MIN = l

    property DESTRAND_MAX:
        def __get__(self): return self.c_univ.DESTRAND_MAX
        def __set__(self, l): self.c_univ.DESTRAND_MAX = l


