PYVERSION = 3.5
PYROOT    = /home/jsaxon/anaconda3
PYINCLUDE = $(PYROOT)/include/python$(PYVERSION)m
PYLIB     = $(PYROOT)/lib/python$(PYVERSION)/config-$(PYVERSION)m/

# location of the Boost Python include files and library
BOOST_INC = $(PYROOT)/include
BOOST_LIB = $(PYROOT)/lib

# compile mesh classes
TARGET = tutorial

$(TARGET).so: $(TARGET).o
	g++ -shared -Wl,--export-dynamic $(TARGET).o -L$(BOOST_LIB) -lboost_python3 -L$(PYLIB) -lpython$(PYVERSION)m -o $(TARGET).so
	 
$(TARGET).o: $(TARGET).cxx
	g++ -I$(PYINCLUDE) -I$(BOOST_INC) -fPIC -c $(TARGET).cxx

