CXX = g++
CXXFLAGS = -O3 -Wall

BOOST_ROOT        := /usr/include/boost/python

PYTHON_INCL       := /usr/include/python3.8
PYTHON_LIBS       := /usr/lib/python3.8

LIBFLAGS          := -lboost_system  -lboost_python38 -lpython3.8


all : metrics

metrics: bin/metrics

clean :
	- rm -r bin

bin :
	mkdir bin

bin/metrics: metrics.cpp | bin
	$(CXX) metrics.cpp $(CXX_FLAGS) $(LIBFLAGS) -I$(PYTHON_INCL) -L$(BOOST_ROOT) -L$(PYTHON_LIBS) -o bin/metrics
