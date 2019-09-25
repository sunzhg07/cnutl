#CC = g++
#PYTHON_VERSION = 3.7
CC=icpc
Inc_python = -I/usr/local/Cellar/python/3.7.4/Frameworks/Python.framework/Versions/3.7/include/python3.7m
Lib_python = -L/usr/local/Cellar/python/3.7.4/Frameworks/Python.framework/Versions/3.7/lib


Inc_boost = -I/usr/local/include
Lib_boost = -L/usr/local/lib

Lib_gsl = -L/usr/local/lib


GSL_LIBS=-lgsl -lgslcblas

SRC = $(wildcard *.cpp)
OBJ = $(addsuffix .o, $(basename $(SRC)))

# this target has to match whatever you have in BOOST_PYTHON_MODULE()
TARGET = Coupling_Module

# make this a so library -Wl,--export-dynamic
$(TARGET).so : $(OBJ)
	$(CC) -shared  $(OBJ) -o $(Inc_python) $(TARGET).so $(Lib_boost) -lboost_python37 $(Lib_python) -lpython3.7m $(Lib_gsl) $(GSL_LIBS)

%.o : %.cpp *.h
	$(CC) -fPIC -c $< $(Inc_python) $(Inc_boost) $(Lib_gsl) $(GSL_LIBS)
#

clean:
	-rm -f *.o *.so
