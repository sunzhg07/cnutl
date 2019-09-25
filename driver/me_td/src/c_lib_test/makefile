EXE=helloworld.so

#CC = icpc -std=c++11 -qopenmp
CC=icpc
Inc_python = -I/usr/local/Cellar/python/3.7.4/Frameworks/Python.framework/Versions/3.7/include/python3.7m/
Lib_python = -L/usr/local/Cellar/python/3.7.4/Frameworks/Python.framework/Versions/3.7/lib
#lib2 = -L$(HOME)/usr/lib

Lib_boost = -L/usr/local/lib
#Inc = -I/Users/zhuo/usr/include
#Inc = -I$(HOME)/usr/include


#GSL_LIBS=-lgsl -lgslcblas

OBJ=helloworld.o

#CFLAGS= -shared -Wl,-soname,#-fPIC -shared #-O3 -Wl,--rpath=$(HOME)/usr/lib

#$(CC) test.cpp $(CFLAGS)  -o test2.so $(Inc_python) $(Inc) $(lib2) -lboost_python
#g++ test.cpp -fPIC -shared -o test.so -I/usr/include/python3.4m -I/home_data/zhuo/usr/include -L/home_data/zhuo/usr/lib -lboost_python
$(EXE):$(OBJ)
	$(CC) -shared -o $(EXE) $(OBJ) $(Lib_python) -lpython3.7m $(Lib_boost) -lboost_python37 #$(GSL_LIBS)
#
helloworld.o : helloworld.cpp
	$(CC) -c -fPIC helloworld.cpp  $(Inc_python)

#g++ -shared -o helloworld.so -fPIC $(Inc_python) helloworld.cpp $(Lib_python) -lpython3.7m $(Lib_boost) -lboost_python37


clean:
	rm *.o *.t *.so *~ *.gch  *.dat ./serial_dat/*.dat #./data/*.dat
