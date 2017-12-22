EIGENPATH = /usr/local/include/eigen3
GSLPATH = /usr/local/Cellar/gsl/2.2.1/include
INCPATH = -I/$(EIGENPATH) -I$(GSLPATH)
LIBPATH = -L/usr/local/lib -lgsl
CXXFLAGS = -c -O3
CXX=g++

superradiantLaser : superradiantLaser.o RNG.o
	$(CXX) -Wall $(LIBPATH) superradiantLaser.o RNG.o -o superradiantLaser

superradiantLaser.o : superradiantLaser.cpp 
	$(CXX) $(CXXFLAGS) $(INCPATH) superradiantLaser.cpp

RNG.o : RNG.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) RNG.cpp

clean :
	rm *.o superradiantLaser