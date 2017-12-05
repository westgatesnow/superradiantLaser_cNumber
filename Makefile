EIGENPATH = /usr/local/include/eigen3
GSLPATH = /usr/local/Cellar/gsl/2.2.1/include
INCPATH = -I/$(EIGENPATH) -I$(GSLPATH)
LIBPATH = -L/usr/local/lib -lgsl
CXXFLAGS = -c -O3
CXX=g++

beamLaser : beamLaser.o RNG.o
	$(CXX) -Wall $(LIBPATH) beamLaser.o RNG.o -o beamLaser

beamLaser.o : beamLaser.cpp 
	$(CXX) $(CXXFLAGS) $(INCPATH) beamLaser.cpp

RNG.o : RNG.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) RNG.cpp

clean :
	rm *.o beamLaser