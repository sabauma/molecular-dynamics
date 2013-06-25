
include ../Makefile.inc

ALL : Simulation.o Statistics.o
	echo Making Simulation

Simulation.o : Simulation.cpp Simulation.h
	$(CPP) $(PRJCPPFLAGS) -c Simulation.cpp

Statistics.o : Statistics.cpp Statistics.h
	$(CPP) $(PRJCPPFLAGS) -c Statistics.cpp

clean :
	echo cleaning up in $(CURDIR)
	-@rm *.o 2>/dev/null || true
