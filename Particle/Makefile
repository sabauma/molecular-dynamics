
include ../Makefile.inc

ALL : Particle.o
	echo Making Particle

Particle.o : Particle.cpp
	$(CPP) $(PRJCPPFLAGS) -c Particle.cpp

clean :
	echo cleaning up in $(CURDIR)
	-@rm *.o 2>/dev/null || true
