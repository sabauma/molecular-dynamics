
include ../Makefile.inc

ALL : Random.o
	echo Making Random

Random.o : Random.cpp Random.h
	$(CPP) $(PRJCPPFLAGS) -c Random.cpp

clean :
	echo cleaning up in $(CURDIR)
	-@rm *.o 2>/dev/null || true
