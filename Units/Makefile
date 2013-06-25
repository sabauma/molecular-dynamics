
include ../Makefile.inc

ALL : Units.o
	echo Making Units

Units.o : Units.cpp
	$(CPP) $(PRJCPPFLAGS) -c Units.cpp

clean :
	echo cleaning up in $(CURDIR)
	-@rm *.o 2>/dev/null || true
