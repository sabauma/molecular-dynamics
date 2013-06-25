
include ../Makefile.inc

ALL : Trace.o
	echo Making Debug

Trace.o: Trace.cpp Trace.h
	$(CPP) $(PRJCPPFLAGS) -c Trace.cpp

clean :
	echo cleaning up in $(CURDIR)
	-@rm *.o 2>/dev/null || true
