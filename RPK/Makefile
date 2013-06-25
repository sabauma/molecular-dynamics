
include ../Makefile.inc

ALL : PressureGauge.o RPK.o
	echo Making RPK


PressureGauge.o : PressureGauge.cpp
	$(CPP) $(PRJCPPFLAGS) -c PressureGauge.cpp

RPK.o : RPK.cpp RPK.h PressureGauge.h
	$(CPP) $(PRJCPPFLAGS) -c RPK.cpp

clean :
	echo cleaning up in $(CURDIR)
	-@rm *.o 2>/dev/null || true
