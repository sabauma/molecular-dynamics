
include ../Makefile.inc

ALL : InfoStruct.o
	echo "Making InfoStruct"

InfoStruct.o : InfoStruct.cpp
	$(CPP) $(PRJCPPFLAGS) -c InfoStruct.cpp

clean :
	echo cleaning up in $(CURDIR)
	-@rm *.o 2>/dev/null || true
