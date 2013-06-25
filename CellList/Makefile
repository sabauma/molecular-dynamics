
include ../Makefile.inc

ALL: CellListFast.o
	echo Making CellList
	$(CPP) $(PRJCPPFLAGS) -c test.cpp
	$(CPP) $(PRJCPPFLAGS) CellListFast.o test.o

# CellList.o : CellList.cpp CellList.h
# 	$(CPP) $(PRJCPPFLAGS) -c CellList.cpp

CellListFast.o : CellListFast.cpp CellListFast.h
	$(CPP) $(PRJCPPFLAGS) -c CellListFast.cpp

clean :
	echo cleaning up in $(CURDIR)
	-@rm *.o 2>/dev/null || true

