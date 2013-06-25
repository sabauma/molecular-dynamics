
include ../Makefile.inc

ALL : EventCalendar.o EventType.o
	echo Making EventCalendar
	$(CPP) $(PRJCPPFLAGS) -c Test.cpp -o Test.o
	$(CPP) $(PRJCPPFLAGS) EventCalendar.o EventType.o Test.o

EventCalendar.o : EventCalendar.cpp
	$(CPP) $(PRJCPPFLAGS) -c EventCalendar.cpp

EventType.o : EventType.cpp
	$(CPP) $(PRJCPPFLAGS) -c EventType.cpp

clean :
	echo cleaning up in $(CURDIR)
	-@rm *.o 2>/dev/null || true
	-@rm a.out 2>/dev/null || true
