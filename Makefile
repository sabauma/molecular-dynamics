
include Makefile.inc

DIRS = CellList Constants Debug EventCalendar Exceptions InfoStruct ODE \
	   Particle Random RPK Simulation Units Vector

OBJS = CellList/CellListFast.o Debug/Trace.o \
	   EventCalendar/EventCalendar.o EventCalendar/EventType.o \
	   InfoStruct/InfoStruct.o Particle/Particle.o RPK/PressureGauge.o \
	   RPK/RPK.o Random/Random.o Simulation/Simulation.o \
	   Simulation/Statistics.o Units/Units.o

SRC = $(OBJS:.o=.cpp)

ALL : main.cpp
	# Compile each module in parallel. Each module compiles serially
	-for d in $(DIRS); do (echo $$d); done | xargs -L 1 --max-procs=16 make -C
	$(CPP) $(PRJCPPFLAGS) -c main.cpp
	$(CPP) $(PRJCPPFLAGS) $(OBJS) main.o -o sono.out $(LINK)

install : ALL
	cp sono.out ~/bin

ode:
	clang++  ode.cpp -o ode -I../ -I./

clean :
	-for d in $(DIRS); do (cd $$d; $(MAKE) clean ); done
	-@rm main.o 2>&1 | true
	-@rm sono.out 2>&1 | true
	-@rm wall_simulation.o 2>&1 | true
	-@rm wall_simulation 2>&1 | true

