
# put common definitions in here
CPP         = g++
# CPP = clang++
# CPP = /usr/lib/gcc-snapshot/bin/g++

# -DBOOST_DISABLE_ASSERTS -DNDEBUG
PRJCPPFLAGS = -DTVMET_OPTIMIZE  -DBOOST_DISABLE_ASSERTS -DNDEBUG\
              -DUSE_PGAUGE -DVSS_PARAMETERS -DADIABATIC_INTERPOLATION \
              -DSPECULAR_WALL \
              -O3 -Wall -Wextra -std=c++11 -march=native -mtune=native \
              -flto -fwhole-program -fvisibility-inlines-hidden -ftree-vectorize -fomit-frame-pointer \
              -ftemplate-depth-100 -I../ -I./ -I./rapidjson/include -I../rapidjson/include

# PRJCPPFLAGS = -DTVMET_OPTIMIZE -DSPECULAR_WALL -DUSE_PGAUGE -DVSS_PARAMETERS \
#               -DADIABATIC_INTERPOLATION -O0 -g -ggdb -Wall -Wextra -std=c++11 \
#               -ftemplate-depth-100 -I../ -I./ -I./rapidjson/include -I../rapidjson/include

LINK = -lm

LD          = g++
RM          = rm
ECHO        = echo

SHELL	= /bin/sh

.SILENT :
