#ifndef __SIMULATION_COLLISION_EVENT_H__
#define __SIMULATION_COLLISION_EVENT_H__

#include "../Vector/Vector.hpp"

struct CollisionRecord
{
    DoubleVector Position;
    double Time;
    double DeltaV;
    double Energy;
    double Distance;
    int Type1;
    int Type2;
};

#endif
