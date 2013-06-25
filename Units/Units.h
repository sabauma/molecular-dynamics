/**
 * Unit conversion code. The Units structure defined below is used to
 * transition real world units to units used in the MD simulation. This
 * defines an (ideally) singleton structure Units whose members
 */

#ifndef __UNITS_UNITS_H__
#define __UNITS_UNITS_H__

#include "../Constants/AtomicProperties.h"
#include "../Constants/Constants.h"

using namespace Constants;

namespace Units
{
    extern const double Temperature;

    extern const double L;

    extern const double Energy;
    extern const double M;
    extern const double T;

    extern const double Pressure;
    extern const double DynamicViscosity;
    extern const double KinematicViscosity;
    extern const double SurfaceTension;
}

#endif /* __UNITS_H__ */
