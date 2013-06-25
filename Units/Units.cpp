
#include "../Constants/AtomicProperties.h"
#include "../Constants/Constants.h"
#include "Units.h"

using namespace Constants;

namespace Units
{
    const double Temperature = 120.0;

#ifdef HARD_SPHERE
    const double L = ATOMIC_DIAMETER_INI[ELEMENT_TYPES - 1];
#else
    const double L = BIRD_D_INI[ELEMENT_TYPES - 1]; // MD unit of length
#endif

    const double Energy = KB * Temperature;
    const double M      = 1.0e-3 * ATOMIC_MASS_INI[ELEMENT_TYPES - 1] / AVOGADRO;
    const double T      = L * sqrt(M / Energy);

    const double Pressure           = M / (L * T * T);
    const double DynamicViscosity   = M / (L * T);
    const double KinematicViscosity = L * L / T;
    const double SurfaceTension     = M / (T * T);
}
