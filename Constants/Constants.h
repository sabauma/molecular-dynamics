/**
 * Constants referenced in the molecular dynamics simulation.
 */

#ifndef __CONSTANTS_CONSTANTS_H__
#define __CONSTANTS_CONSTANTS_H__

#include <boost/math/constants/constants.hpp>
#include <cmath>
#include <limits>
#include <map>

#include "../Vector/Vector.hpp"
#include "./Helper.h"

namespace Constants
{
    /** The number of dimensions in the simulation. */
    constexpr int DIMENSIONS = 3;

    /** Time value for an event that will never occur. */
    const double NEVER = std::numeric_limits<double>::infinity();

    /** The mathematical constant π. */
    const double PI = boost::math::constants::pi<double>();

    /** A particle index representing an invalid particle. */
    constexpr int NULL_PARTICLE = -1;

    /** Atmospheric pressure */
    constexpr double P0 = 1.0e5;                   // atmospheric pressure
    /** Avogadro's number */
    constexpr double AVOGADRO = 6.022045e23;
    /** Boltzmann's constexprant */
    constexpr double KB = 1.380662e-23;            // Boltzmann's constexpr
    /** ??? */
    constexpr double GAMMA = 5.0 / 3.0;
    /** Van Der Waal's constexprant */
    constexpr double ExcludedVolume = 0.00005105;  // Van Der Waal's constexprant

    /** Energy required to initiate fusion in Joules */
    constexpr double FusionBarrier = 4.5e7 * KB * 3.0 / 2.0;
    /** Width of the fusion barrier in meters */
    constexpr double FusionWidth = 1.0e-11;
    /** Plank's constant in m^2 kg / s */
    constexpr double Plank = 6.626068e-34;

    /** Unit Vectors for each direction. */
    const IntVector X(1, 0, 0);
    const IntVector Y(0, 1, 0);
    const IntVector Z(0, 0, 1);

    /** Standard basis vectors in R^3 */
    const IntVector STD_BASIS[3] =
    {
        X,
        Y,
        Z
    };

    /** Negates standard basis vectors in R^3 */
    const IntVector NEG_STD_BASIS[3] =
    {
        IntVector(-X),
        IntVector(-Y),
        IntVector(-Z)
    };

    constexpr int LIQUID_TYPES = 4;

    /** Names of the possible liquids in the simulation */
    const std::string LIQUID_NAMES[LIQUID_TYPES] =
    {
        "lithium",
        "mercury",
        "water",
        "diesel"
    };

    const std::map<std::string, int> LIQUID_FROM_NAME
        = invert_array_mapping<std::string>(LIQUID_NAMES, LIQUID_TYPES);

    /**
     * Kinematic viscosity of the fluids in m^2/s
     */
    const double Viscosity[LIQUID_TYPES] =
    {
        1.245e-6, // Lithium
        0.114e-6, // Mercury
        1.004e-6, // Water
        2.000e-6  // Diesel
    };

    /**
     * Surface tension of the liquids in N/m.
     */
    const double SurfaceTension[LIQUID_TYPES] =
    {
        0.396,   // Lithium
        0.485,   // Mercury
        0.0728,  // Water
        23.8e-3, // Diesel
    };

    /**
     * Density of the liquids in kg/m^3.
     */
    const double Density[LIQUID_TYPES] =
    {
        516,    // Lithium
        5430,   // Mercury
        999.07, // Water
        832.00  // Diesel
    };

    /**
     * Speed of sound of the liquids in m/s.
     */
    const double SpeedOfSound[LIQUID_TYPES] =
    {
        4490, // Lithium
        1450, // Mercury
        1497, // Water
        1250  // Diesel
    };
}

#endif /* __CONSTANTS_H__ */
