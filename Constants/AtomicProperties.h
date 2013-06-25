/**
 * Constants for the program corresponding to atomic/molecular properties of
 * the various particle types in the simulation. All values correspond to real
 * world units unless otherwise stated.
 *
 * Force constants were obtained from "Molecular theory of gasses and liquids"
 * by Hirscfelder, Joseph.
 */

#ifndef __CONSTANTS_ATOMIC_PROPERTIES__
#define __CONSTANTS_ATOMIC_PROPERTIES__

#include <map>
#include <string>

#include "./Helper.h"

namespace Constants
{
    /** The number of different particle types in the simulation. */
    constexpr int ELEMENT_TYPES = 8;

    /** The number of ionization levels tracked. */
    constexpr int IONIZATION_LEVELS = 9;

    /** The names of the elements. */
    const std::string ELEMENT_NAMES[ELEMENT_TYPES] =
    {
        "H"  ,
        "D"  ,
        "H2" ,
        "D2" ,
        "He" ,
        "Ne" ,
        "Ar" ,
        "Xe"
    };

    const std::map<std::string, int> TYPE_FROM_NAME
        = invert_array_mapping<std::string>(ELEMENT_NAMES, ELEMENT_TYPES);

    /** The name of the liquid. */
    const std::string LIQUID_NAME = "H20";

    /** Flag controlling whether or not an element can fuse. */
    constexpr bool FUSABLE[ELEMENT_TYPES] =
    {
        true  , // H
        true  , // D
        false , // H2
        false , // D2
        false , // He
        false , // Ne
        false , // Ar
        false   // Xe
    };

    const bool DISSOCIATES[ELEMENT_TYPES] =
    {
        false , // H
        false , // D
        true  , // H2
        true  , // D2
        false , // He
        false , // Ne
        false , // Ar
        false   // Xe
    };

    constexpr double DISSOCIATION_ENERGY[ELEMENT_TYPES] =
    {
        1.0e100 , // H
        1.0e100 , // D
        436.0   , // H2
        443.6   , // D2
        1.0e100 , // He
        1.0e100 , // Ne
        1.0e100 , // Ar
        1.0e100   // Xe
    };

    const std::string DISSOCIATES_TO_NAME[ELEMENT_TYPES] =
    {
        "",  // H -- No dissociation
        "",  // D -- No dissociation
        "H", // H2 -> 2H
        "D", // D2 -> 2D
        "" , // He -- No dissociation
        "" , // Ne -- No dissociation
        "" , // Ar -- No dissociation
        "" , // Xe -- No dissociation
    };

    /**
     * The entry (i,j) gives the energy reqiured to extract the jth
     * electron from element i. These do not need to be converted into
     * simulation units. They are in kJ/mol, which has a 1-1 conversion
     * to energy units in the simulation. Not sure why that is, off the top
     * of my head, but that is how the numbers work out.
     */
    constexpr double IONIZATION_ENERGY_INI[ELEMENT_TYPES][IONIZATION_LEVELS] =
    {
        {1312.0 , 1.0e100 , 1.0e100 , 1.0e100 , 1.0e100 , 1.0e100 , 1.0e100 , 1.0e100 , 1.0e100} , // H
        {1312.0 , 1.0e100 , 1.0e100 , 1.0e100 , 1.0e100 , 1.0e100 , 1.0e100 , 1.0e100 , 1.0e100} , // D
        {1488.0 , 1.0e100 , 1.0e100 , 1.0e100 , 1.0e100 , 1.0e100 , 1.0e100 , 1.0e100 , 1.0e100} , // H2
        {1488.0 , 1.0e100 , 1.0e100 , 1.0e100 , 1.0e100 , 1.0e100 , 1.0e100 , 1.0e100 , 1.0e100} , // D2 -- approx
        {2372.3 , 5250.5  , 1.0e100 , 1.0e100 , 1.0e100 , 1.0e100 , 1.0e100 , 1.0e100 , 1.0e100} , // He
        {2080.7 , 3952.3  , 6122.0  , 9371.0  , 12177.0 , 1.0e100 , 1.0e100 , 1.0e100 , 1.0e100} , // Ne
        {1520.6 , 2665.8  , 3931.0  , 5771.0  , 7238.0  , 8781.0  , 11995.0 , 13842.0 , 1.0e100} , // Ar
        {1170.4 , 2046.4  , 3099.4  , 4600.0  , 5760.0  , 6930.0  , 9460.0  , 10800.0 , 1.0e100}   // Xe
    };

    /**
     * The atomic mass of particles of the liquid.
     */
    constexpr double LIQUID_ATOMIC_MASS_INI = 18.02;

    /**
     * The atomic radius of the particles of the liquid.
     */
    constexpr double LIQUID_ATOMIC_DIAMETER_INI = 2 * 2.82e-10;

    /**
     * The "initial" atomic mass of each element. This gives the mass of each
     * element in atomic mass units, before they are g to simulation
     * units.
     */
    constexpr double ATOMIC_MASS_INI[ELEMENT_TYPES] =
    {
        1.008, // H
        2.013, // D
        2.016, // H2
        4.028, // D2
        4.003, // He
        20.18, // Ne
        39.95, // Ar
        131.3  // Xe
    };

    /**
     * The diameter of each element in meters. These values are used only in
     * the hard sphere model
     */
    constexpr double ATOMIC_DIAMETER_INI[ELEMENT_TYPES] =
    {
        2 * 1.09e-10, // H  -- Not correct
        2 * 1.09e-10, // D  -- Not correct
        2 * 1.09e-10, // H2 -- Not correct
        2 * 1.09e-10, // D2 -- Not correct
        2 * 1.09e-10, // He
        2 * 1.30e-10, // Ne
        2 * 1.83e-10, // Ar
        2 * 2.46e-10  // Xe
    };

    /**
     * The constant ω (omega) corresponding to each particle type.
     */
    constexpr double BIRD_OMEGA_INI[ELEMENT_TYPES] =
    {
        0.188, // H -- Incorrect, using diatomic constant
        0.199, // D -- Incorrect, using diatomic constant
        0.188, // H2
        0.199, // D2
        0.204, // He
        0.160, // Ne
        0.201, // Ar
        0.302  // Xe
    };

    /**
     * The constant α (alpha) corresponding to each particle type.
     */
    constexpr double BIRD_ALPHA_INI[ELEMENT_TYPES] =
    {
        1.396, // H -- Incorrect, using diatomic constant
        1.420, // D -- Incorrect, using diatomic constant
        1.396, // H2
        1.420, // D2
        1.431, // He
        1.336, // Ne
        1.425, // Ar
        1.651  // Xe
    };

    /**
     * The constant μ (mu) corresponding to each particle type.
     */
    constexpr double BIRD_MU_INI[ELEMENT_TYPES] =
    {
        0.832, // H -- Incorrect, using diatomic constant
        1.185, // D -- Incorrect, using diatomic constant
        0.832, // H2
        1.185, // D2
        1.865, // He
        2.975, // Ne
        2.117, // Ar
        2.107  // Xe
    };

    /**
     * The constant "D" corresponding to each particle type. This corresponds
     * to the diameter of the particle, but is used for the VSS model whereas
     * the other diameter value above is used for the hard sphere model.
     *
     * Technically, this is the expected (i.e. mean) diameter of a particular
     * gas species at the reference temperature.
     *
     * These can be computed by simulating an equilibriated bubble at the
     * reference temperature and taking the mean of the distances between
     * colliding particles of that gas species.
     *
     * These values are only useful for certain logistics involving the cell
     * partitioning. In fact, only the largest gas is truely relevant, so the
     * rest can be nonsense, because they are not really used.
     */
    constexpr double BIRD_D_INI[ELEMENT_TYPES] =
    {
        2.30e-10, // H -- Not Correct
        2.30e-10, // D -- Not Correct
        2.30e-10, // H2 -- Not Correct
        2.30e-10, // D2 -- Not Correct
        2.30e-10, // He
        2.72e-10, // Ne
        4.11e-10, // Ar
        5.65e-10  // Xe
    };

    /**
     * The initial reference temperature.
     */
    const double BIRD_T_REF           = 273.0;
}

#endif /* __ATOMIC_PROPERTIES_H__ */
