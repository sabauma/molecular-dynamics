
#include "./Constants.h"

/** Number of spatial dimenions in the simulation. */
static const int Constants::DIMENSIONS = 3;

/** A time that will never occur (barring error). */
static const int Constants::NEVER = std::numeric_limits<double>::infinity();

/** The number of dimensions in the simulation. */
static const int Constants::DIMENSIONS = 3;

/** Time value for an event that will never occur. */
static const double Constants::NEVER = std::numeric_limits<double>::infinity();

/** The mathematical constant π. */
static const double Constants::PI = boost::math::constants::pi<double>();

/** A particle index representing an invalid particle. */
static const int Constants::NULL_PARTICLE = -1;

/** Liquid viscosity. Presumably of the surrounding fluid. */
static const double Constants::VIS = 1.0e-3;
/** Liquid surface tension */
static const double Constants::SurfaceTension = 7.5e-2;      // liquid surface tension
/** Liquid density */
static const double Constants::Density = 9.990700e+002;      // Liquid density
/** Speed of sound in liquid. */
static const double Constants::C = 1.440700e+003;            // speed of sound in liquid

/** Atmospheric pressure */
static const double Constants::P0 = 1.0e5;                   // atmospheric pressure
/** Avogadro's number */
static const double Constants::AVOGADRO = 6.022045e23;
/** Boltzmann's constant */
static const double Constants::KB = 1.380662e-23;            // Boltzmann's const
/** ??? */
static const double Constants::GAMMA = 5.0 / 3.0;
/** Van Der Waal's constant */
static const double Constants::ExcludedVolume = 0.00005105;  // Van Der Waal's Constant

// Unit Vectors for each direction.
static const IntVector Constants::X(1, 0, 0);
static const IntVector Constants::Y(0, 1, 0);
static const IntVector Constants::Z(0, 0, 1);

static const IntVector Constants::STD_BASIS[3] =
{
    X,
    Y,
    Z
};

static const IntVector Constants::NEG_STD_BASIS[3] =
{
    IntVector(-X),
    IntVector(-Y),
    IntVector(-Z)
};

static const double Constants::BIRD_T_REF           = 273.0;

Constants::Constants()
{ }

Constants Constants::ParseConstants(const std::string& fname)
{
    Constants retval;
    retval.AtomicConstants = ParseAtomicProperties(fname);

    auto biggest = std::max_element(
            retval.AtomicConstants.begin(), retval.AtomicConstants.end(),
            [](const AtomicProperties& i, const AtomicProperties& j)
            {
                return i.AtomicMass < j.AtomicMass;
            });

    Units& units = retval.Conversions;

    units.Temperature = 120.0;
#ifdef VSS_PARAMETERS
    units.L = biggest->DiameterVSS;
#else
    units.L = biggest->AtomicDiameter;
#endif
    units.Energy = Constants::KB * units.Temperature;
    units.M      = 1.0e-3 * bigest->AtomicMass / Constants::AVOGADRO;
    units.T      = units.L * sqrt(units.M);

    units.Pressure       = units.M / (units.L * units.T * units.T);
    units.Viscosity      = units.M / (units.L * units.T);
    units.SurfaceTension = units.M / (units.T * units.T);

    return retval;
}
