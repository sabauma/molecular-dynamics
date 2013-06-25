
#include <algorithm>
#include <stdio.h>
#include <string>
#include <vector>

#include "../Constants/AtomicProperties.h"
#include "../Constants/Constants.h"
#include "../Exceptions/SimulationException.h"
#include "../Math/Functions.h"
#include "../Units/Units.h"
#include "InfoStruct.h"

using namespace std;
using Constants::PI;

static DriverMode
read_driver_mode(const std::string& str)
{
    std::string mode_name(str);
    std::transform(
            mode_name.begin(), mode_name.end(), mode_name.begin(), ::tolower);

    if (mode_name == "constant")
    {
        return Constant;
    }
    else if (mode_name == "acoustic")
    {
        return Acoustic;
    }
    throw SimulationException("Unknown driver mode: '" + str + "'");
}

InfoStruct::InfoStruct()
{ }

InfoStruct
InfoStruct::ToSimulationUnits() const
{
    InfoStruct converted;
    converted.LiquidType = this->LiquidType;

    converted.CopyConstants();

    // Conversion of liquid parameters
    converted.LiquidDensity      = this->LiquidDensity * Pow3(Units::L) / Units::M;
    converted.LiquidViscosity    = this->LiquidViscosity / Units::KinematicViscosity;
    converted.LiquidSpeedOfSound = this->LiquidSpeedOfSound * Units::T / Units::L;
    converted.LiquidSurfaceTension
        = this->LiquidSurfaceTension / Units::SurfaceTension;


    // Conversion of data read directly from the InfoStruct file.
    converted.AtomicParts = this->AtomicParts;
    converted.ConeAngle   = this->ConeAngle * Constants::PI / 180.0; // radians

    converted.Frequency   = this->Frequency * Units::T;

    converted.Pd = this->Pd * Constants::P0 / Units::Pressure;
    converted.R0 = this->R0 / Units::L;

    converted.AmbientPressure
        = this->AmbientPressure * Constants::P0 / Units::Pressure;

    converted.Radi = this->Radi * converted.R0;
    converted.Rini = this->Rini / Units::L;
    converted.Riso = this->Riso * converted.R0;

    converted.Temperature = this->Temperature / Units::Temperature;

    converted.Tini = this->Tini / Units::T;
    converted.Vini = this->Vini / (Units::L / Units::T);

    converted.WallTemp = this->WallTemp / Units::Temperature;

    const double tension_pressure = 2.0 * converted.LiquidSurfaceTension
                                  / converted.R0;

    // Density must factor in the pressure caused by surface tension to obtain
    // a stable bubble.
    converted.Density = (converted.AmbientPressure + tension_pressure)
                      * Pow3(this->R0 / this->Rini) / converted.Temperature ;

    // Convert parameters relevant to fusion
    converted.FusionBarrier = this->FusionBarrier / Units::Energy;
    converted.FusionWidth   = this->FusionWidth / Units::L;
    converted.PlankConstant = this->PlankConstant * Units::T
                            / (Units::L * Units::L * Units::M);

    // Convert AtomicMass
    std::transform(this->AtomicMass, this->AtomicMass + ET,         // Input
                   converted.AtomicMass,                            // Output
                   [this](double m) { return m / this->AtomicMass[ET - 1]; });

    std::transform(this->AtomicDiameter, this->AtomicDiameter + ET,
                   converted.AtomicDiameter,
                   [this](double d) { return d / this->AtomicDiameter[ET-1]; });

    // Convert VWallThermal
    std::transform(converted.AtomicMass, converted.AtomicMass + ET,
                   converted.VWallThermal,
                   [&converted](double mass)
                   { return sqrt(2.0 * converted.WallTemp / mass); });

    // Compute cone tan here, as it is not valid to compute before the
    // conversion.
    converted.ConeTan2 = Sqr(tan(converted.ConeAngle));
	converted.ConeSolidAngle = 2.0 * PI * (1.0 - cos(converted.ConeAngle));

    // Copy appropriate tables to be initialized from the converted data,
    // rather than original SI values.
    converted.InitTables();
    converted.Driver = this->Driver;

    return converted;
}

InfoStruct
InfoStruct::ParseInfoStruct(const std::string& fname)
{
    using namespace rapidjson;

    FILE* infile = fopen(fname.c_str(), "r");
    if (NULL == infile)
    {
        throw SimulationException("Could not open InfoStruct file.");
    }

    FileStream instream(infile);
    Document document;

    if (document.ParseStream<0, UTF8<> >(instream).HasParseError())
    {
        fclose(infile);
        throw SimulationException("Could not parse InfoStruct file.");
    }
    fclose(infile);

    InfoStruct is;
    double sum = 0.0;

    assert(document["Environment"].IsObject());
    Value& environment = document["Environment"];

    // Deterime the liquid type
    std::string liquid_name = environment["Liquid"].GetString();
    std::transform(liquid_name.begin(), liquid_name.end(),
                   liquid_name.begin(), ::tolower);

    // Environment information.
    is.LiquidType      = Constants::LIQUID_FROM_NAME.at(liquid_name);
    is.Driver          = read_driver_mode(environment["DriverMode"].GetString());
    is.ConeAngle       = environment["ConeAngle"].GetDouble();
    is.Temperature     = environment["Temperature"].GetDouble();
    is.WallTemp        = environment["WallTemp"].GetDouble();
    is.R0              = environment["R0"].GetDouble();
    is.Pd              = environment["Pd"].GetDouble();
    is.AmbientPressure = environment["AmbientPressure"].GetDouble();
    is.Frequency       = environment["Frequency"].GetDouble();
    is.Tini            = environment["Tini"].GetDouble();
    is.Rini            = environment["Rini"].GetDouble();
    is.Vini            = environment["Vini"].GetDouble();
    is.Radi            = environment["Radi"].GetDouble();
    is.Riso            = environment["Riso"].GetDouble();

    assert(document["Particles"]["AtomicParts"].IsObject());
    Value& parts = document["Particles"]["AtomicParts"];

    for (int i = 0; i < Constants::ELEMENT_TYPES; ++i)
    {
        if (parts.HasMember(Constants::ELEMENT_NAMES[i].c_str()))
        {
            is.AtomicParts.push_back(
                    parts[Constants::ELEMENT_NAMES[i].c_str()].GetDouble());
        }
        else
        {
            is.AtomicParts.push_back(0.0);
        }
    }

    // Sum the atomic percentages.
    sum = std::accumulate(is.AtomicParts.begin(), is.AtomicParts.end(), 0.0,
                          [](double x, double y) { return x + y; });
    assert(sum > 0.0);
    // Divide by the sum to ensure they total to 1.
    std::transform(is.AtomicParts.begin(), is.AtomicParts.end(),
                   is.AtomicParts.begin(),
                   [sum](double p) { return p / sum; });

    is.Density = is.AmbientPressure * Pow3(is.R0 / is.Rini) / is.Temperature;
    is.ConeTan2 = 0.0;

    // Copy the values from the tables
    is.CopyConstants();
    is.InitTables();

    return is;
}

void
InfoStruct::CopyConstants()
{
    LiquidDensity        = Constants::Density[LiquidType];
    LiquidViscosity      = Constants::Viscosity[LiquidType];
    LiquidSurfaceTension = Constants::SurfaceTension[LiquidType];
    LiquidSpeedOfSound   = Constants::SpeedOfSound[LiquidType];

    FusionBarrier        = Constants::FusionBarrier;
    FusionWidth          = Constants::FusionWidth;
    PlankConstant        = Constants::Plank;

    std::copy(Constants::ATOMIC_DIAMETER_INI,
              Constants::ATOMIC_DIAMETER_INI + ET,
              this->AtomicDiameter);

    std::copy(Constants::ATOMIC_MASS_INI,
              Constants::ATOMIC_MASS_INI + ET,
              this->AtomicMass);
}

void
InfoStruct::InitTables()
{
	// create look-up tables for d2 and half-reduced-mass
	for (int i = 0; i < ET; i++)
    {
		for (int j = i; j < ET; j++)
        {
			this->AtomicDiameterSquared[i][j]
                = this->AtomicDiameterSquared[j][i]
                = Sqr(0.5 * (this->AtomicDiameter[i] + this->AtomicDiameter[j]));

			this->ReducedMass[i][j]
                = this->ReducedMass[j][i]
                = 0.5 * this->AtomicMass[i] * this->AtomicMass[j]
                / (this->AtomicMass[i] + this->AtomicMass[j]);
		}
	}

    // Initialize the Bird constants
	for (int i = 0; i < ET; i++)
    {

		double mu = Constants::BIRD_MU_INI[i] * 1.0e-5 / Units::DynamicViscosity;
		this->BirdConstant[i] = 5.0
            * (Constants::BIRD_ALPHA_INI[i] + 1.0)
            * (Constants::BIRD_ALPHA_INI[i] + 2.0)
            * sqrt(1.0 / PI);

        // Having just temperature in there is valid, as Boltzmann's constant
        // would just be canceled out due to multiplying then dividing by
        // k for the energy conversion.
        this->BirdConstant[i] *=
            pow(Constants::BIRD_T_REF / Units::Temperature,
                Constants::BIRD_OMEGA_INI[i] + 0.5);

        this->BirdConstant[i] /=
            (16.0 * Constants::BIRD_ALPHA_INI[i]
             * exp(lgamma(4.0 - Constants::BIRD_OMEGA_INI[i])) * mu);
	}
}

std::ostream&
operator<<(std::ostream& stream, const InfoStruct& is)
{
    stream  << "ConeAngle: "      << is.ConeAngle      << '\n'
            << "ConeTan2: "       << is.ConeTan2       << '\n'
            << "ConeSolidAngle: " << is.ConeSolidAngle << '\n'
            << "Density: "        << is.Density        << '\n'
            << "Frequency: "      << is.Frequency      << '\n'
            << "Pd: "             << is.Pd             << '\n'
            << "R0: "             << is.R0             << '\n'
            << "Radi: "           << is.Radi           << '\n'
            << "Rini: "           << is.Rini           << '\n'
            << "Riso: "           << is.Riso           << '\n'
            << "Temperature: "    << is.Temperature    << '\n'
            << "Tini: "           << is.Tini           << '\n'
            << "Vini: "           << is.Vini           << '\n'
            << "WallTemp: "       << is.WallTemp       << std::endl;

    stream << "Atomic Parts: " << std::endl;
    for (unsigned int i = 0; i < is.AtomicParts.size(); ++i)
    {
        stream << "   " << Constants::ELEMENT_NAMES[i] << ": "
               << is.AtomicParts[i] << std::endl;
    }

    return stream;
}

void
InfoStruct::Read(const rapidjson::Value& reader)
{
    assert(reader["ConeAngle"].IsDouble());
    assert(reader["ConeTan2"].IsDouble());
    assert(reader["ConeSolidAngle"].IsDouble());
    assert(reader["Density"].IsDouble());
    assert(reader["Frequency"].IsDouble());
    assert(reader["Pd"].IsDouble());
    assert(reader["R0"].IsDouble());
    assert(reader["Radi"].IsDouble());
    assert(reader["Rini"].IsDouble());
    assert(reader["Riso"].IsDouble());
    assert(reader["Temperature"].IsDouble());
    assert(reader["Tini"].IsDouble());
    assert(reader["Vini"].IsDouble());
    assert(reader["WallTemp"].IsDouble());

    this->ConeAngle      = reader["ConeAngle"].GetDouble();
    this->ConeTan2       = reader["ConeTan2"].GetDouble();
    this->ConeSolidAngle = reader["ConeSolidAngle"].GetDouble();
    this->Density        = reader["Density"].GetDouble();
    this->Frequency      = reader["Frequency"].GetDouble();
    this->Pd             = reader["Pd"].GetDouble();
    this->R0             = reader["R0"].GetDouble();
    this->Radi           = reader["Radi"].GetDouble();
    this->Rini           = reader["Rini"].GetDouble();
    this->Riso           = reader["Riso"].GetDouble();
    this->Temperature    = reader["Temperature"].GetDouble();
    this->Tini           = reader["Tini"].GetDouble();
    this->Vini           = reader["Vini"].GetDouble();
    this->WallTemp       = reader["WallTemp"].GetDouble();

    this->CopyConstants();
    this->InitTables();
}

