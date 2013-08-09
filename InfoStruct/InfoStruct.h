/**
 * Defines the structure used to hold the configurable constants of the
 * simulation. These are values read in from a file and must then be g
 * to simulation units.
 *
 * @author Spenser Bauman
 */

#ifndef __INFOSTRUCT_INFOSTRUCT_H__
#define __INFOSTRUCT_INFOSTRUCT_H__

#include "../Constants/AtomicProperties.h"
#include "rapidjson/document.h"
#include "rapidjson/filestream.h"
#include "rapidjson/prettywriter.h"
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

/**
 * Possible driver modes for use in the RPK equation. These modes define what
 * type of pressure function to use for the driving pressure.
 */
enum DriverMode
{
    /* Acoustic is a sinusoidal driver. */
    Acoustic,
    /* Constant is a constant pressure function. */
    Constant
};

class InfoStruct
{
public:

    // Convenient synonym for the number of elements.
    static const int ET = Constants::ELEMENT_TYPES;

    std::vector<double> AtomicParts;

    // The cone angle is half the apex angle of the cone.
    double ConeAngle;
    double ConeTan2;
    double ConeSolidAngle;
    double Density;
    double Frequency;
    double Pd;
    double AmbientPressure;
    double R0;
    double Radi;
    double Rini;
    double Riso;
    double Temperature;
    double Tini;
    double Vini;
    double WallTemp;

    int LiquidType;
    double LiquidDensity;
    double LiquidSpeedOfSound;
    double LiquidViscosity;
    double LiquidSurfaceTension;

    double FusionBarrier;

    double AtomicDiameter[ET];
    double AtomicDiameterSquared[ET][ET];

    double AtomicMass[ET];
    double ReducedMass[ET][ET];

    double VWallThermal[ET];

    double BirdConstant[ET];

    DriverMode Driver;

    double CollisionThreshold;

    /**
     * Empty constructor simply initializes the values to zero. They must be
     * set by the function deserializing the input and stored in this
     * container. Kept private as much of its data needs to be read in from
     * a file before values can be initialized.
     */
    InfoStruct();

    /**
     * Converts the current info struct to simulation units. Uses the
     * information provided by the Units module to convert the given
     * parameters to simulation units and returns a new <code>InfoStruct</code>
     * containing the simulation applicable constants.
     *
     * NOTE: This methdo makes the assumption that the current object is
     *       not in simulation units and is, in fact, using real world units.
     *
     * @return A new <code>InfoStruct</code> where the units are transformed
     *         to simulation units.
     */
    InfoStruct ToSimulationUnits() const;

    /**
     * Parser function to read in a JSON formatted text file and convert it
     * to an info struct. This makes use of the Boost Property Tree library
     * internally.
     *
     * @param fname The name of the file to be read.
     * @return The generated info struct.
     */
    static InfoStruct ParseInfoStruct(const std::string& fname);

    /**
     * Make this easier to output.
     */
    friend std::ostream& operator<<(std::ostream& stream, const InfoStruct& matrix);

    void Read(const rapidjson::Value& reader);

    template <typename Writer>
    void Serialize(Writer& writer) const;

private:

    /**
     * Initializes the tables for the current object using the current
     * values containted therein.
     */
    void InitTables();

    /**
     * Copies the constants defined in the <code>Constants</code> namespace into
     * the current <code>InfoStruct</code>.
     */
    void CopyConstants();
};

template <typename Writer>
void InfoStruct::Serialize(Writer& writer) const
{
    writer.StartObject();

    // Write the atomic parts array.
    writer.String("AtomicParts");
    writer.StartArray();
    {
        std::for_each(AtomicParts.begin(), AtomicParts.end(),
                [&writer](double v)
                {
                    writer.Double(v);
                });
    }
    writer.EndArray();

    // Write all of the other values.
    writer.String("ConeAngle")      , writer.Double(ConeAngle);
    writer.String("ConeTan2")       , writer.Double(ConeTan2);
    writer.String("ConeSolidAngle") , writer.Double(ConeSolidAngle);
    writer.String("Density")        , writer.Double(Density);
    writer.String("Frequency")      , writer.Double(Frequency);
    writer.String("Pd")             , writer.Double(Pd);
    writer.String("R0")             , writer.Double(R0);
    writer.String("Radi")           , writer.Double(Radi);
    writer.String("Rini")           , writer.Double(Rini);
    writer.String("Riso")           , writer.Double(Riso);
    writer.String("Temperature")    , writer.Double(Temperature);
    writer.String("Tini")           , writer.Double(Tini);
    writer.String("Vini")           , writer.Double(Vini);
    writer.String("WallTemp")       , writer.Double(WallTemp);

    writer.EndObject();
}

#endif /* __INFO_STRUCT_H__ */
