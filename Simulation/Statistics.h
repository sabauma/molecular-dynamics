#ifndef __SIMULATION_STATISTICS_H__
#define __SIMULATION_STATISTICS_H__

#include <algorithm>
#include <boost/array.hpp>
#include <boost/circular_buffer.hpp>
#include <boost/multi_array.hpp>

#include "./CollisionRecord.h"
#include "../Constants/AtomicProperties.h"
#include "../InfoStruct/InfoStruct.h"
#include "../Math/Functions.h"
#include "../Units/Units.h"
#include "../Vector/Vector.hpp"

class Statistics
{
public:

    typedef boost::multi_array<double, 2> Real2d;
    typedef boost::multi_array<double, 1> Real1d;
    typedef boost::multi_array<int, 2> Int2d;
    typedef boost::multi_array<int, 1> Int1d;

    typedef boost::multi_array<DoubleVector, 2> Vector2d;

    const InfoStruct& Info;
    const int Bins;
    const int CollisionNum;

    unsigned long long TotalCollisions;             /**< # of collision events           */
    unsigned long long TotalCellCrossings;          /**< # of cell crossing events       */
    unsigned long long TotalWallCollisions;         /**< # of wall collision events      */
    unsigned long long TotalConeBoundaryCollisions; /**< # of boundary collisions        */

    double FusionRate;          /**< # of fusions that have occurred */
    double MaxCollisionEnergy;

    /** Largest collision energy by the types of collision partners. */
    std::map<std::pair<int, int>, double> MaxCollisionEnergySpecies;

    // SNAPSHOT STATISTICS
    double KineticEnergyPerAtom;
    double CurrentRadius;

    boost::circular_buffer<CollisionRecord> Collisions; /**< Collision records.    */

    // Total atoms per binCurrentParticle = -1;
    Int1d TotalAtomsPerBin;

    // Atomic count per species per bin
    Int2d AtomsPerBin;
    // Atomic charge per species per bin
    Int2d ChargePerBin;
    // Velocity-squared per species per bin
    Real2d VelocitySquaredPerBin;
    // Radial velocity per species per bin
    Real2d RadialVelocityPerBin;
    // Largest temperature in each bin
    Real2d MaxTemperature;      /**< Max temp at each bin. */
    // The square of each particle's velocity corrected for the aggregate
    // motion of each bin.
    Real2d ThermalVelocitySquared;

    // Velocity of the bin's center of mass
    Vector2d VelocityCoM;

    Statistics(const InfoStruct& info, const int bins, const int collisions);

    /**
     * Giving a probably p, compute the quantile that corresponds to that
     * probability. That is to say, this returns the values x such that
     *
     * p(X <= x) = p
     *
     * In this case, it gives the collision distance which is larger than
     * p percent of the collision distances.
     */
    double CollisionDistanceQuantile(const double p);

    void RegisterCollision(CollisionRecord& rec);

    template <typename Iter>
    void ComputeStatistics(Iter begin, const Iter end, const double radius);

    template <typename Stream>
    void WriteCollisions(Stream& stream) const;

    template <typename Writer>
    void WriteStatistics(Writer& writer) const;

private:

    void Zero();

    template <typename Writer>
    void WriteTemperature(Writer& writer, const int bin) const;

    template <typename Writer>
    void WriteVelocity(Writer& writer, const int bin) const;

    template <typename Writer>
    void WriteIonization(Writer& writer, const int bin) const;

    template <typename Writer>
    void WritePercentage(Writer& writer, const int bin) const;

    template <typename Writer>
    void WriteMaxTemperature(Writer& writer, const int bin) const;

    template <typename Writer>
    void WriteKineticEnergy(Writer& writer, const int bin) const;

    template <typename Writer>
    void WriteEnergyCoM(Writer& writer, const int bin) const;
};

template <typename Iter>
void Statistics::ComputeStatistics(
        Iter begin,
        const Iter end,
        const double radius)
{
    this->Zero();

    CurrentRadius = radius;

    for (Iter particle = begin; particle != end; ++particle)
    {
        // Atom Type
        const int type  = particle->Type;
        const double v2 = dot(particle->Velocity, particle->Velocity);

        // Total kinetic energy: 0.5*m*v*v
        KineticEnergyPerAtom += 0.5 * Info.AtomicMass[type] * v2;

        // Distance from bubble center
        const double r = norm2(particle->Position);

        // Determine bin # to which this atom belongs
        int bin = int(r * (double) Bins / radius);
        bin = std::min(bin, Bins - 1);

        assert(type < Constants::ELEMENT_TYPES);
        assert(type >= 0);

        // Increase atom count in the bin
        AtomsPerBin[type][bin]++;

        // Increase total count of atoms
        TotalAtomsPerBin[bin]++;

        // Increase charge per bin
        ChargePerBin[type][bin] += particle->Charge;

        // Increase SUM(v*v) per bin
        VelocitySquaredPerBin[type][bin] += v2;

        // Increase sum of radial velocities per bin
        RadialVelocityPerBin[type][bin] +=
            dot(particle->Position, particle->Velocity) / r;

        VelocityCoM[type][bin] += particle->Velocity;
    }

    for (int bin = 0; bin < Bins; ++bin)
    {
        for (int type = 0; type < Constants::ELEMENT_TYPES; ++type)
        {
            VelocityCoM[type][bin] /= AtomsPerBin[type][bin];
        }
    }

    for (Iter particle = begin; particle != end; ++particle)
    {
        const int type  = particle->Type;

        // Distance from bubble center
        const double r = norm2(particle->Position);

        // Determine bin # to which this atom belongs
        int bin = int(r * (double) Bins / radius);
        bin = std::min(bin, Bins - 1);

        const DoubleVector vCoM(particle->Velocity - VelocityCoM[type][bin]);
        ThermalVelocitySquared[type][bin] += dot(vCoM, vCoM);
    }

    for (int bin = 0; bin < Bins; ++bin)
    {
        for (int type = 0; type < Constants::ELEMENT_TYPES; ++type)
        {
            if (AtomsPerBin[type][bin] > 0)
            {
                const double radial =
                    RadialVelocityPerBin[type][bin] / AtomsPerBin[type][bin];
                const double squared =
                    VelocitySquaredPerBin[type][bin] / AtomsPerBin[type][bin];

                // Calculate species temperature
                const double temperature =
                    Info.AtomicMass[type] / 3.0 * (squared - Sqr(radial));

                MaxTemperature[type][bin] =
                    std::max(MaxTemperature[type][bin], temperature);
            }
        }
    }
}

template <typename Stream>
void Statistics::WriteCollisions(Stream& stream) const
{
    for (auto it = Collisions.begin(); it != Collisions.end(); ++it)
    {
        stream << it->Time * Units::T << ','
               << it->Position[0] * Units::L << ','
               << it->Position[1] * Units::L << ','
               << it->Position[2] * Units::L << ','
               << it->DeltaV * Units::L / Units::T << ','
               << it->Energy * Units::Energy  << ','
               << it->Distance * Units::L << ','
               << it->Type1 << ','
               << it->Type2 << '\n';
    }
}

template <typename Writer>
void Statistics::WriteTemperature(Writer& writer, const int bin) const
{
    writer.String("temperature");
    writer.StartObject();
    for (int type = 0; type < Constants::ELEMENT_TYPES; ++type)
    {
        const std::string& type_name = Constants::ELEMENT_NAMES[type];

        writer.String(type_name.c_str());
        if (AtomsPerBin[type][bin] > 0)
        {
            // Calculate species temperature
            const double temperature =
                Units::Temperature * Info.AtomicMass[type] / 3.0
                * (VelocitySquaredPerBin[type][bin] / AtomsPerBin[type][bin]
                   - Sqr(RadialVelocityPerBin[type][bin] / AtomsPerBin[type][bin]));

            writer.Double(temperature);
        }
        else
        {
            writer.Double(0.0);
        }
    }
    writer.EndObject(); // End current temperature object
}

template <typename Writer>
void Statistics::WriteVelocity(Writer& writer, const int bin) const
{
    writer.String("velocity");
    writer.StartObject();
    for (int type = 0; type < Constants::ELEMENT_TYPES; ++type)
    {
        const std::string& type_name = Constants::ELEMENT_NAMES[type];

        writer.String(type_name.c_str());
        if (AtomsPerBin[type][bin] > 0)
        {
            const double velocity =
                RadialVelocityPerBin[type][bin] / AtomsPerBin[type][bin];
            writer.Double(velocity * Units::L / Units::T);
        }
        else
        {
            writer.Double(0.0);
        }
    }
    writer.EndObject(); // End current velocity object
}

template <typename Writer>
void Statistics::WriteIonization(Writer& writer, const int bin) const
{
    writer.String("ionization");
    writer.StartObject();
    for (int type = 0; type < Constants::ELEMENT_TYPES; ++type)
    {
        const std::string& type_name = Constants::ELEMENT_NAMES[type];

        writer.String(type_name.c_str());
        if (AtomsPerBin[type][bin] > 0)
        {
            writer.Double(double(ChargePerBin[type][bin]) / AtomsPerBin[type][bin]);
        }
        else
        {
            writer.Double(0.0);
        }
    }
    writer.EndObject(); // End current ionization object
}

template <typename Writer>
void Statistics::WritePercentage(Writer& writer, const int bin) const
{
    writer.String("percentage");
    writer.StartObject();
    for (int type = 0; type < Constants::ELEMENT_TYPES; ++type)
    {
        const std::string& type_name = Constants::ELEMENT_NAMES[type];

        writer.String(type_name.c_str());
        if (AtomsPerBin[type][bin] > 0)
        {
            writer.Double((double) AtomsPerBin[type][bin] / (double) TotalAtomsPerBin[bin]);
        }
        else
        {
            writer.Double(0.0);
        }
    }
    writer.EndObject(); // End current ionization object
}

template <typename Writer>
void Statistics::WriteMaxTemperature(Writer& writer, const int bin) const
{
    writer.String("max_temperature");
    writer.StartObject();
    for (int type = 0; type < Constants::ELEMENT_TYPES; ++type)
    {
        const std::string& type_name = Constants::ELEMENT_NAMES[type];

        writer.String(type_name.c_str());
        writer.Double(MaxTemperature[type][bin] * Units::Temperature);
    }
    writer.EndObject(); // End current max temperature object
}

template <typename Writer>
void Statistics::WriteKineticEnergy(Writer& writer, const int bin) const
{
    writer.String("kinetic_energy");
    writer.StartObject();
    for (int type = 0; type < Constants::ELEMENT_TYPES; ++type)
    {
        const std::string& type_name = Constants::ELEMENT_NAMES[type];

        writer.String(type_name.c_str());
        if (AtomsPerBin[type][bin] > 0)
        {
            const double ke =
                0.5 * Info.AtomicMass[type] * VelocitySquaredPerBin[type][bin];
            writer.Double(ke * Units::Energy);
        }
        else
        {
            writer.Double(0.0);
        }
    }
    writer.EndObject(); // End current temperature object
}

template <typename Writer>
void Statistics::WriteEnergyCoM(Writer& writer, const int bin) const
{
    writer.String("kinetic_energy_com");
    writer.StartObject();

    for (int type = 0; type < Constants::ELEMENT_TYPES; ++type)
    {
        const std::string& type_name = Constants::ELEMENT_NAMES[type];

        writer.String(type_name.c_str());
        if (AtomsPerBin[type][bin] > 0)
        {
            const double thermal =
                0.5 * Info.AtomicMass[type] * ThermalVelocitySquared[type][bin];
            writer.Double(thermal * Units::Energy);
        }
        else
        {
            writer.Double(0.0);
        }
    }
    writer.EndObject(); // End current temperature object
}

template <typename Writer>
void Statistics::WriteStatistics(Writer& writer) const
{
    writer.StartArray();
    for (int bin = 0; bin < Bins; ++bin)
    {
        // Bin volume (from cone volume)
        const double volume = Info.ConeSolidAngle / 3.0
                            * Pow3(CurrentRadius / Bins * Units::L)
                            * (Pow3(bin + 1) - Pow3(bin));

        const double density = TotalAtomsPerBin[bin] / volume;

        writer.StartObject();

        const double radius = (bin + 0.5) * CurrentRadius * Units::L / Bins;
        writer.String("r"), writer.Double(radius);
        writer.String("n"), writer.Double(TotalAtomsPerBin[bin]);
        writer.String("density"), writer.Double(density);

        // Write temperature information for each element
        this->WriteTemperature(writer, bin);
        this->WriteVelocity(writer, bin);
        this->WriteIonization(writer, bin);
        this->WritePercentage(writer, bin);
        this->WriteMaxTemperature(writer, bin);
        this->WriteKineticEnergy(writer, bin);
        this->WriteEnergyCoM(writer, bin);

        writer.EndObject(); // End the current bin object.
    }

    writer.EndArray();

    writer.String("Fusionrate"); writer.Uint64(FusionRate);
    writer.String("MaxCollisionEnergy");
    writer.Double(MaxCollisionEnergy * Units::M * Sqr(Units::L / Units::T));

    writer.String("CollsionEnergies");
    writer.StartObject();
    for (auto it = MaxCollisionEnergySpecies.begin();
            it != MaxCollisionEnergySpecies.end(); ++ it)
    {
        const std::pair<int, int>& key = it->first;
        const double energy            = it->second;

        std::string key_name = Constants::ELEMENT_NAMES[key.first] + "-"
                             + Constants::ELEMENT_NAMES[key.second];

        writer.String(key_name.c_str());
        writer.Double(energy * Units::Energy);
    }
    writer.EndObject();
}

#endif
