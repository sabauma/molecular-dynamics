
#include <vector>
#include <boost/multi_array.hpp>

#include "./Statistics.h"

Statistics::Statistics(const InfoStruct& info, const int bins, const int collisions) :
    Info(info),
    Bins(bins),
    CollisionNum(collisions),
    TotalCollisions(0),
    TotalCellCrossings(0),
    TotalWallCollisions(0),
    TotalConeBoundaryCollisions(0),
    FusionRate(0),
    MaxCollisionEnergy(0),
    Collisions(collisions),
    TotalAtomsPerBin(boost::extents[bins]),
    AtomsPerBin(boost::extents[Constants::ELEMENT_TYPES][bins]),
    ChargePerBin(boost::extents[Constants::ELEMENT_TYPES][bins]),
    VelocitySquaredPerBin(boost::extents[Constants::ELEMENT_TYPES][bins]),
    RadialVelocityPerBin(boost::extents[Constants::ELEMENT_TYPES][bins]),
    MaxTemperature(boost::extents[Constants::ELEMENT_TYPES][bins]),
    ThermalVelocitySquared(boost::extents[Constants::ELEMENT_TYPES][bins]),
    VelocityCoM(boost::extents[Constants::ELEMENT_TYPES][bins])
{
    // We do not want to zero out the max temperature after every snapshot,
    // since it is the max of all snapshots.
    std::fill(
            MaxTemperature.data(),
            MaxTemperature.data() + MaxTemperature.num_elements(),
            0);

    this->Zero();
}

void Statistics::Zero()
{
    KineticEnergyPerAtom = 0;
    CurrentRadius = 0.0;

    std::fill(TotalAtomsPerBin.begin(), TotalAtomsPerBin.end(), 0);

    std::fill(
            AtomsPerBin.data(),
            AtomsPerBin.data() + AtomsPerBin.num_elements(),
            0);

    std::fill(
            ChargePerBin.data(),
            ChargePerBin.data() + AtomsPerBin.num_elements(),
            0);

    std::fill(
            VelocitySquaredPerBin.data(),
            VelocitySquaredPerBin.data() + VelocitySquaredPerBin.num_elements(),
            0.0);

    std::fill(
            ThermalVelocitySquared.data(),
            ThermalVelocitySquared.data() + ThermalVelocitySquared.num_elements(),
            0.0);

    std::fill(
            RadialVelocityPerBin.data(),
            RadialVelocityPerBin.data() + RadialVelocityPerBin.num_elements(),
            0.0);

    std::fill(
            VelocityCoM.data(),
            VelocityCoM.data() + VelocityCoM.num_elements(),
            0.0);

    std::fill(
            ThermalVelocitySquared.data(),
            ThermalVelocitySquared.data() + ThermalVelocitySquared.num_elements(),
            0.0);
}

double Statistics::CollisionDistanceQuantile(const double p)
{
    const size_t size  = Collisions.size();
    const size_t quant = (size_t) ((double) size * p);
    std::vector<double> dists(size, 0.0);

    for (size_t i = 0; i < Collisions.size(); ++i)
    {
        dists[i] = Collisions[i].Distance;
    }

    std::nth_element(dists.begin(), dists.begin() + quant, dists.end());

    return dists[quant];
}

void Statistics::RegisterCollision(CollisionRecord& rec)
{
    // Make sure type1 is the smaller element. Just for convenience.
    if (rec.Type1 > rec.Type2)
    {
        std::swap(rec.Type1, rec.Type2);
    }

    if (rec.Energy >= Info.CollisionThreshold)
    {
        this->Collisions.push_back(rec);
    }

    const std::pair<int, int> key = std::make_pair(rec.Type1, rec.Type2);
    MaxCollisionEnergySpecies[key] =
        std::max(MaxCollisionEnergySpecies[key], rec.Energy);

    MaxCollisionEnergy = std::max(MaxCollisionEnergy, rec.Energy);
}
