
#include "./InfoStruct/InfoStruct.h"
#include "./Math/Functions.h"
#include "./Particle/Particle.h"
#include "./Units/Units.h"
#include "./Vector/Vector.hpp"
#include "./Random/Random.h"
#include "./Simulation/Statistics.h"

#include <iostream>
#include <boost/math/special_functions/fpclassify.hpp>

static const double gas_temperature    = 300.0;
static const double liquid_temperature = 300.0;

static const int type = 0;

static const DoubleVector normal(0.0, 0.0, -1.0);

int main(int argc, char** argv)
{
    if (argc <= 1)
    {
        std::cout << "Need an infostruct to execute." << std::endl;
    }

    std::vector<Particle> particles;

    DoubleVector p_wall(0.0, 0.0, 0.0);
    DoubleVector v_wall(0.0, 0.0, 0.0);

    Particle pa;
    pa.Position = DoubleVector(0.0, 0.0, 1.0);
    pa.Type = type;

    InfoStruct Info = InfoStruct::ParseInfoStruct(argv[1]).ToSimulationUnits();

    const double offset =
        0.5 * (Info.AtomicDiameter[type] + Info.LiquidAtomicDiameter);

    p_wall = pa.Position - normal * offset;


    //Find a velocity vector for the wall particle that will result in it
    //colliding with the gas particle.
    for (int i = 0, iters = 0; i < 10000; ++i, iters = 0)
    {
        double int_time = 0.0;
        maxwell_boltzman_velocity(gas_temperature / Units::Temperature, Info.AtomicMass[type], pa.Velocity);
        do
        {
            maxwell_boltzman_velocity(liquid_temperature / Units::Temperature, Info.LiquidAtomicMass, v_wall);

            int_time = intersection_time(pa.Position, pa.Velocity, p_wall, v_wall, offset);
        }
        while (isinf(int_time) && ++iters < 100000);

        if (!isinf(int_time))
        {
            process_collision(pa.Position, pa.Velocity, Info.AtomicMass[type],
                              p_wall, v_wall, Info.LiquidAtomicMass);
            particles.push_back(pa);
        }

    }

    Statistics Stats(Info, 1, 100);
    Stats.ComputeStatistics(particles.begin(), particles.end(), 1.0);

    const double temp
        = Info.AtomicMass[type] / 3.0
        * (Stats.VelocitySquaredPerBin[type][0] / Stats.AtomsPerBin[type][0]);

    std::cout << "Wall Temperature: " << temp * Units::Temperature << std::endl;
    std::cout << "Radial Velocity: " << Stats.RadialVelocityPerBin[type][0] * Units::L / Units::T << std::endl;

}
