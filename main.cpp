
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <vector>

#include "Constants/AtomicProperties.h"
#include "Constants/Constants.h"
#include "Debug/Trace.h"
#include "EventCalendar/EventType.h"
#include "InfoStruct/InfoStruct.h"
#include "Particle/Particle.h"
#include "Simulation/Simulation.h"
#include "Units/Units.h"
#include "Vector/Vector.hpp"

//static const std::string DELETE_LINE(240, '\b');
static const std::string DELETE_LINE = "\r";

int main(int argc, char** argv)
{
    std::string filename;

    if (argc < 2)
    {
        std::cout << "No filename given." << std::endl;
        exit(EXIT_FAILURE);
    }
    else
    {
        filename = argv[1];
    }

    InfoStruct is =
        InfoStruct::ParseInfoStruct(filename).ToSimulationUnits();

    std::cout << filename << std::endl;
    Simulation simulation(is);

    printf("Initilization Information\n\n");

    std::cout << is << std::endl;

    printf("Units::L = %g m\n", Units::L);
    printf("Units::T = %g s\n", Units::T);
    printf("Units::M = %g kg\n", Units::M);

    printf("Units::Pressure = %g Pa\n", Units::Pressure);
    printf("Units::DynamicViscosity = %g Pa*s\n", Units::DynamicViscosity);
    printf("Units::SurfaceTension = %g ???\n", Units::SurfaceTension);

    printf("AtomicDiameterSquared\n");
    for (int i = 0; i < Constants::ELEMENT_TYPES; ++i)
    {
        for (int j = 0; j < Constants::ELEMENT_TYPES; ++j)
        {
            printf("%12g", is.AtomicDiameterSquared[i][j]);
        }
        printf("\n");
    }

    printf("BirdConstant\n");
    for (int j = 0; j < Constants::ELEMENT_TYPES; ++j)
    {
        printf("%g\n", is.BirdConstant[j]);
    }

    printf("AtomicMass\n");
    for (int i = 0; i < Constants::ELEMENT_TYPES; ++i)
    {
        printf("%12g\n", is.AtomicMass[i]);
    }

    printf("ReducedMass\n");
    for (int i = 0; i < Constants::ELEMENT_TYPES; ++i)
    {
        for (int j = 0; j < Constants::ELEMENT_TYPES; ++j)
        {
            printf("%12g", is.ReducedMass[i][j]);
        }
        printf("\n");
    }

    printf("VWallThermal\n");
    for (int i = 0; i < Constants::ELEMENT_TYPES; ++i)
    {
        printf("%12g\n", is.VWallThermal[i]);
    }

    int after_rebound = 0;
    for (unsigned long long i = 0, updates = 0; ; ++i)
    {
        EventType type = simulation.NextStep();
        //std::cerr << getEventName(type) << std::endl;

        if (i % 1024000 == 0)
        {
            std::cout << DELETE_LINE;
            printf("Processing Event %llu at %0.12e of type %s",
                   i, simulation.GetCurrentTime(), getEventName(type).c_str());
            fflush(stdout);
        }

        if (type == UpdateSystemEvent)
        {
            char fname[256] = {'\0'};
            sprintf(fname, "./SavePoint%.6llu.txt", updates + 1);
            std::ofstream particle_file(fname);

            std::vector<Particle>::const_iterator it = simulation.begin();

            for (; it != simulation.end(); ++it)
            {
                particle_file << it->Position[0] << ','
                              << it->Position[1] << ','
                              << it->Position[2] << ','
                              << it->Velocity[0] << ','
                              << it->Velocity[1] << ','
                              << it->Velocity[2] << ','
                              << it->Type        << '\n';
            }

            particle_file.close();
            ++updates;

            if (simulation.HasRebounded() && ++after_rebound == 10) break;
        }
    }

    return 0;
}

