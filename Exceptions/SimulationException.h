#ifndef __SIMULATION_EXCEPTION_H__
#define __SIMULATION_EXCEPTION_H__

#include <stdexcept>

class SimulationException: public std::runtime_error
{
public:
    explicit SimulationException(const std::string& msg) : std::runtime_error(msg)
    { }
};

#endif
