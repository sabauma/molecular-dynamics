/**
 * This module defines the RPK class. This class is used for the integration
 * of the RPK equation. It defines the methods necessarry to step the equation
 * and determine its relevant state. In addition, it handles the computation
 * of collisions between the bubble wall and gas particles in the system.
 *
 * @author Spenser Bauman
 */

#ifndef __RPK_RPK_H__
#define __RPK_RPK_H__

#include "../Constants/Constants.h"
#include "../InfoStruct/InfoStruct.h"
#include "../Math/Functions.h"
#include "../Particle/Particle.h"
#include "../Units/Units.h"
#include "../Vector/Vector.hpp"
#include "PressureGauge.h"

#include "rapidjson/document.h"
#include "rapidjson/prettywriter.h"

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <math.h>
#include <stdio.h>

const int DIM = 2;

namespace odeint = boost::numeric::odeint;

typedef boost::numeric::ublas::vector<double> state_type;
typedef boost::numeric::ublas::matrix<double> matrix_type;
typedef odeint::runge_kutta_fehlberg78<state_type> error_stepper_type;
typedef odeint::controlled_runge_kutta<error_stepper_type> controlled_stepper_type;

// Numerical solution to Rayleigh-Plesset equation for bubble collapse
class RPK
{
private:

    const double& CurrentTime;

    double kR;       // Interpolation coefficient, aka velocity = R'

    state_type Time;
    state_type PrevRadius;
    state_type NextRadius;

    std::vector<double> Times;
    std::vector<double> Radii;
    std::vector<double> Velocities;

    controlled_stepper_type ControlledSolver;

    // Constants computed using the InfoStruct data. These are kept around for
    // performance reasons.
    double R03;
    double a3;
    double P0;
    double Density;
    double Viscosity;
    double SurfaceTension;
    double SpeedOfSound;
    double Omega;

    // Coefficients for the linear interpolation of the polytropic exponent.
    double A;
    double B;

    double TimeStep;

    // Pressure interpolation coefficients
    double ComputedPressure;    // Pressure used by the RPK function (ideal or pgauge)
    double ComputedDpdt;

    PressureGauge PGauge;

    InfoStruct Info;

    double StepSize;

    FILE* Actual;

    /**
     * Initialization function to setup the data using the passed InfoStruct
     * the object was constructed with.
     */
    void Init();

public:

    /**
     * Create a new RPK equation using the given <code>InfoStruct</code> for
     * the various constants needed. This assumes that the given info struct
     * is alread in simulation units.
     *
     * @param currentTime A reference to the variable holding the current
     *                    simulation time.
     * @param info An <code>InfoStruct</code> whose values are in simulation
     *             units.
     */
    RPK(const double& currentTime, const InfoStruct& info);

    /**
     * Clean up the dynamically allocated memory for the solver methods.
     */
    ~RPK();

    /**
     * @return The end time of the step made.
     */
    double NextStep();

    /**
     * Add a sample to the pressure gauge.
     */
    void AddSample(const double time, const double pressure);

    /**
     * @param time The time to compute the radius.
     * @return Bubble radius: R(Time)
     */
    double R(const double time) const;

    /**
     * @param The time to compute the radius.
     * @return Bubble velocity: V(t) = R'(t) = kR
     */
    double V(const double time) const;

    /**
     * Smarter version of the velocity calculation that uses the
     * computed velocities for each end of the current portion of the
     * RPK equation and interpolates over those, rather than guessing
     * by using the slope of the radius.
     */
    double VSmart(const double time) const;

    /**
     * @return Bubble Acceleration: R''(t)
     */
    double d2Rdt2() const;

    /**
     * @return Next time step
     */
    double NextStepTime();

    /**
     * Reads the pressure gauge to update the current pressure and dpdt value.
     * This should be called after a fixed amount of simulation time.
     */
    void ReadPressureGauge();

    /**
     * @param n The particle to compute intersection time with the bubble
     *          wall. NOTE: This function assumes that the <code>Update</code>
     *          method was performed used to update the particle position to
     *          where it should be at the current time.
     *
     * @return Time necessary for the given atom to impact the wall based on
     *         the estimated interpolation of the bubble wall position.
     */
    double FindCrossingTime(const Particle& n) const;

    /**
     * New serialization method that makes use of the ptree library.
     * Rather than write to a file, this function writes to the property
     * tree given as a parameter.
     */
    template <typename Writer>
    void Serialize(Writer& writer) const;

    /**
     * Deserialization method that can read in the data produced by the
     * serialize method.
     */
    // void Read(const boost::property_tree::ptree& tree);

    /**
     * Performs an evaluation of the RPK function. Given r, r', and t, this
     * will return the estimated next r' and r''.
     */
    state_type operator()(const state_type& x, state_type& dxdt, const double t) const;

    /**
     * Compute the acoustic pressure at the given time.
     */
    double AcousticPressure(const double time) const;

    /**
     * Compute the derivative of the acoustic pressure at the given time.
     */
    double AcousticPressureDerivative(const double time) const;
};

template <typename Writer>
void RPK::Serialize(Writer& writer) const
{

    writer.String("kR")             , writer.Double(kR);
    writer.String("Time")           , writer.Double(Time);
    writer.String("Time")           , writer.Double(Time);
    writer.String("TimeStep")       , writer.Double(TimeStep);
    writer.String("A")              , writer.Double(A);
    writer.String("B")              , writer.Double(B);
    writer.String("R03")            , writer.Double(R03);
    writer.String("a3")             , writer.Double(a3);
    writer.String("P0")             , writer.Double(P0);
    writer.String("Density")        , writer.Double(Density);
    writer.String("Viscosity")      , writer.Double(Viscosity);
    writer.String("SurfaceTension") , writer.Double(SurfaceTension);
    writer.String("SpeedOfSound")   , writer.Double(SpeedOfSound);
    writer.String("Omega")          , writer.Double(Omega);

    writer.String("PGauge");
    PGauge.Serialize<Writer>(writer);
}

#endif

