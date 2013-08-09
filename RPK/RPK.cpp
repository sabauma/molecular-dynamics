
#include "../Exceptions/SimulationException.h"
#include "RPK.h"

#include <boost/numeric/odeint.hpp>
#include <boost/ref.hpp>
#include <limits>

using namespace boost::numeric::odeint;

static const unsigned int MAX_ITER = 100000;

static const double REL_ERROR = 1.0e-15;
static const double ABS_ERROR = 1.0e-15;
static const double A_X       = 1.0;
static const double A_DXDT    = 1.0;

/**
 * Performs a linear interpolation of the bubble radius to determine the
 * radius when a pressure sample was taken.
 *
 * @param times The times of each step in the RPK equation.
 * @param radii The radii of each step in the RPK equation.
 * @param sample_time The time of a pressure sample for which the radius is
 *                    computed.
 * @return The radius at time sample_time.
 */
static double
compute_radius(const std::vector<double>& times,
               const std::vector<double>& radii,
               const double sample_time)
{
    assert(!times.empty() && !radii.empty());
    assert(times.size() == radii.size());

    size_t index  = times.size() - 1;

    // Find the first RPK step before the pressure sample.
    for (; times[index] > sample_time && index > 0; --index) ;

    assert(times[index] <= sample_time && times[index+1] >= sample_time);

    const double dr    = radii[index + 1] - radii[index];
    const double dt    = times[index + 1] - times[index];
    const double slope = dr / dt;

    double interpolated = radii[index] + slope * (sample_time - times[index]);

    // Performa  linear interpolation of the radii to get a close match.
    return interpolated;
}

static inline bool isvalid(const double x)
{
    return !isnan(x) && !isinf(x);
}

RPK::RPK(const double& currentTime, const InfoStruct& info):
    CurrentTime(currentTime), Time(2), PrevRadius(2), NextRadius(2),
    ControlledSolver(default_error_checker<double>(ABS_ERROR, REL_ERROR, A_X, A_DXDT)),
    Info(info)
{
	Actual = fopen("RPactual.dat", "w");

    this->Init();

    StepSize = 1.0e-3 / Info.Frequency;
}

RPK::~RPK()
{
    fclose(Actual);
}

void RPK::AddSample(const double time, const double pressure)
{
#ifdef USE_PGAUGE
    this->PGauge.AddSample(time, pressure);
#endif
}

// Bubble radius: R(Time)
double RPK::R(const double time) const
{
    // R(t) = R0 + kR*(t - t0)
    return kR * (time - Time[0]) + PrevRadius[0];
}

// Bubble velocity: V(t) = R'(t) = kR
double RPK::V(const double /* time */) const
{
    return kR;
}

double RPK::VSmart(const double time) const
{
    return PrevRadius[1] + d2Rdt2() * (time - Time[0]);
}

// Bubble Acceleration: R''(t)
double RPK::d2Rdt2() const
{
    // R''(t) = (v1 - v0) / (t1 - t0)
    return (NextRadius[1] - PrevRadius[1]) / (Time[1] - Time[0]);
}

// Next time step
double RPK::NextStepTime()
{
    return Time[1];
}

// Determine next time step in RP numerical solution
double RPK::NextStep()
{

#ifdef USE_PGAUGE

    if (PGauge.CanSample())
    {
        const double sample_time = PGauge.PressureTime();
        const double radius      = compute_radius(Times, Radii, sample_time);

        // PGauge.GetPressure(this->CurrentTime, CurrentPressure, dpdt);
        PGauge.SampleGauge(radius);
    }

    const bool use_pgauge = PGauge.IsReady();

#else
    const bool use_pgauge = false;
#endif

	// Advance initial conditions to the next time step
	Time[0]       = Time[1];
    PrevRadius[0] = NextRadius[0];
    PrevRadius[1] = NextRadius[1];

    controlled_step_result status = fail;

    state_type state_backup(2);
    double time_backup;
    double step_backup;

    while (status == fail)
    {
        std::copy(NextRadius.begin(), NextRadius.end(), state_backup.begin());
        time_backup = Time[1];
        step_backup = StepSize;

        status = ControlledSolver.try_step(
                boost::ref(*this), NextRadius, Time[1], StepSize);

        // If an error results in nan, we took too large of a step. Restore
        // the state and reduce step size by a factor of 0.5;
        if (!isvalid(NextRadius[0]) || !isvalid(NextRadius[1]) || !isvalid(Time[1]))
        {
            status   = fail;
            Time[1]  = time_backup;
            StepSize = step_backup * 0.95;

            std::copy(state_backup.begin(), state_backup.end(), NextRadius.begin());
        }
    }

    // Store the time, radius, and velocity of each step of the RPK equation.
    Times.push_back(Time[1]);
    Radii.push_back(NextRadius[0]);
    Velocities.push_back(NextRadius[1]);

    printf("StepSize: %e\n", StepSize);

    if (status != success || !isvalid(NextRadius[0]) || !isvalid(NextRadius[1]))
    {
        fprintf(stderr, "ERROR: Non-stiff Problem at x = %25.17e ~ %d\n",
                Time[1], (int) status);
        throw SimulationException("Error Stepping the RPK solution.");
    }

	// Determine interpolation coefficients. See R(t) and V(t) for explanation
    if (NextRadius[0] != PrevRadius[0] && Time[0] != Time[1])
    {
        kR = (NextRadius[0] - PrevRadius[0]) / (Time[1] - Time[0]);
    }

    double p    = 0.0;
    double dpdt = 0.0;

    if (use_pgauge)
    {
        p    = PGauge.ComputePressure(Time[1], PrevRadius[0], PrevRadius[1]);
        dpdt = PGauge.ComputeDpdt(Time[1], PrevRadius[0], PrevRadius[1]);
    }

	// Write actual R(t) and V(t)
	fprintf(Actual, "%.15e,\t%.15e,\t%.15e,\t%.15e,\t%.15e\n",
            Time[1] * Units::T,
            NextRadius[0] * Units::L,
            NextRadius[1] * Units::L/Units::T,
            p * Units::Pressure,
            dpdt * Units::Pressure / Units::T);
    fflush(Actual);

	return Time[1];
}

state_type RPK::operator()(const state_type& y, state_type& dydt, const double time) const
{
    dydt[0] = y[1];

	// R^2
	const double R2 = y[0] * y[0];
	// R^3
	const double R3 = R2 * y[0];

	// Gas pressure within the bubble.
	double Pg;
    // Derivative of the gas pressure.
    double Pg_dot;

	// PGauge is working?
    if (this->PGauge.IsReady())
	{
		// Use both pressure and derivative dpdt: linear interpolation
        Pg     = this->PGauge.ComputePressure(time, y[0], y[1]);
        Pg_dot = this->PGauge.ComputeDpdt(time, y[0], y[1]);
	}
	else
	{
        const double surface   = 2.0 * SurfaceTension / Info.R0;
        const double bubble_P0 = this->P0 + surface;

		// Adiabatic index
		// Isothermal expansion or adiabatic? Interpolate for in-between
        double k = y[0] > Info.Riso ? 1
                 : (y[0] > Info.Radi
                         ? (this->B + this->A * y[0])
                         : Constants::GAMMA);

        Pg     = bubble_P0 * pow(this->R03 / (R3 - a3), k);
        Pg_dot = -3.0 * k * Pg * R2 * y[1] / R3;
	}



#ifdef USE_RP

	// Acoustic driving pressure and driving pressure derivative
	const double Pa     = this->AcousticPressure(time);
	const double Pa_dot = this->AcousticPressureDerivative(time);

    // Note that we use the kinematic viscosity here, so the viscosity term
    // lacks an accompanying division by density as most forumulations give,
    // as kinematic viscosity is the dynamic viscosity divided by density.
	const double t1 = (Pg - this->P0 - Pa) / this->Density;
	const double t2 = -4.0 * this->Viscosity * y[1] / y[0];
	const double t3 = -2.0 * this->SurfaceTension / this->Density / y[0];
	const double t4 = y[0] * (Pg_dot - Pa_dot) / this->Density / this->SpeedOfSound;
	const double t5 = -1.5 * y[1] * y[1];

	dydt[1] = (t1 + t2 + t3 + t4 + t5) / y[0];

#elif NAIVE_RP

	const double Pa     = this->AcousticPressure(time);

    // A naive formulation of the RP equation which lacks terms for viscosity,
    // surface tension, and acoustic damping.
    dydt[1] = ((Pg - this->P0 - Pa) / this->Density - 3.0 / 2.0 * y[1] * y[1])
            / y[0];

#else

    // Need to use the dynamic viscosity, but we have kinematic viscosity, so
    // multiply by the density.
    const double VIS = this->Viscosity * this->Density;

    const double tr  = y[0] / this->SpeedOfSound;

	// Acoustic driving pressure and driving pressure derivative. Here, the
    // acoustic pressure is offset by the traversal time 'tr', but the
    // derivative is at the current time (this is not a mistake).
	const double Pa = this->AcousticPressure(time + tr);

    // Mach number of the bubble wall at the current time. Used to correct for
    // damping at high velocities.
    const double Mach   = y[1] / this->SpeedOfSound;

    const double Pb     = Pg - 2.0 * this->SurfaceTension / y[0]
                        - 4.0 * VIS * y[1] / y[0];
    const double denom  = (1.0 - Mach) * y[0] * this->Density
                        + 4.0 * VIS * tr / y[0];

    dydt[1] = ((1 + Mach) * (Pb - Pa - this->P0)
               - 3.0 / 2.0 * (1 - Mach / 3.0) * y[1] * y[1] * this->Density
               + tr * Pg_dot
               + 2.0 * this->SurfaceTension * tr / R2
               + 4.0 * VIS * y[1] * y[1] * tr / R2) / denom;

#endif /* ifdef USE_RP */

    state_type retval(2);
    retval[0] = Pg;
    retval[1] = Pg_dot;

    return retval;
}

double RPK::AcousticPressure(const double time) const
{
    switch (Info.Driver)
    {
        case Acoustic:
            return -Info.Pd * sin(this->Omega * (time + Info.Tini));
        case Constant:
            return Info.Pd;
    }

    throw SimulationException("RPK::AcousticPressure: Unknown driver mode.");
}

double RPK::AcousticPressureDerivative(const double time) const
{
    switch (Info.Driver)
    {
        case Acoustic:
            return -Info.Pd * this->Omega * cos(this->Omega * (time + Info.Tini));
        case Constant:
            return 0.0;
    }

    throw SimulationException("RPK::AcousticPressureDerivative: Unknown driver mode");
}

double RPK::FindCrossingTime(const Particle& n) const
{
	/*
     * Here we solve quadratic equation (r + vt)^2 = R(t)^2, where r and v are
     * position vectors and R(t) is the bubble's radius R(t) is linearly
     * interpolated between adjacent data points scalar products (r,r), (r,v)
     * and (v,v), which were computed at time sp_time[n], are stored in xx[n],
     * xv[n] and vv[n], correspondingly.
     */

    const double sp_time = n.Updated;

    // R at the time the scalar products were computed
	const double sp_R = PrevRadius[0] + kR * (sp_time - Time[0]);

    // Coefficients of the quadratic equation to solve.
    const double a = dot(n.Velocity, n.Velocity) - Sqr(kR);
    const double b = dot(n.Position, n.Velocity) - sp_R * kR;
    const double c = dot(n.Position, n.Position) - Sqr(sp_R);

	// Infinite time by default
	double time = std::numeric_limits<double>::infinity();

    // There is no solution to the
	/*if (c > 0 && b >= 0) {
		printf("A terrible thing has happened!\n");
	}*/

	if (a > 0)
    {
        if (b > 0)
        {
            time = -c / (b + sqrt(b*b - a*c));
        }
        else
        {
            time = (-b + sqrt(b*b - a*c)) / a;
        }
    }
	else
    {
        if (b > 0)
        {
            time = -c / (b + sqrt(b*b - a*c));
        }
    }

	time += sp_time;

	return (time <= Time[1]) ? time : std::numeric_limits<double>::infinity();
}

/*
 *void RPK::Read(const ptree& tree)
 *{
 *    kR              = tree.get<double>("kR");
 *    Time[0]         = tree.get<double>("Time_0");
 *    Time[1]         = tree.get<double>("Time_1");
 *    TimeStep        = tree.get<double>("TimeStep");
 *    A               = tree.get<double>("A");
 *    B               = tree.get<double>("B");
 *    R03             = tree.get<double>("R03");
 *    a3              = tree.get<double>("a3");
 *    P0              = tree.get<double>("P0");
 *    Density         = tree.get<double>("Density");
 *    Viscosity       = tree.get<double>("Viscosity");
 *    SurfaceTension  = tree.get<double>("SurfaceTension");
 *    SpeedOfSound    = tree.get<double>("SpeedOfSound");
 *    Omega           = tree.get<double>("Omega");
 *
 *    PGauge.Read(tree.get_child("PGauge"));
 *}
 */

void RPK::Init()
{
    // Density of the carrier liquid, kg/m^3
    // Density = Constants::Density * Pow3(Units::L) / Units::M;
    Density        = Info.LiquidDensity;
    Viscosity      = Info.LiquidViscosity;
    SpeedOfSound   = Info.LiquidSpeedOfSound;
    SurfaceTension = Info.LiquidSurfaceTension;

    // Viscosity, Pa*s
    // Viscosity = Constants::VIS / Units::Viscosity;

    // Speed of sound in liquid, m/s
    // SpeedOfSound = Constants::C * Units::T / Units::L;

    // Acoustic drive circular Frequency, 1/s
    Omega = 2.0 * Constants::PI * Info.Frequency;

    // R03 = R0^3
    R03 = Pow3(Info.R0);

    // Surface tension, N/m
    // SurfaceTension = Constants::SurfaceTension / Units::SurfaceTension;

    // Ambient pressure
    // P0 = Constants::P0 / Units::Pressure;
    P0 = Info.AmbientPressure;

    const double surface = 2.0 * SurfaceTension / Info.R0;

    // Van der Waals total excluded volume
    a3 = R03 * (P0 + surface) * Constants::ExcludedVolume
        / (Constants::AVOGADRO * Info.Temperature * Pow3(Units::L));

    // Coefficients for the linear interpolation of the polytropic exponent
    A = (1.0 - Constants::GAMMA) / (Info.Riso - Info.Radi);
    B = (Constants::GAMMA * Info.Riso - Info.Radi) / (Info.Riso - Info.Radi);

    const int DIM = 2;

    double x_end = 1.0 / Info.Frequency;

    /* Variables used for evolving */

    double h_nonstiff = 1.0e-9 / Info.Frequency;
    double x_nonstiff = 0.0;

    state_type y_nonstiff(DIM);
    y_nonstiff[0] = Info.Rini;
    y_nonstiff[1] = Info.Vini;

    /* Definitions to determine ODE solvers */
    controlled_stepper_type solver_nonstiff(
            default_error_checker<double>(1.0e-12, 1.0e-12, 1.0, 1.0));

    FILE* out = fopen("RPguess2.dat", "w");
    /* Integration with stepsize control */
    fprintf(out, "Non-stiff Problem:\n");

    fprintf(out, "     x                    R[0]                     V[1]       P[2]\n");

    while(x_nonstiff < x_end)
    {
        controlled_step_result status = fail;

        for (unsigned int i = 0; i < MAX_ITER && status == fail; ++i)
        {
            state_type y_backup(2);
            y_backup[0] = y_nonstiff[0];
            y_backup[1] = y_nonstiff[1];
            double x_backup = x_nonstiff;

            status = solver_nonstiff.try_step(
                    boost::ref(*this), y_nonstiff, x_nonstiff, h_nonstiff);

            if (!isvalid(y_nonstiff[0]) || !isvalid(y_nonstiff[1]))
            {
                status        = fail;
                y_nonstiff[0] = y_backup[0];
                y_nonstiff[1] = y_backup[1];
                x_nonstiff    = x_backup;
                h_nonstiff   /= 2.0;
            }
        }

        state_type dydx(2);
        state_type pressure = (*this)(y_nonstiff, dydx, x_nonstiff);

        // true_solution(true_y, x_nonstiff);
        // vec_rel_error(relerr_y, y_nonstiff, true_y, DIM);
        fprintf(out, "%15.7e, %25.17e, %25.17e, %25.17e, %25.17e\n",
                x_nonstiff * Units::T,
                y_nonstiff[0] * Units::L,
                y_nonstiff[1] * Units::L / Units::T,
                pressure[0] * Units::Pressure,
                pressure[1] * Units::Pressure / Units::T);
    }

    fclose(out);

    // Reset initial conditions
    Time[0]       = Time[1]       = 0.0;
    PrevRadius[0] = NextRadius[0] = Info.Rini;
    PrevRadius[1] = NextRadius[1] = Info.Vini;

    Times.push_back(Time[0]);
    Radii.push_back(NextRadius[0]);

    // Write header for actual RP solution
    // out = fopen("RPactual.dat", "w");
    fprintf(Actual, "Time\tR\tV\tPressureNow\tdpdt\n");
    fprintf(Actual, "%.15e,\t%.15e,\t%.15e,\t%.15e,\t%.15e\n",
            Time[0] * Units::T,
            PrevRadius[0] * Units::L,
            PrevRadius[1] * Units::L/Units::T,
            0.0, 0.0);
}

