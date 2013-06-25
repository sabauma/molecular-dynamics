/**
 * PressureGauge.cpp - Pressure guage implementation.
 *
 * @author Spenser Bauman
 */

#include "../Constants/Constants.h"
#include "../Math/Functions.h"
#include "../Units/Units.h"
#include "rapidjson/prettywriter.h"
#include "PressureGauge.h"

#include <algorithm>
#include <boost/array.hpp>
#include <boost/circular_buffer.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <deque>
#include <functional>
#include <numeric>
#include <stdio.h>

/* The size of the buffer used for the pressure gauge. Too many values will
 * result in a time lag in the pressure, whereas too few will make the
 * pressure gauge unstable. */
static const size_t BUFFER_SIZE = 1000;

/**
 * Comparator function for sorting pairs of doubles.
 */
static bool
comp(const std::pair<double, double>& a, const std::pair<double, double>& b)
{
    return a.first < b.first;
}

static std::pair<double, double>
means(const PressureGauge::BufferType& dq)
{
    const double n = dq.size();
    double v1 = 0.0;
    double v2 = 0.0;

    for (auto it = dq.begin(); it != dq.end(); ++it)
    {
        v1 += it->first;
        v2 += it->second;
    }

    return std::make_pair(v1 / n, v2 / n);
}

PressureGauge::PressureGauge()
    : Buffer(BUFFER_SIZE),
      SampleRadius(0.0), SampleTime(0.0),
      LastPressure(0.0), LastDpdt(0.0), LastTime(0.0),
      CurrentPressure(0.0), CurrentDpdt(0.0), CurrentTime(0.0),
      Ready(false)
{ }

// Add new sample pair to sample arrays
void
PressureGauge::AddSample(const double time, const double pressure)
{
#ifdef USE_PGAUGE

    assert(!isnan(pressure) && !isnan(time) && !isinf(pressure) && !isinf(time));
    Buffer.push_back(std::make_pair(time, pressure));

#endif /* ifdef USE_PGAUGE */
}

/*
 *void
 *PressureGauge::Read(const ptree& tree)
 *{
 *    LastPressure    = tree.get<double>("LastPressure");
 *    LastDpdt        = tree.get<double>("LastDpdt");
 *    LastTime        = tree.get<double>("LastTime");
 *    CurrentPressure = tree.get<double>("CurrentPressure");
 *    CurrentDpdt     = tree.get<double>("CurrentDpdt");
 *    CurrentTime     = tree.get<double>("CurrentTime");
 *
 *    ptree first  = tree.get_child("Buffer.first");
 *    ptree second = tree.get_child("Buffer.second");
 *
 *    const std::vector<double> first_data  = ParseArray<double>(first);
 *    const std::vector<double> second_data = ParseArray<double>(second);
 *
 *    for (size_t i = 0; i < first_data.size(); ++i)
 *    {
 *        this->AddSample(first_data[i], second_data[i]);
 *    }
 *}
 */

bool
PressureGauge::IsReady() const
{
#ifdef USE_PGAUGE
    return Ready;
#else
    return false;
#endif
}

bool
PressureGauge::CanSample() const
{
    return Buffer.size() > 10;
}

void
PressureGauge::SetPressureStats(const double t, const double p, const double dpdt)
{
    this->CurrentTime     = t;
    this->CurrentPressure = p;
    this->CurrentDpdt     = dpdt;
}

double
PressureGauge::PressureTime() const
{
    return means(Buffer).first;
}

// Pressure & time samples are populated and processed in a round-robin fasion
double
PressureGauge::SampleGauge(const double radius)
{
    if (!this->CanSample())
    {
        return CurrentTime;
    }

    Ready = true;

    std::sort(Buffer.begin(), Buffer.end(), comp);

    const double t_start = Buffer.front().first;
    const double t_end   = Buffer.back().first;

    double y_sum = 0.0;
    double t_sum = 0.0;

    // Current data becomes the previous data.
    LastPressure = CurrentPressure;
    LastDpdt     = CurrentDpdt;
    LastTime     = CurrentTime;

    for (size_t i = 0; i < Buffer.size(); ++i)
    {
        t_sum += Buffer[i].first;
        y_sum += Buffer[i].second;
    }

    const double N  = Buffer.size();
    const double dt = t_end - t_start;

    assert(dt > 0.0);
    assert(y_sum >= 0.0);

    // Computed pressure is the change in momentum per area per time
    CurrentPressure = y_sum / dt;

    // Computed times is the mean time of the used pressure events.
    CurrentTime     = t_sum / N;

    // dpdt is computed as the slope of the line between the current and last
    // pressure.
    if (CurrentTime != LastTime)
    {
        CurrentDpdt = (CurrentPressure - LastPressure) / (CurrentTime - LastTime);
    }
    else
    {
        CurrentDpdt = LastDpdt;
    }

    SampleRadius = radius;

    Buffer.clear();

    return CurrentTime;
}

double
PressureGauge::ComputePressure(
        const double currentTime,
        const double radius,
        const double drdt) const
{
#ifdef ADIABATIC_INTERPOLATION
    return CurrentPressure * pow(SampleRadius / radius, 3.0 * Constants::GAMMA);
#else
    return CurrentPressure + (currentTime - CurrentTime) * CurrentDpdt;
#endif
}

double
PressureGauge::ComputeDpdt(
        const double currentTime,
        const double radius,
        const double drdt) const
{
#ifdef ADIABATIC_INTERPOLATION
    return -3.0 * Constants::GAMMA
        * this->ComputePressure(currentTime, radius, drdt)
        * drdt / radius;
#else
    return CurrentDpdt;
#endif
}

