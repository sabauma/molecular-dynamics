
#ifndef __FUNCTIONS_H__
#define __FUNCTIONS_H__

#include <vector>
#include <functional>
#include <numeric>

#include "../Constants/Constants.h"
#include "../Particle/Particle.h"

/**
 * Produces the value with the largest magnitude of the two given parameters.
 */
inline double
max_mag(const double t1, const double t2)
{
    return std::abs(t1) >= std::abs(t2) ? t1 : t2;
}

/**
 * Extract the sign of the (hopefully) numerical value given.
 */
template <typename T>
inline T sign(const T n)
{
    return ((T) (n > 0));
}

template <typename T>
inline double Sqr(const T x)
{
    return x * x;
}

template <typename T>
inline double Pow3(const T x)
{
    return x * x * x;
}

template <typename T>
T& identity(T& x)
{
    return x;
}

inline double
mean_velocity(const std::vector<Particle>& arr)
{
    assert(arr.size() > 0);
    double result = std::accumulate(arr.begin(), arr.end(), 0.0,
            [](double d, Particle p) { return d + norm2(p.Velocity); });

    return result / (double) arr.size();
}

inline double
median_velocity(const std::vector<Particle>& arr)
{
    assert(arr.size() > 0);
    const size_t size = arr.size();
    const size_t mid = size / 2;
    std::vector<double> vel(size);

    std::transform(arr.begin(), arr.end(), vel.begin(),
                   [](const Particle& p)
    {
        return norm2(p.Velocity);
    });

    if (size % 2 == 0)
    {
        double v1;
        nth_element(vel.begin(), vel.end()+mid, vel.end());
        v1 = vel[mid];
        nth_element(vel.begin(), vel.end()+mid-1, vel.end());

        return v1 + vel[mid-1] / 2.0;
    }
    else
    {
        nth_element(vel.begin(), vel.end() + mid, vel.end());
        return vel[mid];
    }
}

/**
 * Computes the time to intersection given a pair of positions, velocities,
 * and the minimum pass distance.
 *
 * @param p1 Position of the first particle.
 * @param v1 Velocity of the first particle.
 * @param p2 Position of the second particle.
 * @param v2 Velocity of the second particle.
 * @param diameter The minimum approach diameter of the two particles.
 *
 * @return The time until intersection.
 */
inline double
intersection_time(const DoubleVector& p1, const DoubleVector& v1,
                  const DoubleVector& p2, const DoubleVector& v2,
                  double diameter)
{
    const DoubleVector dr(p1 - p2);
    const DoubleVector dv(v1 - v2);
    const double b = dot(dr, dv);

    // We are interested in the minimum squared diameter.
    diameter *= diameter;

    if (b < 0.0)
    {
        const double rr = dot(dr, dr);
        const double vv = dot(dv, dv);

        if (rr < diameter)
        {
            diameter = rr * 0.99;
        }

        const double d = Sqr(b) - vv * (rr - diameter);

        if (d >= 0.0)
        {
            return -(sqrt(d) + b) / vv;
        }
    }

    return Constants::NEVER;
}

/**
 * Computes the new velocity vectors for a pair of colliding particles with
 * the given positions, velocities and masses.
 */
inline void
process_collision(const DoubleVector& p1, DoubleVector& v1, const double m1,
                  const DoubleVector& p2, DoubleVector& v2, const double m2)
{
    const DoubleVector dr(p1 - p2);
    const DoubleVector dv(v1 - v2);
    const double Minv = 1.0 / (m1 + m2);
    const double fac  = dot(dr, dv) / dot(dr, dr);

    v1 -= (2.0 * m2 * fac * Minv) * dr;
    v2 += (2.0 * m1 * fac * Minv) * dr;
}

#endif
