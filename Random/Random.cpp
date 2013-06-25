
#include <boost/random.hpp>
#include <math.h>
#include <time.h>

#include "../Vector/Vector.hpp"
#include "Random.h"

using namespace boost;

#define SEED time(0)
// #define SEED 1

static variate_generator< mt19937, normal_distribution<> >
    norm_rnd(mt19937(SEED), normal_distribution<>(0.0, 1.0));

static variate_generator< mt19937, uniform_01<> >
    uniform_rnd(mt19937(SEED), uniform_01<>());

double uniform_rand()
{
    return uniform_rnd();
}

void rand_vector(DoubleVector& v)
{
    v[0] = norm_rnd();
    v[1] = norm_rnd();
    v[2] = norm_rnd();
}

void maxwell_boltzman_velocity(const double temp, const double mass, DoubleVector& v)
{
    rand_vector(v);
    v *= sqrt(temp / mass);
}

double rand_maxwell_2d()
{
    /*
     * Generates a positive number with p(v) = 2v*exp(-v^2), v > 0 This serves
     * dual purpose:
     *
     * a) this is the distribution law for 2D Maxwell
     * b) this is also the distribution for the normal
     *    component of the particle coming
     *    off the wall Multiply by sqrt (2*T)
     */
    return sqrt(-log(1-uniform_rand()));
}

double cosine_distribution()
{
    double r = uniform_rand();
    return asin(2.0 * r - 1.0);
}
