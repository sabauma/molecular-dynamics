
#ifndef __RANDOM_RANDOM_H__
#define __RANDOM_RANDOM_H__

#include "../Vector/Vector.hpp"

/**
 * Generates uniformly distributed random numbers in the range [0, 1)
 */
double uniform_rand();

/**
 * Generates a vector whose components have a standard normal distribution.
 */
void rand_vector(DoubleVector& v);

/**
 * Generates a velocity vector according to the Maxwell-Boltzman Distribution.
 */
void maxwell_boltzman_velocity(const double temp, const double mass, DoubleVector& v);

/**
 * Generates samples from the 2D Maxwell Distribution.
 */
double rand_maxwell_2d();

/**
 * Generates a random number accoring to the cosine distribution
 * f(x) = 0.5 * cos(x) for -pi/2 <= x <= pi/2.
 */
double cosine_distribution();

#endif
