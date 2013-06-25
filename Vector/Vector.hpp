/**
 * The vector types used in this simulation are specializations of Boost's
 * uBLAS vectors. In particular, we use ones composed of doubles and integers
 * for the coordinates of the simulations.
 *
 * Integral valued vectors are used for book keeping in the various data
 * structures.
 */
#ifndef __VECTOR_VECTOR_HPP__
#define __VECTOR_VECTOR_HPP__

#include <tvmet/Vector.h>

using namespace tvmet;
using namespace tvmet::element_wise;

/**
 * Vector type for representing particle positions and velocities as
 * double precision floating point values.
 */
typedef tvmet::Vector<double, 3UL> DoubleVector;

/**
 * A vector of integer values to hold a particle's position in the
 * grid of cells.
 */
typedef tvmet::Vector<int, 3UL> IntVector;

#endif /* __VECTOR_H__ */

