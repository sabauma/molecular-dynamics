/**
 * Definition of a basic particle type. This module defines the particle type
 * used in the simulation. We make use of uBLAS vectors to provide
 * multidimensional vectors.
 *
 * Particles store, along with position and velocity, the time at which these
 * values were last updated (in system units). This lazy update strategy saves
 * the need of having to update the entire system after each event.
 *
 * @author Spenser Bauman
 */

#ifndef __PARTICLE_PARTICLE_H__
#define __PARTICLE_PARTICLE_H__

#include <iostream>

#include "rapidjson/prettywriter.h"
#include "rapidjson/document.h"
#include "../Vector/Vector.hpp"

struct Particle
{
    DoubleVector Position; /**< Position of the particle in simulation units.   */
    DoubleVector Velocity; /**< Velocity of the particle in simulation units.   */
    double       Updated;  /**< Time at which the simulation was last updated.  */
    int          Charge;   /**< Ionization level of the particle.               */
    int          Type;     /**< Index representing the element of the particle. */

    /**
     * Create a "null" particle with default values of zero. The position and
     * velocity vectors are still initialized to length 3.
     */
    Particle();

    /**
     * Initialize the particle with position and velocity to the given vectors
     * and the type to the given type.
     *
     * @param position The position of the particle in simulation units.
     * @param velocity The velocity of the particle in simulation units.
     * @param updated The time at which the particle was updated in simulation
     *                units.
     * @param type The index of the element of the particle.
     */
    Particle(const DoubleVector& position,
             const DoubleVector& velocity,
             const double updated,
             const int type);

    /**
     * Copy constructor. Perform a verbatim copy of the particle's
     * information.
     */
    Particle(const Particle& other);

    /**
     * Update the particle to the current system time. Computes the time
     * differential and then advances the particle that far based on its
     * velocity.
     */
    void UpdateParticle(const double currentTime);

    friend std::ostream& operator<<(std::ostream& stream, const Particle& particle);

    template <typename Writer>
    void Serialize(Writer& pt) const;
};

template <typename Writer>
void Particle::Serialize(Writer& writer) const
{
    writer.StartObject();

    writer.String("Updated") , writer.Double(Updated);
    writer.String("Charge")  , writer.Integral(Charge);
    writer.String("Type")    , writer.Integral(Type);

    writer.String("Position");
    writer.StartArray();
    {
        writer.Double(this->Position[0]);
        writer.Double(this->Position[1]);
        writer.Double(this->Position[2]);
    }
    writer.EndArray();

    writer.String("Velocity");
    writer.StartArray();
    {
        writer.Double(this->Velocity[0]);
        writer.Double(this->Velocity[1]);
        writer.Double(this->Velocity[2]);
    }
    writer.EndArray();

    writer.EndObject();
}

#endif /* __PARTICLE_H__ */
