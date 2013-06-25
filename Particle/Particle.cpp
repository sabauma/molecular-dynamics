/**
 * Implementation of the particle data structure.
 *
 * @author Spenser Bauman
 */

#include "../Constants/Constants.h"
#include "../Vector/Vector.hpp"
#include "Particle.h"

#include <cassert>
#include <iostream>

Particle::Particle()
    : Position(0.0), Velocity(0.0), Updated(0.0), Charge(0), Type(0)
{ }

Particle::Particle(const DoubleVector& position,
                   const DoubleVector& velocity,
                   const double updated,
                   const int type)
    : Position(position), Velocity(velocity), Updated(updated), Charge(0), Type(type)
{ }

Particle::Particle(const Particle& other) :
    Position(other.Position), Velocity(other.Velocity),
    Updated(other.Updated), Charge(other.Charge), Type(other.Type)
{ }

void
Particle::UpdateParticle(const double currentTime)
{
    this->Position += this->Velocity * (currentTime - this->Updated);
    this->Updated = currentTime;
}

std::ostream& operator<<(std::ostream& stream, const Particle& particle)
{
    return stream << "Position: " << particle.Position << std::endl
                  << "Velocity: " << particle.Velocity << std::endl
                  << "Type: "     << particle.Type     << std::endl
                  << "Charge: "   << particle.Charge;
}

/*
 *void Particle::Read(const boost::property_tree::ptree& tree)
 *{
 *    Updated = tree.get<double>("Updated");
 *    Charge  = tree.get<int>("Charge");
 *    Type    = tree.get<int>("Type");
 *
 *    ParseArray<DoubleVector>(tree.get_child("Position"),
 *                             Position.begin(), Position.end());
 *
 *    ParseArray<DoubleVector>(tree.get_child("Velocity"),
 *                             Velocity.begin(), Velocity.end());
 *}
 */

