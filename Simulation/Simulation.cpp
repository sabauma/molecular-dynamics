/**
 * Implementation of the simulation class. Provides all the operations needed
 * to manipulate and analyze the simulation.
 *
 * @author Spenser Bauman
 */

#include "../CellList/CellListFast.h"
#include "../Constants/AtomicProperties.h"
#include "../Constants/Constants.h"
#include "../Debug/Trace.h"
#include "../Math/Functions.h"
#include "../Random/Random.h"
#include "../Vector/Vector.hpp"
#include "Simulation.h"
#include "rapidjson/writer.h"
#include "rapidjson/filewritestream.h"
#include "rapidjson/prettywriter.h"

#include <algorithm>
#include <boost/date_time/local_time/local_time.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <cassert>
#include <fstream>
#include <iostream>
#include <limits>
#include <math.h>
#include <string>
#include <time.h>
#include <vector>

// Controls the error rate in the spatial partitioning algorithm.
// This rate effectively chooses the quantile to use in the distiribution
// of collision cross sections. The closer it is to 1, the smaller the
// error will be, though setting it to 1 does not guarantee zero error,
// since the actual partition size is derived from a sample of all collisions.
static const double ERROR_RATE = 0.990;

// Shrink the region when the VOLUME is cut in half, as opposed to the radius,
// which could result in a high cell occupancy.
static const double SHRINK_RADIUS = pow(0.5, 1.0 / 3.0);

static const int CONST_CELL_FACTOR = 2;
static const size_t BUFFER_SIZE = 65536;

/**
 * Performs a rotation of the vector v around the vector n based of the
 * given angle (in radians). All vectors are assumed to have a dimension of
 * 3, otherwise this will not be a valid transformation.
 *
 * @param v The vector to rotate.
 * @param n The vector to about which the rotation takes place.
 * @param theta The angle of rotation in degrees.
 */
static inline void
RotateVector(DoubleVector& v, const DoubleVector& n, const double theta)
{
    DoubleVector vn(n);
    DoubleVector vt(cross(v, n));

    // The projection of v onto n and remove the projection from n.
    // n is now orthogonal to v.
    vn *= dot(n, v);
    v  -= vn;

    // Rotate the orthogonal components.
    vt *= sin(theta);
    v  *= cos(theta);

    // Reconstitue the vector in rotated form.
    v  += vn;
    v  += vt;
}

/**
 * For now, a stub function to act as the diameter of each particle. This will
 * actually return the square of the diameter of the pair of particles given
 * as t1 and t2.
 */
double
Simulation::Diameter2(const double vv, const Particle& p1, const Particle& p2)
{
#ifdef VSS_PARAMETERS
    const int type1 = p1.Type;
    const int type2 = p2.Type;
    const double ene = Info.ReducedMass[type1][type2] * vv;

    const double d1 = sqrt(Info.BirdConstant[type1] * sqrt(Info.AtomicMass[type1])
                           / pow(ene, Constants::BIRD_OMEGA_INI[type1]));
    const double d2 = sqrt(Info.BirdConstant[type2] * sqrt(Info.AtomicMass[type2])
                           / pow(ene, Constants::BIRD_OMEGA_INI[type2]));

    return Sqr(0.5 * (d1 + d2));
    //return 0.5 * ((Info.BirdConstant[type1]
                  /// pow(vv, Constants::BIRD_OMEGA_INI[type1] - 0.5))
                  //+ (Info.BirdConstant[type2]
                  /// pow(vv, Constants::BIRD_OMEGA_INI[type2] - 0.5)));
#else
    return Info.AtomicDiameterSquared[p1.Type][p2.Type];
#endif // VSS_PARAMETERS
}

int Simulation::CellCrossingTimes(const int na, DoubleVector& times)
{
    // The current particle being consulted.
    Particle& pa = this->Particles[na];

    const IntVector cell(Partitions.ParticlePosition(na));

    assert(pa.Updated == CurrentTime);

    // Vector operators are elementwise, so use the operations to compute all
    // three times at once.
    times = (((cell + (pa.Velocity > 0)) * Extent / Cells)
            - pa.Position - Offset) / pa.Velocity;

    // assert(all_elements(times >= 0.0));

    // Return the axis with the shortest time to intersection.
    return times[1] < times[2]
        ? (times[0] < times[1] ? 0 : 1)
        : (times[0] < times[2] ? 0 : 2);
}

double Simulation::ConeCollisionTime(const int na)
{
    // The particle for the simulation, its velocity and position.
    Particle& pa = this->Particles[na];

    const DoubleVector& vel = pa.Velocity;
    const DoubleVector& pos = pa.Position;

    // Terms in the quadratic equation for cone collision time.
    const double ca = Sqr(vel[0])
                    + Sqr(vel[1])
                    - Sqr(vel[2]) * Info.ConeTan2;

    const double cb = vel[0] * pos[0]
                    + vel[1] * pos[1]
                    - vel[2] * pos[2] * Info.ConeTan2;

    const double cc = Sqr(pos[0])
                    + Sqr(pos[1])
                    - Sqr(pos[2]) * Info.ConeTan2;

    // Solve the quadratic equation in a way that reduces the numerical error
    // as much as is possible.
    if (ca > 0)
    {
        return cb > 0
            ? -cc / (cb + sqrt(Sqr(cb) - ca * cc))
            : (-cb + sqrt(Sqr(cb) - ca * cc)) / ca;
    }
    else
    {
        return cb > 0
            ? -cc / (cb + sqrt(Sqr(cb) - ca * cc))
            : Constants::NEVER;
    }
}

double Simulation::WallCollisionTime(const int na)
{
    assert(na >= 0);
    assert((unsigned) na < Particles.size());

    // Make use of the RP equations to interpolate.
    return RPKEquation.FindCrossingTime(Particles[na]);
}

/**
 * Code is similar to what is found in <code>intersection_time</code>, but
 * is optimized to work on particle types. This allows for the delay of the
 * computation of the diameter parameter.
 */
double Simulation::ComputeIntersectionTime(const Particle& pa, const Particle& pb)
{
    // Compute the differences in velocities and posisions.
    const DoubleVector dr(pa.Position - pb.Position);
    const DoubleVector dv(pa.Velocity - pb.Velocity);

    // Use the dot product to ensure the dr and dv are not pointing in
    // opposite directions (i.e. they are getting closer to each other).
    const double b = dot(dr, dv);

    // Check that they are getting closer.
    if (b < 0.0)
    {
        const double rr  = dot(dr, dr);
        const double vv  = dot(dv, dv);
        double diameter2 = Diameter2(vv, pa, pb);

        // Check that the particles are not inside of eachother due to
        // numerical error.
        if (rr < diameter2)
        {
            diameter2 = rr * 0.99;
        }

        const double d = Sqr(b) - vv * (rr - diameter2);

        // If the equation has a solution, return the time computed.
        if (d >= 0.0)
        {
            return -(sqrt(d) + b) / vv;
        }
    }

    // Otherwise, there will never be an intersection.
    return Constants::NEVER;
}

void Simulation::ProcessCellCrossing()
{
    const int na   = Calendar->GetCurrentObjectA();
    const int axis = Calendar->GetCurrentObjectB();

    Particle& pa = Particles[na];
    pa.UpdateParticle(CurrentTime);

    // Positive or negative direction along that axis.
    if (pa.Velocity[axis] >= 0.0)
    {
        Partitions.MoveParticle(na, Constants::STD_BASIS[axis]);
    }
    else
    {
        Partitions.MoveParticle(na, Constants::NEG_STD_BASIS[axis]);
    }

    // Repredict events based on the last event being a cell crossing.
    // Predict collisions with all particles, but cell crossing events are
    // optimized in the <code>PredictEvents</code> method to check only
    // certain voxels along the particle's path.
    PredictEvents<AllParticles>(na, axis, CellCrossingEvent);
}

void Simulation::ProcessBubbleWall()
{
    const int na = Calendar->GetCurrentObjectA();

    Particle& pa = Particles[na];
    pa.UpdateParticle(CurrentTime);

    // Normal vector at the point of intersection on the bubble wall.
    const DoubleVector normal(-normalize(pa.Position));
    const DoubleVector v_ini(pa.Velocity);

    // TODO: possibly use the actual computed RPK values of velocity, rather
    // than interpolating the position.
    const double vWall = RPKEquation.V(CurrentTime);

#ifdef SPECULAR_WALL
    pa.Velocity -= 2.0 * (dot(pa.Velocity, normal) + vWall) * normal;
#else

    DoubleVector perp;
    if (std::abs(normal[0]) <= std::abs(normal[1]))
    {
        if (std::abs(normal[0]) <= std::abs(normal[2]))
        {
            perp = DoubleVector(0, -normal[2], normal[1]);
        }
        else
        {
            perp = DoubleVector(-normal[1], normal[0], 0);
        }
    }
    else
    {
        if (std::abs(normal[1]) <= std::abs(normal[2]))
        {
            perp = DoubleVector(-normal[2], 0, normal[0]);
        }
        else
        {
            perp = DoubleVector(-normal[1], normal[0], 0);
        }
    }

    perp = normalize(perp);
    RotateVector(perp, normal, 2.0 * Constants::PI * uniform_rand());

    pa.Velocity = (Info.VWallThermal[pa.Type] * rand_maxwell_2d() - vWall) * normal +
                  (Info.VWallThermal[pa.Type] * rand_maxwell_2d()) * perp;

#endif

    const double dV   = std::abs(dot(pa.Velocity - v_ini, normal));
    const double r    = RPKEquation.R(CurrentTime);
    const double area = Info.ConeSolidAngle * r * r;

    // We add one pressure sample as the change in momentum per unit
    // of area.
    const double momentum = Info.AtomicMass[pa.Type] * dV / area;

    // Add the force of the particle to the pressure gauge.
    RPKEquation.AddSample(CurrentTime, momentum);
    // Predict new events for the particle.
    PredictEvents<AllParticles>(na, -1, BubbleWallEvent);
}

void Simulation::ProcessCollision()
{
    const int na = Calendar->GetCurrentObjectA();
    const int nb = Calendar->GetCurrentObjectB();

    Particle& pa = Particles[na];
    Particle& pb = Particles[nb];

    // Update the particles to their positions at CurrentTime
    pa.UpdateParticle(CurrentTime);
    pb.UpdateParticle(CurrentTime);

    const DoubleVector dr(pa.Position - pb.Position);
    const DoubleVector dv(pa.Velocity - pb.Velocity);

    const int typeA   = pa.Type;
    const int typeB   = pb.Type;
    const double mA   = Info.AtomicMass[typeA];
    const double mB   = Info.AtomicMass[typeB];
    const double Minv = 1.0 / (mA + mB);
    const double fac  = dot(dr, dv) / dot(dr, dr);

    // Compute the change in velocities.
    pa.Velocity -= (2.0 * mB * fac * Minv) * dr;
    pb.Velocity += (2.0 * mA * fac * Minv) * dr;

    /*****************************************************************************\
    *  Steve's Comment: If collision energy>ionization energy we have work to do: *
    *  1. Compute the velocity of the center of mass, vCoM                        *
    *  and the relative velocities of particle A and B                            *
    *  2. Compute the kinetic energy of A (energyA) and B (energyB)               *
    *  3. If the total energy (energyA+energyB) is enough to                      *
    *  strip an electron off one of the atoms, strip it off the                   *
    *  one that requires less energy to do so.  Subtract enough                   *
    *  energy from the system to take into account of the ionization.             *
    *  This energy subtraction step is done in a totally ad-hoc way.              *
    *  We use something that assures no negative energies.                        *
    *  4. Now we have 3 particles in the system, A, B and the electron.           *
    *  Assign 1/3 of the energy to the electron and drop it from                  *
    *  the calculation.  Multiply the energy for particles A and                  *
    *  B by 2/3.                                                                  *
    *  CORRECTION: Now we do not assign any energy to the electron (!)            *
    \*****************************************************************************/

    const double mReduced = Info.ReducedMass[typeA][typeB];
    // K = 0.5 * m * v^2
    double ECoM = mReduced * dot(dv, dv); // total energy in CoM

    const DoubleVector CoM(0.5 * (pa.Position + pb.Position));

    CollisionRecord rec;
    rec.Position = CoM;
    rec.Time     = CurrentTime;
    rec.DeltaV   = norm2(dv);
    rec.Energy   = ECoM;
    rec.Distance = norm2(dr);
    rec.Type1    = pa.Type;
    rec.Type2    = pb.Type;

    //Stats.Collisions.push_back(rec);
    //Stats.MaxCollisionEnergy = std::max(Stats.MaxCollisionEnergy, ECoM);

    Stats.RegisterCollision(rec);

    const double CONST_ENERGY_TO_KEEP = 1.0;

    // Check for dissociations. When a dissociation occurs, the collision
    // processing returns, rather than checking for ionizations and fusions.
    // Dissociation requires special logic for event detection.
    if (Constants::DISSOCIATION_ENERGY[typeA]
            < Constants::DISSOCIATION_ENERGY[typeB])
    {
        if (Constants::DISSOCIATION_ENERGY[typeA] < ECoM)
        {
            const DoubleVector vCoM((mA * pa.Velocity + mB * pb.Velocity) * Minv);
            const double dis = Constants::DISSOCIATION_ENERGY[typeA];
            const double k   = sqrt(CONST_ENERGY_TO_KEEP * (ECoM - dis) / ECoM);

            pa.Velocity = vCoM - mB*Minv*k*dv;
            pb.Velocity = vCoM + mA*Minv*k*dv;

            DissociateParticle(na);
            PredictEvents<AllParticles>(nb, na, CollisionEvent);

            return;
        }
    }
    else
    {
        if (Constants::DISSOCIATION_ENERGY[typeB] < ECoM)
        {
            const DoubleVector vCoM((mA * pa.Velocity + mB * pb.Velocity) * Minv);
            const double dis = Constants::DISSOCIATION_ENERGY[typeB];
            const double k   = sqrt(CONST_ENERGY_TO_KEEP * (ECoM - dis) / ECoM);

            pa.Velocity = vCoM - mB*Minv*k*dv;
            pb.Velocity = vCoM + mA*Minv*k*dv;

            DissociateParticle(nb);
            PredictEvents<AllParticles>(na, nb, CollisionEvent);

            return;
        }
    }

    // Check for ionizations
    if (Constants::IONIZATION_ENERGY_INI[typeA][pa.Charge]
            < Constants::IONIZATION_ENERGY_INI[typeB][pb.Charge])
    {
        if (Constants::IONIZATION_ENERGY_INI[typeA][pa.Charge] < ECoM)
        {
            const DoubleVector vCoM((mA * pa.Velocity + mB * pb.Velocity) * Minv);
            const double ion = Constants::IONIZATION_ENERGY_INI[typeA][pa.Charge];
            const double k   = sqrt(CONST_ENERGY_TO_KEEP * (ECoM - ion) / ECoM);

            pa.Velocity = vCoM - mB*Minv*k*dv;
            pb.Velocity = vCoM + mA*Minv*k*dv;
            ++pa.Charge;
            ECoM -= ion;
        }
    }
    else
    {
        if (Constants::IONIZATION_ENERGY_INI[typeB][pb.Charge] < ECoM)
        {
            const DoubleVector vCoM((mA * pa.Velocity + mB * pb.Velocity) * Minv);
            const double ion = Constants::IONIZATION_ENERGY_INI[typeB][pb.Charge];
            const double k   = sqrt(CONST_ENERGY_TO_KEEP * (ECoM - ion) / ECoM);

            pa.Velocity = vCoM - mB*Minv*k*dv;
            pb.Velocity = vCoM + mA*Minv*k*dv;
            ++pb.Charge;
            ECoM -= ion;
        }
    }

    // Predict all events for the first particle.
    PredictEvents<AllParticles>(na, nb, CollisionEvent);
    // Predict all events for nb skipping na, as we already checked it.
    PredictEvents<Skipping>(nb, na, CollisionEvent);

    // Fusion occurs only between two particles of the lighest elements.
    // This assumes that type index 0 corresponds to Hydrogen.
    if (typeA == typeB && Constants::FUSABLE[typeA])
    {
        const double T = ECoM * 2.0 / Constants::KB * Units::Energy;

        // If the temperature is higher then Coulomb barier penetration
        if (T > 4.5e7)
        {
            ++Stats.FusionRate;
        }
    }
}

void Simulation::DissociateParticle(const int na)
{
    assert(ParticleCount < MaxParticles);

    const int nb = ParticleCount++;

    Particle& pa = Particles[na];
    Particle& pb = Particles[nb];

    const std::string& new_name = Constants::DISSOCIATES_TO_NAME[pa.Type];
    const int new_type          = Constants::TYPE_FROM_NAME.at(new_name);

    // Assign the new type to each of the particles.
    pa.Type = pb.Type = new_type;

    pa.Updated = pb.Updated = CurrentTime;

    // Each particle has the velocity of the original
    pb.Velocity = pa.Velocity;

    const DoubleVector radial(normalize(pa.Position));

    // The second particle is placed beside the first particle, but translated
    // radially inward. This prevents problems with crossing the bubble
    // or cone wall.
    pb.Position = pa.Position - radial * Info.AtomicDiameter[pa.Type];

    // If we translate through the origin, set us at the origin. This is
    // possible, since we move it radially inward by an atomic diameter.
    if (dot(pb.Position, Constants::Z) < 0.0)
    {
        pb.Position = 0.0;
    }

    // Register the new particle with the cell list
    const IntVector inCell((pb.Position + Offset) * Cells / Extent);
    Partitions.SetParticle(nb, inCell);

    // Recognize events for the new particle.
    PredictEvents<AllParticles>(na, nb, CollisionEvent);
    PredictEvents<Skipping>(nb, na, CollisionEvent);
}

void Simulation::ProcessCone()
{
    const int na = Calendar->GetCurrentObjectA();

    Particle& pa = Particles[na];
    pa.UpdateParticle(CurrentTime);

    DoubleVector normal(pa.Position);
    normal[2] = -sqrt((Sqr(normal[0]) + Sqr(normal[1])) * Info.ConeTan2);

    // Specular reflection
    pa.Velocity -= 2.0 * dot(normal, pa.Velocity) / dot(normal, normal) * normal;

    PredictEvents<AllParticles>(na, -1, ConeWallEvent);
}

void Simulation::ProcessUpdateSystem()
{
    printf("\nParticle Count: %d\n", ParticleCount);
    for (int i = 0; i < ParticleCount; ++i)
    {
        Particles[i].UpdateParticle(CurrentTime);
    }

    // Current bubble radius from RP-equation
    const double currentRadius   = RPKEquation.R(CurrentTime);
    const double currentVelocity = RPKEquation.V(CurrentTime);

    if (currentVelocity > 0 && LastVelocity < 0 && !Rebounded)
    {
        Rebounded = true;
    }

    LastVelocity = currentVelocity;

    Stats.ComputeStatistics(begin(), end(), currentRadius);

    const double collision_distance = Stats.Collisions.size() > 500
                                    ? Stats.CollisionDistanceQuantile(ERROR_RATE)
                                    : 1.0;

    printf("Collision Distance: %e\n", collision_distance);

    SnapshotCount++;

    WriteSnapshot();
    WriteCollisionEnergy();

    double nextStepTime = RPKEquation.NextStepTime();

    if (CurrentTime >= nextStepTime)
    {
        bool rebin = false;
        // Calculate the next step for RP
        nextStepTime = RPKEquation.NextStep();
        printf("### currentRadius: %0.12e ~ LastRadius %0.12e -> %0.12e\n",
               currentRadius, LastRebinRadius, currentRadius / LastRebinRadius);

        // Rebin every time the radius drops by 0.8. This means the total
        // volume decreased by about half. Otherwise, we can attempt to shrink
        // when the last rebin radius was set to zero, which indicates that it
        // was unable to shrink due to cell length concerns. These may not be
        // valid, as the energy in the simulation may have increased and made
        // the necessarry cell length shorter.
        if (currentRadius < SHRINK_RADIUS * LastRebinRadius || LastRebinRadius == 0.0)
        {
            printf("Re-binning(shrinking)\n");

            // Shrink regions by currentRadius/LastRebinRadius. The 1.001
            // factor makes the region slightly larger particles do not end up
            // outside the bounding box and to compensate for numerical error.
            const double shrink_factor = LastRebinRadius != 0.0
                                       ? 1.001 * currentRadius / LastRebinRadius
                                       : 0.0;

            Extent *= shrink_factor;
            Offset *= shrink_factor;

            // Do not shrink below the size of the largest particle in the
            // system.
            if (Extent[0] < Cells[0] * collision_distance * 1.01)
            {
                printf("*** Preventing over shrinkage\n");
                Extent    = Cells * collision_distance * 1.01;
                Offset    = 0.5 * Extent;
                Offset[2] = 0.0;

                // Terminate rebinning when Region shrinks under cell size
                LastRebinRadius = 0.0;
            }
            else
            {
                LastRebinRadius = currentRadius;
            }

            rebin = true;

            printf("End Re-binning\n");
        }

        const double nextRadius = RPKEquation.R(nextStepTime);

        // Expand Region (bubble is growing or we shrunk too far).
        if (currentRadius >= Extent[2] || nextRadius >= Extent[2])
        {
            printf("Re-binning(expanding)\n");

            const double newRadius = std::max(currentRadius, nextRadius) * 1.01;
            // Maintain the proper ratio such that the cells are square and
            // snugly fit the cone.
            // Extent = newRadius;

            Extent    = newRadius * Cells / Cells[2];
            Offset    = 0.5 * Extent;
            Offset[2] = 0.0;
            // Extent *= 2.0;
            // Offset *= 2.0;

            LastRebinRadius = 0.0;

            rebin = true;

            printf("End Re-binning\n");
        }

        if (rebin)
        {
            BinSystem();
        }

        // Perform this operation all the time rather than mucking around
        // with UpdateWallCollisionTimes. This is not that expensive of
        // an operation. This will rebin the simulation and clear/repopulate
        // the event calendar with new events.
        Calendar->ClearCalendar();
        PopulateCalendar();

        // Update wall collision times as the RPK solver was updated. This
        // will provide more accurate times estimates.
        // UpdateWallCollisionTimes();
    }

    Calendar->ScheduleEvent(nextStepTime, UpdateSystemEvent);
}

void Simulation::BinSystem()
{
    Partitions.Clear();
    for (int i = 0; i < ParticleCount; ++i)
    {
        Particle& pi = Particles[i];
        pi.UpdateParticle(CurrentTime);
        const IntVector inCell((pi.Position + Offset) * Cells / Extent);

        Partitions.SetParticle(i, inCell);
    }
}

void Simulation::StartRun()
{
    BinSystem();
    Calendar->ClearCalendar();
    PopulateCalendar();
}

int Simulation::ComputeParticlesInSystem() const
{
    DoubleVector c;
    DoubleVector gap(Extent / InitUCell);
    DoubleVector rp;

    int n = 0;

    for (int nZ = 0; nZ < InitUCell[2]; nZ++)
    {
        for (int nY = 0; nY < InitUCell[1]; nY++)
        {
            for (int nX = 0; nX < InitUCell[0]; nX++)
            {
                c = gap * DoubleVector(nX + 0.25, nY + 0.25, nZ + 0.25) - Offset;

                for (int j = 0; j < 4; j++) // Not sure what this loop does.
                {
                    for ( int k = 0; k < Constants::DIMENSIONS; k++)
                    {
                        if (j == k || j == 3)
                            { rp[k] = c[k]; }
                        else
                            { rp[k] = c[k] + gap[k] * 0.5; }
                    }

                    double in = Sqr(rp[0]) + Sqr(rp[1]) - Sqr(rp[2]) * Info.ConeTan2;

                    if (in < 0 && dot(rp, rp) < Sqr(Info.Rini))
                    {
                        n++;
                    }
                }
            }
        }
    }

    return n;
}

void Simulation::InitializeParticles()
{
    DoubleVector gap(Extent / InitUCell);
    DoubleVector c;
    DoubleVector rp;

    const int MaxAtoms = ParticleCount;

    // Count of atoms of particular species
    double nAtomOfType[Constants::ELEMENT_TYPES];

    for (int i = 0; i < Constants::ELEMENT_TYPES; i++)
    {
        nAtomOfType[i] = MaxAtoms * Info.AtomicParts[i];
    }

    int nAtomLeft = MaxAtoms;
    int ndx = 0;

    for ( int nZ = 0; nZ < InitUCell[2]; nZ++)
    {
        for ( int nY = 0; nY < InitUCell[1]; nY++)
        {
            for ( int nX = 0; nX < InitUCell[0]; nX++)
            {
                c = gap * DoubleVector(nX + 0.25, nY + 0.25, nZ + 0.25) - Offset;

                for (int j = 0; j < 4; j++)
                {
                    for (int k = 0; k < Constants::DIMENSIONS; k++)
                    {
                        if (j == k || j == 3)
                        {
                            rp[k] = c[k];
                        }
                        else
                        {
                            rp[k] = c[k] + gap[k] * 0.5;
                        }
                    }

                    double in = Sqr(rp[0]) + Sqr(rp[1]) - Sqr(rp[2])*Info.ConeTan2;

                    if  (in < 0 && dot(rp, rp) < Sqr(Info.Rini))
                    {
                        // Choose atom type
                        double dice = (nAtomLeft--) * uniform_rand();
                        double sum = 0;

                        int type = 0;
                        while ((sum += nAtomOfType[type]) < dice)
                        {
                            type++;
                        }

                        nAtomOfType[type] -= 1.0;
                        Particles[ndx].Type = type;
                        Particles[ndx].Position = rp;
                        ++ndx;
                    }
                }
            }
        }
    }

    double thermalVelocity[Constants::ELEMENT_TYPES];

    // Each element is characterized by thermal velocity
    for (int i = 0; i < Constants::ELEMENT_TYPES; i++)
    {
        // Temperature in MD units
        thermalVelocity[i] = sqrt(Info.Temperature / Info.AtomicMass[i]);
    }

    for (int i = 0; i < ParticleCount; ++i)
    {
        rand_vector(Particles[i].Velocity);
        Particles[i].Velocity *= thermalVelocity[Particles[i].Type];
    }
}

void Simulation::PopulateCalendar()
{
    // Predict collisions between all particles less than the current one.
    for (int i = 0; i < ParticleCount; ++i)
    {
        PredictEvents<LessThan>((int) i, (int) i, CollisionEvent);
    }
}

/*********************************************************************************/
/*                                    Interface                                  */
/*********************************************************************************/
Simulation::Simulation(const InfoStruct& info) :
    Info(info),
    Partitions(0, 0, 0, 0),
    Particles(0),
    ConeTime(0, 0),
    WallTime(0, 0),
    CurrentTime(0),
    SnapshotCount(0),
    RPKEquation(CurrentTime, info),
    LastVelocity(0.0),
    LastSnapshotTime(0.0),
    Rebounded(false),
    Stats(Info, MAX_BINS, COLLISION_RECORDS),
    PerformanceClock(0)
{
    // Calculate partition information.
    InitUCell[0]
        = InitUCell[1]
        = int(ceil(2.0 * Info.Rini * sin(Info.ConeAngle) *
                   pow(0.25 * Info.Density, 1.0 / 3.0)));

    InitUCell[2] = int(ceil(Info.Rini * pow(0.25 * Info.Density, 1.0 / 3.0)));

    Extent = InitUCell / pow(Info.Density / 4.0, 1.0/3.0);

    Offset[0] = Offset[1] = 0.5 * Extent[0];
    Offset[2] = 0.0;

    Cells = CONST_CELL_FACTOR * InitUCell;

    // Get the number of particles in the system and init the structures that
    // depend on that information. Reserve twice the number of particle as
    // there are at initialization for dissociation.
    ParticleCount = ComputeParticlesInSystem();
    MaxParticles  = ParticleCount * 2;

    Particles.resize(MaxParticles);
    ConeTime.resize(MaxParticles);
    WallTime.resize(MaxParticles);

    std::fill(ConeTime.begin(), ConeTime.end(), 0.0);
    std::fill(WallTime.begin(), WallTime.end(), 0.0);

    Partitions.Resize(Particles.size(), Cells);

    Calendar = std::unique_ptr<EventCalendar>(new EventCalendar(MaxParticles));

    // Initialize the coordinates of each particle.
    InitializeParticles();

    // Must call next step before the system can be populated with
    // events, otherwise cell wall crossing will be wrong.
    RPKEquation.NextStep();

    LastRebinRadius = Info.Rini;
    LastVelocity    = RPKEquation.V(CurrentTime);

    // Binning
    StartRun();

    Calendar->ScheduleEvent(0, UpdateSystemEvent);
    // this->Calendar->ScheduleEvent(nextStepTime, UpdateSystemEvent);
}

bool Simulation::HasRebounded() const
{
    return Rebounded;
}

EventType Simulation::NextStep()
{
    assert(Calendar->HasMoreEvents());

    Calendar->NextEvent();

    // Time is monotonically increasing.
    // assert(this->CurrentTime <= this->Calendar->GetCurrentTime());
    CurrentTime = Calendar->GetCurrentTime();

    const EventType et = Calendar->GetCurrentEventType();

    switch (et)
    {
        case CellCrossingEvent:
            ProcessCellCrossing();
            ++Stats.TotalCellCrossings;
            break;
        case BubbleWallEvent:
            ProcessBubbleWall();
            ++Stats.TotalWallCollisions;
            break;
        case ConeWallEvent:
            ProcessCone();
            ++Stats.TotalConeBoundaryCollisions;
            break;
        case CollisionEvent:
            ProcessCollision();
            ++Stats.TotalCollisions;
            break;
        case UpdateSystemEvent:
            ProcessUpdateSystem();
            break;
        default:
            // This should not happen.
            printf("Unknown event type %d\nContinuing Anyway...", (int) et);
            break;
    }

    return et;
}

void Simulation::PurgeEvents()
{
    while (Calendar->HasMoreEvents())
    {
        Calendar->NextEvent();
    }
}

std::vector<Particle>::const_iterator
Simulation::begin()
{
    return Particles.begin();
}

std::vector<Particle>::const_iterator
Simulation::end()
{
    return Particles.begin() + ParticleCount;
}

double
Simulation::GetCurrentTime() const
{
    return CurrentTime;
}

std::ostream& operator<<(std::ostream& stream, const Simulation& sim)
{
    return stream << "Particles: "   << sim.ParticleCount << '\n'
           << "CurrentTime: " << sim.CurrentTime      << '\n'
           << "Extent: "      << sim.Extent           << '\n'
           << "Offset: "      << sim.Offset           << '\n'
           << "Cells: "       << sim.Cells            << '\n'
           << "InitUCell: "   << sim.InitUCell        << std::endl;
}

unsigned long long Simulation::GetCollisionCount() const
{
    return Stats.TotalCollisions;
}

unsigned long long Simulation::GetCellCrossingCount() const
{
    return Stats.TotalCellCrossings;
}

unsigned long long Simulation::GetWallCollisionCount() const
{
    return Stats.TotalWallCollisions;
}

unsigned long long Simulation::GetConeBoundaryCount() const
{
    return Stats.TotalConeBoundaryCollisions;
}

void Simulation::WriteCollisionEnergy()
{
    // Initialize snapshot file name
    char fileName[128] = {'\0'};
    sprintf(fileName, "CollisionEnergy%08d.csv", SnapshotCount);

    std::ofstream outfile(fileName);
    outfile.precision(15);

    Stats.WriteCollisions(outfile);
    Stats.Collisions.clear();
    outfile.close();
}

void Simulation::WriteSnapshot()
{
    using namespace rapidjson;

    // Buffer used for writing
    char writeBuffer[BUFFER_SIZE];
    // Initialize snapshot file name
    char fileName[128] = {'\0'};
    sprintf(fileName, "Snapshot%08d.dat", SnapshotCount);

    FILE* outfile = fopen(fileName, "w");
    assert(NULL != outfile);
    FileWriteStream os(outfile, writeBuffer, sizeof(writeBuffer));
    PrettyWriter<FileWriteStream> writer(os);
    writer.StartObject(); // Start the JSON structure

    const unsigned int MaxAtoms = ParticleCount;

    // Bubble radius offset ~1/2 an atomic diameter.
    const double currentRadius = RPKEquation.R(CurrentTime) + 0.5;

    const double delta_t =
        (double) (clock() - PerformanceClock) / CLOCKS_PER_SEC;

    const double me =
        Stats.MaxCollisionEnergy * Units::Energy * 6.24150974e15;

    // Print Summary
    printf("\nSnapshot #%d\n dT=%g s\n Time=%g s\n Radius=%g μm\n Collisions=%llu\n"
           " Cell Crossings=%llu\n WallCount=%llu\n ConeCount=%llu\n"
           " EnergyPerAtom=%g KeV\n CollisionEnergy=%0.12e KeV\n"
           " FusionRate=%.12e\n",
           SnapshotCount,
           double(clock() - PerformanceClock) / CLOCKS_PER_SEC,
           CurrentTime * Units::T,
           currentRadius * Units::L * 1.0e6,
           Stats.TotalCollisions,
           Stats.TotalCellCrossings,
           Stats.TotalWallCollisions,
           Stats.TotalConeBoundaryCollisions,
           Stats.KineticEnergyPerAtom / MaxAtoms * Units::Energy * 6.24150974e15,
           me,
           Stats.FusionRate);

    PerformanceClock = clock();

    writer.String("snap_shot") , writer.Int(SnapshotCount);
    writer.String("t")         , writer.Double(CurrentTime * Units::T);
    writer.String("dt")        , writer.Double((CurrentTime - LastSnapshotTime) * Units::T);
    writer.String("r")         , writer.Double(currentRadius * Units::L);
    writer.String("real_time") , writer.Double(delta_t);
    writer.String("cone_angle"), writer.Double(Info.ConeAngle * 180 / Constants::PI);

    // Write each bin of the simulation statistics.
    writer.String("bins");
    Stats.WriteStatistics(writer);

    // More performance data
    writer.String("TotalCollisions")    , writer.Uint64(Stats.TotalCollisions);
    writer.String("TotalCellCrossings") , writer.Uint64(Stats.TotalCellCrossings);
    writer.String("TotalWallCollisions"), writer.Uint64(Stats.TotalWallCollisions);

    writer.String("TotalConeBoundaryCollisions");
    writer.Uint64(Stats.TotalConeBoundaryCollisions);

    writer.EndObject(); // End the JSON structure
    os.Flush();
    fclose(outfile);

    LastSnapshotTime = CurrentTime;
}

/*
 *void Simulation::Read(const ptree& tree)
 *{
 *    using std::vector;
 *
 *    Info.Read(tree.get_child("Info"));
 *    RPKEquation.Read(tree.get_child("RPKEquation"));
 *
 *    CurrentTime         = tree.get<double>("CurrentTime");
 *    TotalCollisions     = tree.get<unsigned long long>("TotalCollisions");
 *    TotalCellCrossings  = tree.get<unsigned long long>("TotalCellCrossings");
 *    TotalWallCollisions = tree.get<unsigned long long>("TotalWallCollisions");
 *    TotalConeBoundaryCollisions =
 *        tree.get<unsigned long long>("TotalConeBoundaryCollisions");
 *    FusionRate = tree.get<unsigned long long>("FusionRate");
 *
 *    print = tree.get<bool>("print");
 *
 *    ptree particles = tree.get_child("Particles");
 *    Particles.resize(particles.size());
 *
 *    ReadArray<vector<Particle> >(particles, Particles.begin(), Particles.end());
 *}
 */

