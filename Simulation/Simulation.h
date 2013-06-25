/**
 * Definition of the simulation class. This will hold the relevant information
 * for the simulation. This includes the event tree, the cell partitions, and the
 * particles in the system. This provides an interface to execute the
 * simulation along with methods to access the data.
 *
 * The interface is provided through a series of iterators that allow access
 * to the particles with contain all of the information about the system.
 *
 * The simulation that this code describes is a complete rewrite of the work of
 * Alexander Bass. His original implementation of an event driven molecular
 * dynamics simulation of bubble cavitation is described in his Ph.D.
 * dissertation "Molecular Dynamics Simulation of Sonoluminescence."
 *
 * This describes an MD simulation which leverages symmetry to reduce the
 * simulation space to a narrow conic region of the original bubble. This
 * makes the simulation computationally feasible. This is an attempt at a
 * production quality version of the simulation that is more easily verifiable
 * and easier to modify.
 *
 * @author Spenser Bauman
 */

#ifndef __SIMULATION_SIMULATION_H__
#define __SIMULATION_SIMULATION_H__

#include "../CellList/CellListFast.h"
#include "../EventCalendar/EventCalendar.h"
#include "../EventCalendar/EventType.h"
#include "../Particle/Particle.h"
#include "../RPK/RPK.h"
#include "../Vector/Vector.hpp"

#include "./CollisionRecord.h"
#include "./Statistics.h"

#include <math.h>
#include <vector>


/**
 * Function object to regulate the ScheduleParticleCollisions function. This
 * will cause the function to examine all particles that are in adjacent
 * cells.
 */
class AllParticles
{
public:

    /**
     * Build the object. Nothing to be done really.
     */
    explicit AllParticles(const int)
    { }

    /**
     * Test true for all given particle indices.
     * @param pa The index to test against the upper limit.
     */
    inline bool operator()(const int)
    {
        return true;
    }
};

/**
 * Function object to regulate the ScheduleParticleCollisions function. This
 * class will cause the function to examine all particles with index values
 * less than the index this function was constructed with.
 */
class LessThan
{
private:
    const int Upper;    // Strict upper bound of particles to accept.

public:
    /**
     * Build the object with the upper bound.
     */
    explicit LessThan(const int i) : Upper(i)
    { }

    /**
     * Test for acceptance.
     * @param pa The index to test against the upper limit.
     */
    inline bool operator()(const int pa)
    {
        return pa < Upper;
    }
};

/**
 * Function object to regulate the ScheduleParticleCollisions function. This
 * will cause the function to examine all particles that are in adjacent cells
 * that do not have index values corresponding to the one this object was
 * constructed with.
 */
class Skipping
{
private:
    const int Skip;     // Index value to skip.

public:
    /**
     * Create the object with the value to skip.
     * @param s The value to skip.
     */
    explicit Skipping(const int s) : Skip(s)
    { }

    /**
     * Reject only when the the Skip value is given.
     * @param pa The particle index to test against.
     */
    inline bool operator()(const int pa)
    {
        return pa != Skip;
    }
};

/**
 * The Simulation class is the environment in which the simulation occurs. It
 * internally tracks all of the relevant information to the system and
 * provides an interface for accessing this information for analysis purposes.
 */
class Simulation
{
private:

    // The number of radial bins to compute statistic for when writing a snapshot.
    static const int MAX_BINS = 48;

    // The number of collisions to be recorded.
    static const size_t COLLISION_RECORDS = 100000;

    typedef std::unique_ptr<EventCalendar> CalendarPtr;

    InfoStruct Info;            /**< Simulation information.                   */

    CellList Partitions;  /**< The cells for nearest neighbor search. */
    CalendarPtr Calendar; /**< The events occurring in the system.    */

    std::vector<Particle> Particles;   /**< The particles in the system.       */
    std::vector<double>   ConeTime;    /**< Cached cell crossing times         */
    std::vector<double>   WallTime;    /**< Cached bubble wall crossing times. */
    double                CurrentTime; /**< Cached current simulation time.    */

    int ParticleCount; /**< # of particles currently in the system. */
    int MaxParticles;  /**< max # of particle after dissociation. */

    int SnapshotCount;                 /**< The number of snapshots taken. */

    DoubleVector Extent;    /**< The size of the simulation space on each axis. */
    DoubleVector Offset;    /**< Offset of the bounding box w.r.t. the origin.  */
    IntVector    Cells;     /**< Number of cells along each axis.               */
    IntVector    InitUCell; /**< Number of cells used to initialize the system. */

    RPK        RPKEquation;     /**< The RPK equation representing the bubble. */
    double     LastRebinRadius; /**< Radius when last rebin happened.          */
    double     LastVelocity;    /**< Velocity at the last system update */
    double     LastSnapshotTime; /**< Time of the last snapshot */
    bool       Rebounded;       /**< Flag indicating that the simulation has rebounded */

    Statistics Stats;

    clock_t PerformanceClock;

    /**
     * Predict all events based on the given previous event. This is
     * a optimized version of all events that can elide some
     * of the work.
     *
     * The template parameter is used to control the collision prediction
     * mode. When predicting collisions, an object of type
     * <code>Controller</code> is created with <code>nb</code> as a parameter
     * to its constructor.
     *
     * @param na The particle to predict events for as the particle's index.
     * @param ev The event that just occurred, which is used to elide some of
     *           the prediction calculations based on the invariance of some
     *           event types w.r.t. the occurrence of other events.
     */
    template <class Controller>
    inline void PredictEvents(const int na, const int nb, const EventType ev);

    /**
     * Predict the cell crossing event for a single particle. This determines
     * the crossing time and stores the information in the given DoubleVector.
     *
     * The return value is the integer representation of the dimension with
     * the nearest time to intersection.
     *
     * @param na The particle to find crossing times for.
     * @param times Storage for the relative intersection time for each axis.
     * @return The axis with the nearest intersection time.
     */
    int CellCrossingTimes(const int na, DoubleVector& times);

    /**
     * Predict the time until a cone collision occurs for the given particle.
     * This function does not alter the ConeTime array at all, it only
     * performs the time calculation. If the actual event is being predicted,
     * the calculated time needs to be inserted into the ConeTime array to
     * allow efficient calculation.
     *
     * @param na The particle to predict the collision time for.
     * @return The time until the collision occurs in simulation units.
     */
    double ConeCollisionTime(const int na);

    /**
     * Computes the time until a bubble wall collision occurs for the
     * particle.
     *
     * @param na The particle to predict the collision time for.
     * @return The time until the collision occus in simulation units.
     */
    double WallCollisionTime(const int na);

    /**
     * The time until a particle collision occurs. This will check collision
     * times with all particle in adjacent cells and return the nearest
     * collision time. If no collision occurs, Constants::NEVER is returned.
     * The second parameter will be used to store the resulting index of
     * the particle colliding with this one.
     *
     * @param na The particle to predict the collision time for.
     * @param con The function object that controls the particle collision
     *            scheduling loop. This object is simply a class with the
     *            function call operator overridden to act as a preciate.
     *            Only particles which this predicate tests true for will be
     *            examined for intersection with the particle corresponding
     *            to na.
     * @param low The low bounds on each axis for the search.
     * @param high The high bounds on each axis for the search.
     */
    template <class Controller>
    void ScheduleParticleCollisions(const int na, Controller& con,
                                    const IntVector& low,
                                    const IntVector& high);

    /**
     * Functions the same as the previous method, but does not take a
     * pair of bounds on which to search through. This could call the
     * other version with default values, but performance is a concern
     * and this will reduce potential copying.
     */
    template <class Controller>
    void ScheduleParticleCollisions(const int na, Controller& con);

    /**
     * Given a pair of particles, this computes their intersection time in
     * simulation seconds.
     *
     * @param pa The first particle of the potential collision.
     * @param pb The second particle of the potential collision.
     * @return The time until the collsion in simulation seconds. If there is
     *         no collision, <code>Constants::NEVER</code> is returned.
     */
    double ComputeIntersectionTime(const Particle& pa, const Particle& pb);

    /**
     * Processes the cell wall crossing event. This function assumes that
     * a cell wall crossing is the event that just occurred in the simulation
     * and that <code>CurrentTime</code> was advanced to the time of that
     * event.
     */
    void ProcessCellCrossing();

    /**
     * Handles bubble wall collision events. This function assumes that
     * a bubble wall impact occurs at the time specified in the
     * <code>CurrentTime</code> variable.
     */
    void ProcessBubbleWall();

    /**
     * Handles particle-particle collision events. This function assumes that
     * the current event is a collision event and that it occurs at the time
     * specified by the <code>CurrentTime</code> variable.
     */
    void ProcessCollision();

    /**
     * Used in the <code>ProcessCollision</code> method. This method takes
     * the index of a particle and dissociates it, if that is possible.
     *
     * This method handles all the physics involved in molecular dissociation
     * as well as the details of allocating new particles.
     * The particle at the provided index is converted into one of the two
     * particles being produced, and a new particle is allocated for the
     * second particle.
     *
     * @param na The index of the particle to dissociate.
     */
    void DissociateParticle(const int na);

    /**
     * Processes a cone collision event. This function assumes that the
     * current event is a cone intersection event and that it occurs at the
     * time specified by the <code>CurrentTime</code> variable.
     */
    void ProcessCone();

    /**
     * Performs the system update event. This function assumes that the
     * current event is a system update event. Unlike the other events, there
     * is no negative effect of performing this operation, other than a drop
     * in performance. As such, it should only be used when needed, though it
     * will cause no harm to the simulation (in terms of simulation
     * correctness).
     */
    void ProcessUpdateSystem();

    /**
     * Has the RPKEquation driver read the pressure gauge. Afterwards, a new
     * one is scheduled at the appropriate time. However, it is important to
     * note that anything that clears the calendar must reschedule this event.
     * In general, it should be sufficient ot call the
     * <code>PopulateCalendar</code> method.
     */
    void ProcessReadPGauge();

    /**
     * Computes the square of the minimum approach distance between two atoms.
     *
     * @param vv The dot product (with itself) of the difference in
     *           velocity vectors between the two particles.
     * @param p1 The first particle.
     * @param p2 The second particle.
     * @return The square of the minimum approach distance between the two
     *         atoms.
     */
    double Diameter2(const double vv, const Particle& p1, const Particle& p2);

    /**
     * Perform a rebinning of the simulation. This is makes use of the
     * values of Extent, Offset, and Cells to determine the location of each
     * particle in the Partitions structure. Calling this will invalidate any
     * coordinates and iterators obtained from Partitions before its use.
     */
    void BinSystem();

    /**
     * Clears out the partitions and the event calendar. The system is
     * then rebinned and the event calendar is repopulated.
     */
    void StartRun();

    /**
     * Counts the number of sites in the grid which can hold a particle. One
     * particle is placed on the face of each grid. This ensures that particles
     * are only placed inside of the cone and will reside within the bubble
     * given the initial bubble radius.
     */
    int ComputeParticlesInSystem() const;

    /**
     * Initializes the positions and velocities of all the particles in the
     * system as well as the types of each particle in the system.
     */
    void InitializeParticles();


    /**
     * Detect events for all particles and schedule them.
     */
    void PopulateCalendar();

    /**
     * Writes snapshot information to a file. This is used by the update
     * system function to determine when another update event should occur.
     */
    void WriteSnapshot();


    /**
     * Write the collision energy information to disk. This produces a file
     * in CSV format in the current directory named
     * "CollisionEnergy########.csv".
     */
    void WriteCollisionEnergy();

public:

    friend std::ostream& operator<<(std::ostream& stream, const Simulation& sim);

    /**
     * A templatized version of the constructor that takes a particle factory
     * as a type. This makes it fairly easy to swap out particle generators
     * for testing purposes.
     *
     * @param info The info struct providing needed constants and other
     *             configurable parameters of the simulation.
     * @param n The number of particle to generate.
     */
    explicit Simulation(const InfoStruct& info);

    /**
     * Prints all the events in the calendar to standard out if we are in
     * debug mode.
     */
    void PurgeEvents();

    /**
     * Returns whether or not the bubble wall has rebounded yet. This is
     * useful for determining when to stop advancing the simulation.
     *
     * @returns True if the bubble wall has rebounded, otherwise false.
     */
    bool HasRebounded() const;

    /**
     * Advance the simulation by one step. This mean executing one event and
     * then updating the calendar and cell list appropriately.
     */
    EventType NextStep();

    /**
     * Beginning iterator for the particles in the system. This is just a
     * standard iterator over a vector of particles.
     */
    std::vector<Particle>::const_iterator begin();

    /**
     * Ending iterator used to determine the stopping criteria for the
     * particle iterator.
     */
    std::vector<Particle>::const_iterator end();

    /**
     * @return The current time in simulation units. Mostly just to view the
     * progress of the simulation.
     */
    double GetCurrentTime() const;

    /**
     * Get the collision count from the simulation.
     */
    unsigned long long GetCollisionCount() const;

    /**
     * Get the number of cell crossing events so far.
     */
    unsigned long long GetCellCrossingCount() const;

    /**
     * Get the number of cell wall events so far.
     */
    unsigned long long GetWallCollisionCount() const;

    /**
     * Get the number of cone boundary events so far.
     */
    unsigned long long GetConeBoundaryCount() const;

    /**
     * Serialize the simulation to the writer being given.
     */
    template <typename Writer>
    void Serialize(Writer& writer) const;

    // void Read(const boost::property_tree::ptree& tree);
};

template <typename Controller>
void Simulation::PredictEvents(const int na, const int nb, const EventType ev)
{
    Controller con(nb);
    Particle& pa = Particles[na];
    pa.UpdateParticle(this->CurrentTime);

    DoubleVector cellTime(0.0, 0.0, 0.0);

    const int evCode  = this->CellCrossingTimes(na, cellTime);
    cellTime[evCode] += this->CurrentTime;

    // Do not repredict cone wall and bubble wall events after a
    // cell crossing event. The previous values are not invalidated by
    // these events.
    if (ev != CellCrossingEvent)
    {
        this->ConeTime[na] = this->ConeCollisionTime(na) + this->CurrentTime;

        // Do not repredict Bubble wall crossings after cone
        // events, as cone wall events, supposedly, do not alter the
        // radial velocity of the particles
        if (ev != ConeWallEvent)
        {
            this->WallTime[na] = this->WallCollisionTime(na);
        }
    }

    // Find the nearest event for this particle and schedule.
    if (cellTime[evCode] < this->ConeTime[na])
    {
        if (cellTime[evCode] < this->WallTime[na])
        {
            // Schedule a cell crossing event.
            this->Calendar->ScheduleEvent(
                    cellTime[evCode], CellCrossingEvent, na, evCode);
        }
        else
        {
            // Schedule a bubble wall event. The second particle in bubble
            // wall events is ignored.
            this->Calendar->ScheduleEvent(
                    this->WallTime[na], BubbleWallEvent, na, -1);
        }
    }
    else
    {
        if (this->ConeTime[na] < this->WallTime[na])
        {
            // Schedule a cone wall event. The second particle in cone
            // wall events is ignored.
            this->Calendar->ScheduleEvent(
                    this->ConeTime[na], ConeWallEvent, na, -1);
        }
        else
        {
            // Schedule a bubble wall event. The second particle in bubble
            // wall events is ignored.
            this->Calendar->ScheduleEvent(
                    this->WallTime[na], BubbleWallEvent, na, -1);
        }
    }


    // If the last event was a cell crossing, we can limit the search to
    // those in the direction we are heading, as we will have retained
    // the collisions from particles within the current and previous
    // cell.
    IntVector low(-1, -1, -1);
    IntVector high(1, 1, 1);

    // If a cell crossing just occurred, we can optimize by providing
    // more restricted bounds on the neighborhood search.
    if (ev == CellCrossingEvent)
    {
        if (pa.Velocity[nb] >= 0.0)
        {
            low[nb] = 1;
        }
        else
        {
            high[nb] = -1;
        }
    }

    this->ScheduleParticleCollisions<Controller>(na, con, low, high);
}

template <typename Controller>
void Simulation::ScheduleParticleCollisions(const int na, Controller& con)
{
    static const IntVector low(-1, -1, -1);
    static const IntVector high(1, 1, 1);
    this->ScheduleParticleCollisions<Controller>(na, con, low, high);
}

/*
 * Alright, this seems like it deserves some explanation, as it is somewhat of
 * a complex design.
 *
 * This is to cheat the need for having multiple versions of this function or
 * having a convoluted data encoding in the parameters/global variables.
 * Controller is a function object with a function call operator
 * that acts as a predicate function.
 *
 * If the predicate tests true, then the given particle is tested for
 * collision with the particle designated by <code>na</code>.
 *
 * In this sense, the parameter <code>con</code> acts as a closure containing
 * all the data needed to determine whether or not to use this particle.
 */
template <typename Controller>
void Simulation::ScheduleParticleCollisions(
        const int na, Controller& con,
        const IntVector& low,
        const IntVector& high)
{
    assert(na >= 0);
    assert((unsigned) na < Particles.size());

    Particle& pa = Particles[na];

    CellList::NeighborCellIterator it = Partitions.BeginNeighborCell(na, low, high);

    for (; it.HasNext(); ++it)
    {
        if (*it != na && con(*it))
        {
            Particle& pb = Particles[*it];
            pb.UpdateParticle(this->CurrentTime);

            const double time = ComputeIntersectionTime(pa, pb);

            if (time != Constants::NEVER)
            {
                this->Calendar->ScheduleEvent(
                        this->CurrentTime + time, CollisionEvent, na, *it);
            }
        }
    }
}

template <typename Writer>
void Simulation::Serialize(Writer& writer) const
{
    writer.String("Info")        , Info.Serialize(writer);
    writer.String("RPKEquation") , RPKEquation.Serialize(writer);

    writer.String("CurrentTime")         , writer.Double(CurrentTime);
    writer.String("TotalCollisions")     , writer.UInt64(Stats.TotalCollisions);
    writer.String("TotalCellCrossings")  , writer.UInt64(Stats.TotalCellCrossings);
    writer.String("TotalWallCollisions") , writer.UInt64(Stats.TotalWallCollisions);

    writer.String("TotalConeBoundaryCollisions");
    writer.UInt64(Stats.TotalConeBoundaryCollisions);

    writer.String("FusionRate"), writer.UInt64(Stats.FusionRate);

    // Serialize the array of particles to the document
    writer.String("Particles");
    writer.StartArray();
    std::for_each(Particles.begin(), Particles.end(), [&writer](const Particle& p)
            {
                p.Serialize(writer);
            });
    writer.EndArray();
}

#endif /* __SIMULATION_H__ */
