/**
 * @author Spenser Bauman
 * @section Description
 *
 * Enumeration representing the different types of events that can occur in
 * a simulation. These are used to schedule the events with the calendar as
 * well as to act upon them when their time is due.
 */

#ifndef __EVENT_TYPE_H__
#define __EVENT_TYPE_H__

#include <string>

#define GUARD_TYPE(typ) (assert((typ) == UpdateSystemEvent || \
                                (typ) == BubbleWallEvent   || \
                                (typ) == CellCrossingEvent || \
                                (typ) == ConeWallEvent     || \
                                (typ) == CollisionEvent))


/**
 * Enumeration of the various event types in the system. The names correspond
 * to the events exactly as expected.
 */
enum EventType
{
    // Associated with one particle
    CellCrossingEvent      = 0,
    BubbleWallEvent        = 1,
    ConeWallEvent          = 2,

    // Associated with two particles
    CollisionEvent         = 3,

    // Associated with no particles
    UpdateSystemEvent      = 4,
};

/**
 * Function for translating an event type to an appropriate string represting
 * the event.
 *
 * @param type An event.
 * @retval A string representation of the given event.
 */
std::string getEventName(const EventType type);

/**
 * Produces the number of objects associated with the event.
 * This is always in the range [0,2].
 *
 * @param type The type of the event.
 * @retval The number of particles associated with the given event type.
 */
int associatedObjects(const EventType type);

/**
 * Check function to determine whether or not the occurrence of a given event
 * will invalidate ones already queued. Currently, the only such event type
 * is the <code>CellCrossingEvent</code>. Technically speaking,
 * the <code>UpdateSystemEvent</code> does not invalidate its associated
 * particles, but it does not have any, so logically, either response works,
 * and this reduces it to one comparision operation.
 */
bool
invalidatesAssociated(const EventType type);

#endif
