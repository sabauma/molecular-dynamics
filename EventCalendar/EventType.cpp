/**
 * @author Spenser Bauman
 * @section Description
 *
 * Implementation of functionality for the event types.
 */

#include <assert.h>
#include <boost/lexical_cast.hpp>
#include <string>

#include "EventType.h"
#include "../Exceptions/SimulationException.h"

using boost::lexical_cast;

/** Strings for the names and possible errors. */
static const std::string BUBBLE_WALL_STR = std::string("Bubble Wall Event");
static const std::string COLLISION_STR   = std::string("Collision Event");
static const std::string CONE_WALL_STR   = std::string("Cone Wall Event");
static const std::string CROSSING_STR    = std::string("Crossing Event");
static const std::string UPDATE_STR      = std::string("Update Event");

static const std::string UNKNOWN_EVENT =
    std::string("ERROR: unknown event type at line ");

#define LOOKUP_ERROR \
    (UNKNOWN_EVENT + lexical_cast<std::string>(__LINE__) + " in " + __func__)

static const int ASSOCIATED_OBJECTS[] = { 1, 1, 1, 2, 0 };

/**
 * Function for translating an event type to an appropriate string represting
 * the event.
 *
 * @param type An event.
 * @retval A string representation of the given event.
 */
std::string
getEventName(const EventType type)
{
    GUARD_TYPE(type);

    switch (type)
    {
        case UpdateSystemEvent:
            return UPDATE_STR;
        case BubbleWallEvent:
            return BUBBLE_WALL_STR;
        case CellCrossingEvent:
            return CROSSING_STR;
        case ConeWallEvent:
            return CONE_WALL_STR;
        case CollisionEvent:
            return COLLISION_STR;
        default:
            throw SimulationException(LOOKUP_ERROR);
    }
}

/**
 * Produces the number of objects associated with the event.
 * This is always in the range [0,2].
 *
 * @param type The type of the event.
 * @retval The number of particles associated with the given event type.
 */
int
associatedObjects(const EventType type)
{
    GUARD_TYPE(type);
    assert(type >= CellCrossingEvent && type <= UpdateSystemEvent);
    return ASSOCIATED_OBJECTS[(int) type];

    //switch (type)
    //{
        //case UpdateSystemEvent:
        //case ReadPressureGaugeEvent:
            //return 0;
        //case BubbleWallEvent:
        //case CellCrossingEvent:
        //case ConeWallEvent:
            //return 1;
        //case CollisionEvent:
            //return 2;
        //default:
            //throw SimulationException(LOOKUP_ERROR);
    //}
}

bool
invalidatesAssociated(const EventType type)
{
    return type != CellCrossingEvent;
}
