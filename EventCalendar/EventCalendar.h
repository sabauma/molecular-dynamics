/**
 * @author Spenser Bauman
 *
 * Interface for the Event Calendar class.
 *
 * A node in the event calendar. Each node represents one event in the tree.
 * It is structured using a binary search tree ordered by time. Each event
 * involves two objects in the system. For collisions between two particles,
 * it is the integer IDs of the two particles. For cell crossing events, the
 * first is the ID of the particle and the second is the ID of the cell wall
 * being crossed.
 *
 * Each node uses six pointers to allow for quick updating of the calendar.
 * Left and Right track the children as a typical BST would. Parent is the
 * node's parent in the tree (where NULL indicates we are at the root).
 *
 * CircleAL/CircleAR are the left and right pointers in the circular list of
 * events involving object A. CircleBL/CircleBR do the same for object B.
 *
 * An important invariant of this code is the invariance of the Type field for
 * each <code>CalendarNode</code>. Once the Type field is set, just after
 * allocation, it is not changed until <code>FreeNode</code> is called. In
 * effect, the <code>Type</code> field is const for the entire time each node
 * is validly allocated.
 */

#ifndef __EVENT_CALENDAR_EVENT_CALENDAR_H__
#define __EVENT_CALENDAR_EVENT_CALENDAR_H__

#include <vector>
#include "EventType.h"
#include "../Constants/Constants.h"

class EventCalendar
{
private:

    /**
     * Representation of the nodes in the tree. Contains all of the
     * information about each event, along with all of the pointers needed to
     * manage the tree and the list of associated events.
     */
    struct CalendarNode
    {
        /** Tree structure */
        CalendarNode* Parent;   /**< The parent of the current node.                  */
        CalendarNode* Left;     /**< The left child of the current node.              */
        CalendarNode* Right;    /**< The right child of the current node.             */
        CalendarNode* CircleAL; /**< Left pointer in the circular linked list for A.  */
        CalendarNode* CircleAR; /**< Right pointer in the circular linked list for A. */
        CalendarNode* CircleBL; /**< Left pointer in the circular linked list for B.  */
        CalendarNode* CircleBR; /**< Right pointer in the circular linked list for B. */

        /** Event information. */
        double Time;    /**< The time of the event. Used to schedule the event. */
        int ObjectA;    /**< The ID of the first object in the event.           */
        int ObjectB;    /**< The ID of the second object in the event.          */
        EventType Type; /**< The Type of the event.                             */


        /** Constructors */
        CalendarNode();
        CalendarNode(const CalendarNode& other);

        /** Insert a new node into the tree. */
        void Insert(CalendarNode* ins);

        /**
         * Deletes a node in the tree. Ensures that the binary tree property is
         * preserved. It also ensures that the linked list structure is preserved.
         *
         * There are four cases that need to be covered to ensure the tree is
         * restructured properly. In the worst case, the node to be deleted has
         * two children and extra work needs to be performed, as we cannot just
         * overwrite node contents as in a naive approach.
         *
         * @return A pointer to the node that replaces the one being deleted in the
         *         tree. If we delete the root of the tree, this should
         *         become the new root node.
         */
        CalendarNode* Delete();

        /**
         * Computes the pair of values of the min and max distance to a leaf
         * from the calling node.
         */
        std::pair<unsigned long, unsigned long> MinMax() const;
    };

    std::vector<CalendarNode> ListAnchor;  // Entry point for circularly linked lists

    CalendarNode* FreeList;     /**< List of free calendar nodes.        */
    CalendarNode* Root;         /**< Root of the tree.                   */

    bool HasStepped;            /**< True if an event has been executed. */

    double CurrentEventTime;    /**< Time of the current event           */
    int CurrentEventObjectA;    /**< First object of the current event   */
    int CurrentEventObjectB;    /**< Second object of the current event  */
    EventType CurrentEventType; /**< Type of the current event           */

    /** Allocates one node from the free list. */
    CalendarNode* AllocateNode();

    /**
     * Returns the given node to the pool. This invalidates all of the pointers in
     * the node, resulting in an error if it referenced after being freed.
     *
     * @param node The node to be returned to the pool.
     */
    void FreeNode(CalendarNode* node);

    /**
     * Safely deletes a node from the tree. This performs the delete operation but
     * also prevents the root node from becoming invalid. If the node that is
     * deleted is the root, we replace the root with the node returned by the
     * delete function, which is the node that takes the place of the
     * deleted one.
     */
    CalendarNode* SafeDelete(CalendarNode* node);

    /**
     * Checks that a given node is in the tree by checking its pointer values.
     * This function does not actually traverse the tree to determine if it
     * can be reached, due to the computational cost of doing so.
     */
    bool NodeInTree(CalendarNode const * const node);

public:

    /**
     * Initialize the calendar for a given number of particles.
     *
     * @param particleCount The number of particles in the system that the
     *                      event calendar is tracking.
     */
    EventCalendar(const unsigned particleCount);

    /** Clean up the dynamically allocated memory. */
    ~EventCalendar();

    /** Schedule an event with the calendar. */
    void ScheduleEvent(const double time, const EventType type,
                       const int idA=Constants::NULL_PARTICLE,
                       const int idB=Constants::NULL_PARTICLE);

    /**
     * Advance the calendar one step and put the next event into the current
     * event place holders.
     */
    void NextEvent();

    /**
     * Test whether the calendar has more events in it.
     *
     * @return True if there are more events in the calendar, otherwise false.
     */
    bool HasMoreEvents();

    /**
     * Empties all events from the event calendar and returns the allocated
     * tree nodes to the free list. The memory is kept, as it is expected that
     * it will be used again.
     */
    void ClearCalendar();

    /**
     * @return The time of the current event in the simulation.
     */
    double GetCurrentTime() const;

    /**
     * @return The id of the first object in the simulation.
     */
    int GetCurrentObjectA() const;

    /**
     * @return The id of the second object in the simulation.
     */
    int GetCurrentObjectB() const;

    /**
     * Computes the ratio of the longest path to a leaf over the shortest
     * path to a leaf. This is for useful diagnostics about tree balancing.
     */
    double BalanceRatio() const;

    /**
     * @return The type of the current event in the simulation.
     */
    EventType GetCurrentEventType() const;
};

#endif  /* __EVENT_CALENDAR_H__ */

