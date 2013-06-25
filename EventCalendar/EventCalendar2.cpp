/**
 * @author Spenser Bauman
 * @section Description
 *
 * Implementation of the EventCalendar class. This is a representation of an
 * event calendar similar to the one found in "The Art of Molecular Dynamics"
 * by Rapaport.
 *
 * The data structure is a highly modified binary search tree. Each node is
 * part of two doubly linked circular lists which are used to associate all
 * events that involve the same particles.
 *
 * Nodes for the tree are allocated dynamically as needed, but the space
 * allocated is stored in a pool of nodes for later use. This is allows
 * the calendar to allocate only as much space as it needs and maintains
 * good performance by not repeatedly calling <code>new</code>
 * and <code>delete</code>. The calendar will quickly allocate the maximum
 * number of nodes it needs and not need to allocate more after that, as
 * the total number of nodes in use is usually some constant factor of the
 * number of particles in the system.
 *
 * Another thing to note. When deallocting nodes, it is imperative that these
 * nodes not be referenced anywhere in the tree or already be in the free
 * list. This will cause bad things to happen.
 */

#include <assert.h>
#include <cstdlib>
#include <limits>
#include <stdio.h>
#include <vector>

#include "EventCalendar.h"
#include "EventType.h"

using namespace std;

EventCalendar::CalendarNode::CalendarNode() :
    Parent(NULL),
    Left(NULL),
    Right(NULL),
    CircleAL(NULL),
    CircleAR(NULL),
    CircleBL(NULL),
    CircleBR(NULL),
    Time(0.0),
    ObjectA(0),
    ObjectB(0),
    Type(CollisionEvent)
{ }

EventCalendar::CalendarNode::CalendarNode(const CalendarNode& other) :
    Parent(other.Parent),
    Left(other.Left),
    Right(other.Right),
    CircleAL(other.CircleAL),
    CircleAR(other.CircleAR),
    CircleBL(other.CircleBL),
    CircleBR(other.CircleBR),
    Time(other.Time),
    ObjectA(other.ObjectA),
    ObjectB(other.ObjectB)
{ }

/**
 * Insert a new node into the tree. Ensures the binary tree property is
 * satisified based on the Time field of the node. If the node given is NULL,
 * nothing is done, but this makes other operations simpler.
 *
 * @param ins The node to be inserted into the tree. Note, the pointer given
 *            cannot be NULL.
 */
void
EventCalendar::CalendarNode::Insert(CalendarNode* ins)
{
    EventCalendar::CalendarNode* curr = NULL; // Current node
    EventCalendar::CalendarNode* next = this; // Next node

    // Descend to a leaf
    while (NULL != next)
    {
        curr = next;
        if (ins->Time < next->Time)
        {
            next = next->Left;
        }
        else
        {
            next = next->Right;
        }
    }

    if (ins->Time < curr->Time) // Insert into the left node
    {
        curr->Left = ins;
    }
    else                        // Insert into the right node
    {
        curr->Right = ins;
    }

    // Setup pointers to surrounding node
    ins->Parent = curr;
}

std::pair<unsigned long, unsigned long>
EventCalendar::CalendarNode::MinMax() const
{
    std::pair<unsigned long, unsigned long> left(0, 0);
    std::pair<unsigned long, unsigned long> right(0, 0);
    std::pair<unsigned long, unsigned long> res(0,0);

    if (NULL != this->Left)
    {
        left = this->Left->MinMax();
    }

    if (NULL != this->Right)
    {
        right = this->Right->MinMax();
    }

    res.first  = std::min(left.first, right.first) + 1;
    res.second = std::max(left.second, right.second) + 1;

    return res;
}

EventCalendar::CalendarNode*
EventCalendar::CalendarNode::Delete()
{
    const bool atRoot = this->Parent == NULL;

    EventCalendar::CalendarNode** fromParent = NULL;
    EventCalendar::CalendarNode* retval = NULL;

    // Used to update the pointer from the parent node to the current node.
    if (!atRoot)
    {
        fromParent = this->Parent->Left == this
                     ? &this->Parent->Left
                     : &this->Parent->Right;
    }

    if (NULL == this->Left)
    {
        // Left child is empty
        retval = this->Right;   // Right child is new root
    }
    else if (NULL == this->Right)
    {
        //Right child is empty
        retval = this->Left;    // Left child is new root
    }
    else
    {
        // Neither child is empty
        retval = this->Left;

        // Insert left->right grand child into the right tree. This clears up
        // the left->right grand child position for the current right child.
        if (NULL != this->Left->Right)
        {
            this->Right->Insert(this->Left->Right);
        }

        this->Left->Right   = this->Right;
        this->Right->Parent = this->Left;
    }

    if (!atRoot)
    {
        *fromParent = retval;
    }

    // Set the parent pointer, if we have a new child to take this one's
    // place.
    if (NULL != retval)
    {
        retval->Parent = this->Parent;
    }

    this->Left = this->Right = this->Parent = NULL;

    // Only fix up the list if we are not at a cell crossing event.  Cell
    // crossing events anchor the list, so they should never be removed from
    // the list. Otherwise, base it off the number of objects associated with
    // the event type.
    if (this->Type != CellCrossingEvent)
    {
        // NOTE: This uses the natural fall through behaviour of the switch
        // statement. Notice the lack of a break for case 2.
        switch (associatedObjects(this->Type))
        {
            case 2:
                this->CircleBL->CircleBR = this->CircleBR;
                this->CircleBR->CircleBL = this->CircleBL;
            case 1:
                this->CircleAL->CircleAR = this->CircleAR;
                this->CircleAR->CircleAL = this->CircleAL;
                break;
        }

        this->CircleAL = this->CircleAR = this->CircleBL = this->CircleBR = NULL;
    }

    return retval;
}

EventCalendar::EventCalendar(const unsigned int particleCount)
    : Root(NULL),
      HasStepped(false),
      CurrentEventTime(0.0),
      CurrentEventObjectA(0),
      CurrentEventObjectB(0),
      CurrentEventType(CollisionEvent)
{
    const size_t PREALLOC_SIZE = 2 * particleCount;
    EventCalendar::CalendarNode* prealloc_nodes = NULL;
    // Reserve space
    ListAnchor.resize(particleCount);

    std::for_each(ListAnchor.begin(), ListAnchor.end(),
            [](CalendarNode& node)
            {
                node.Type     = CellCrossingEvent;
                node.CircleAL = node.CircleAR
                              = node.CircleBL
                              = node.CircleBR
                              = &node;
            });

    FreeList = NULL;

    // Allocate a large number of nodes in the hope that they are contiguous.
    // The hope is to improve cache behaviour.
    for (size_t i = 0; i < PREALLOC_SIZE; ++i)
    {
        EventCalendar::CalendarNode* next = this->AllocateNode();
        next->Right = prealloc_nodes;
        prealloc_nodes = next;
    }

    // Return nodes to the free list
    while (prealloc_nodes != NULL)
    {
        EventCalendar::CalendarNode* curr = prealloc_nodes;
        prealloc_nodes = prealloc_nodes->Right;
        this->FreeNode(curr);
    }
}

EventCalendar::~EventCalendar()
{
    // Return all dynamically allocated nodes to the free list.
    this->ClearCalendar();

    CalendarNode* curr = NULL;
    CalendarNode* next = FreeList;

    // Traverse the free list and delete each element.
    while (next != NULL)
    {
        curr = next;
        next = next->Right;
        delete curr;
    }
}

void EventCalendar::ClearCalendar()
{
    // Delete all events from the tree, returning all the nodes to
    // the free list.
    while (this->HasMoreEvents())
    {
        this->NextEvent();
    }

    std::for_each(ListAnchor.begin(), ListAnchor.end(),
            [](CalendarNode& node)
            {
                node.Type     = CellCrossingEvent;
                node.CircleAL = node.CircleAR
                              = node.CircleBL
                              = node.CircleBR
                              = &node;
            });

    CurrentEventTime    = 0.0;
    CurrentEventObjectA = -1;
    CurrentEventObjectB = -1;
    CurrentEventType    = UpdateSystemEvent;
}

/**
 * Allocates a new node from the pool. This will fail if the pool has been
 * completely allocated.
 *
 * @retval A pointer to the newly allocated node. All of the pointer fields
 *         are set to NULL.
 */
EventCalendar::CalendarNode*
EventCalendar::AllocateNode()
{
    // If the free list is empty, allocate a new node and put it on the
    // free list.
    if (NULL == FreeList)
    {
        return new CalendarNode();
    }

    CalendarNode* next = FreeList;
    FreeList = FreeList->Right;

    // next->Left  = next->Right    = next->CircleAL = next->CircleAR
    //             = next->CircleBL = next->CircleBR = NULL;
    next->Right = NULL;

    return next;
}

void
EventCalendar::FreeNode(EventCalendar::CalendarNode* node)
{
    assert(node != NULL);
    node->Left = node->Parent
               = node->CircleAL
               = node->CircleAR
               = node->CircleBL
               = node->CircleBR = NULL;

    node->Time = -1.0;
    // node->ObjectA = node->ObjectB = -1;
    node->ObjectB = -1;
    node->Type = UpdateSystemEvent;

    node->Right = FreeList;
    FreeList = node;
}

/**
 * Schedule an event in the simulation.
 *
 * @param time The time the event will occur in simulation units.
 * @param type The type of the event.
 * @param idA The id of the first object in the event. This referes to the the
 *            first particle unless this is an UpdateSystemEvent.
 * @param idB The id of the second object in the event. This may correspond to
 *            other things (i.e. the cell wall being crossed) depending on the
 *            type of the event.
 */
void
EventCalendar::ScheduleEvent(const double time, const EventType type,
                             const int idA, const int idB)
{
    assert(time >= 0.0);
    GUARD_TYPE(type);

    EventCalendar::CalendarNode* newNode = NULL;

    // Use the ListAnchor for cell crossing events.
    const int assocCount = associatedObjects(type);

    newNode = this->AllocateNode();

    switch (associatedObjects(type))
    {
        case 2:
            ObjectAssociations[idB].second.push_back(newNode);
        case 1:
            ObjectAssociations[idA].first.push_back(newNode);
            break;
        default:
            break;
    }

    // If there is a collision, register the second particle with
    // appropriate event list.
    newNode->Time    = time;
    newNode->Type    = type;
    newNode->ObjectA = idA;
    newNode->ObjectB = idB;

    // Anchor the node based on the particles in the event.
    if (NULL == Root)
    {
        this->Root    = newNode;
        newNode->Left = newNode->Right = newNode->Parent = NULL;
    }
    else
    {
        this->Root->Insert(newNode);
    }
}

void
EventCalendar::NextEvent()
{
    assert(NULL != this->Root);

    EventCalendar::CalendarNode* curr = this->Root;

    this->HasStepped = true;

    // Go to the most recent event (i.e. the leftmost tree node).
    while (curr->Left != NULL)
    {
        curr = curr->Left;
    }

    // Copy event information to the current data buffer.
    this->CurrentEventTime    = curr->Time;
    this->CurrentEventType    = curr->Type;
    this->CurrentEventObjectA = curr->ObjectA;
    this->CurrentEventObjectB = curr->ObjectB;

    const int associated = associatedObjects(CurrentEventType);

    // Remove the most recent event from the tree and fixes the pointers.
    this->SafeDelete(curr);

    // Delete associated events as they are made invalid. Cell crossing
    // events do not invalidate other.
    if (associated >= 1 && invalidatesAssociated(this->CurrentEventType))
    {
        std::vector<CalendarNode*>& assocA = this->ObjectAssociations[CurrentObjectA].first;
        std::vector<CalendarNode*>& assocB = this->ObjectAssociations[CurrentObjectA].second;

        std::for_each(assoc.begin(), assoc.end(),
                [](CalendarNode* it)
                {
                    this->SafeDelete(it);
                    this->FreeNode(it);
                });

        this->ObjectAssociations[CurrentObjectA].second.clear();

        // Remove the second pair of particle lists corresponding to ObjectB
        if (associated == 2)
        {
            particle = &ListAnchor[CurrentEventObjectB];

            while (particle->CircleAL != particle)
            {
                EventCalendar::CalendarNode* next = particle->CircleAL;
                assert(next != curr);
                assert(next->Type != CellCrossingEvent);
                this->SafeDelete(next);
                this->FreeNode(next);
            }

            while (particle->CircleBL != particle)
            {
                EventCalendar::CalendarNode* next = particle->CircleBL;
                assert(next != curr);
                assert(next->Type != CellCrossingEvent);
                this->SafeDelete(next);
                this->FreeNode(next);
            }

            // The list should be empty afterwards.
            assert(particle->CircleAL == particle);
            assert(particle->CircleAR == particle);
            assert(particle->CircleBL == particle);
            assert(particle->CircleBR == particle);

            // Remove the cell crossing event from the calendar too. Only do this
            // if the particle isn't the current event that was just free, though.
            if (particle != curr && this->NodeInTree(particle))
            {
                this->SafeDelete(particle);
            }
        }
    }

    // Cell crossings are kept in the ListAnchor, so they do not need freed.
    if (CurrentEventType != CellCrossingEvent)
    {
        this->FreeNode(curr);
    }
}

bool
EventCalendar::HasMoreEvents()
{
    return this->Root != NULL;
}

double
EventCalendar::GetCurrentTime() const
{
    assert(HasStepped);
    return this->CurrentEventTime;
}

int
EventCalendar::GetCurrentObjectA() const
{
    assert(HasStepped);
    return this->CurrentEventObjectA;
}

int
EventCalendar::GetCurrentObjectB() const
{
    assert(HasStepped);
    return this->CurrentEventObjectB;
}

EventType
EventCalendar::GetCurrentEventType() const
{
    assert(HasStepped);
    return this->CurrentEventType;
}

EventCalendar::CalendarNode*
EventCalendar::SafeDelete(EventCalendar::CalendarNode* node)
{
    EventCalendar::CalendarNode* newRoot = node->Delete();

    if (this->Root == node)
    {
        this->Root = newRoot;
    }

    return newRoot;
}

