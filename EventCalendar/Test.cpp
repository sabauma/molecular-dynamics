
#include <stdio.h>

#include "EventType.h"
#include "EventCalendar.h"

const int EVENT_COUNT = 10000;

int main()
{
    EventCalendar ev(50000);
    ev.ScheduleEvent(2.0, CellCrossingEvent, 1, 2);
    // ev.ScheduleEvent(2.0, CellCrossingEvent, 1, 3);
    ev.ScheduleEvent(5.0, CellCrossingEvent, 2, 2);

    //for (int i = 0; i < EVENT_COUNT; ++i)
    //{
        //printf("Scheduling (%5.0f, %d, %d)\n", (double) i, i+1, i);
        //ev.ScheduleEvent((double) i, CellCrossingEvent, i, 0);
    //}

    for (int i = 0; i < EVENT_COUNT; ++i)
    {
        printf("Scheduling (%5.0f, %d, %d)\n", (double) i, 0, i+1);
        ev.ScheduleEvent((double) i, CollisionEvent, 0, i+1);
    }

    while (ev.HasMoreEvents())
    {
        ev.NextEvent();
        printf("Executing (%5.0f, %d, %d, %s)\n",
                ev.GetCurrentTime(),
                ev.GetCurrentObjectA(),
                ev.GetCurrentObjectB(),
                getEventName(ev.GetCurrentEventType()).c_str());
    }

    return 0;
}
