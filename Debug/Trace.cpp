
#include "Trace.h"

void print_trace()
{
    void* array[10];
    size_t size;

    // get void*'s for all entries on the stack
    size = backtrace(array, 10);

    // print out all the frames to stderr
    backtrace_symbols_fd(array, size, 2);
}
