/**
 * Debugging output based on a compiler flag. If the macro <code>TRACE</code>
 * is not defined, then provided trace function will not do anything.
 * Otherwise it will be substitued for the printf function.
 */
#ifndef __TRACE_H__
#define __TRACE_H__

#include <execinfo.h>
#include <iostream>
#include <stdio.h>

#ifdef TRACE

#define trace(msg, ...) \
    printf(msg, ##__VA_ARGS__)

#define dout std::cout

#else

#define trace(msg, ...) \
    do {} while(false)

#define dout 0 && std::cout

#endif /* TRACE */

void print_trace();

#endif /* __TRACE_H__ */
