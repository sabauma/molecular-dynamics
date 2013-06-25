
#ifndef __CONSTANTS_HELPER_H__
#define __CONSTANTS_HELPER_H__

#include <map>

namespace Constants
{

    /**
     * Function to invert an array mapping. Given an array and its length, this
     * function will produce a map that takes the values in the array and returns
     * the corresponding index. This assumes that each element of the array is
     * unique.
     *
     * @param source The array to convert.
     * @param len    The number of elements in the array.
     */
    template <typename source>
    inline std::map<source, int> invert_array_mapping(
            source const * const arr, const int len)
    {
        std::map<source, int> retval;

        for (int i = 0; i < len; ++i)
        {
            retval[arr[i]] = i;
        }

        return retval;
    }
}

#endif
