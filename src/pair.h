#pragma once

#include <cstdint>

namespace coacd
{
    union pair
    {
        int m[2];
        uint64_t i;

        inline pair() {}
        inline pair(const int first, const int second)
        : m{ first, second }
        {
            // Nothing to do.
        }
    };

    inline bool operator==(const pair& p1, const pair& p2)
    {
        return p1.i == p2.i;
    }

    inline bool operator!=(const pair& p1, const pair& p2)
    {
        return p1.i != p2.i;
    }
}