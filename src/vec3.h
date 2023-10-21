#pragma once

namespace coacd
{
    struct vec3d
    {
        double m[3];

        inline vec3d() {}
        inline vec3d(double x, double y, double z) : m{ x, y, z } {}

        inline double& operator[](const int index) { return m[index]; }
        inline double operator[](const int index) const { return m[index]; }
    };

    struct vec3i
    {
        int m[3];

        inline vec3i() {}
        inline vec3i(int x, int y, int z) : m{ x, y, z } {}

        inline int& operator[](const int index) { return m[index]; }
        inline int operator[](const int index) const { return m[index]; }
    };

    struct vec3s
    {
        short m[3];

        inline vec3s() {}
        inline vec3s(short x, short y, short z) : m{ x, y, z } {}

        inline short& operator[](const int index) { return m[index]; }
        inline short operator[](const int index) const { return m[index]; }
    };
}