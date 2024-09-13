#pragma once

#include "vec3.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <string>
#include <vector>
#include <cstdlib>
#include <time.h>
#include <assert.h>
#include <math.h>

#include <list>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <array>
#include <limits>

using std::vector;
using std::array;
using std::min;
using std::max;
using std::string;
using std::pair;
using std::make_pair;
using std::ofstream;
using std::runtime_error;

namespace coacd
{
    constexpr double INF = std::numeric_limits<double>::max();

    struct HalfEdge
    {
        int head;
        std::list<HalfEdge>::iterator opposite;
    };

    class Edge
    {
    public:
        vec3d p0, p1;
        Edge() {}
        Edge(vec3d _p0, vec3d _p1)
        {
            p0 = _p0;
            p1 = _p1;
        }
    };

    class Plane
    {
    public:
        double a, b, c, d;
        vec3d p0, p1, p2; // three point form
        bool pFlag;       // whether three point form exists
        inline short CutSide(vec3d p0, vec3d p1, vec3d p2, Plane plane);
        inline short BoolSide(const vec3d& p) const;
        inline short Side(const vec3d& p, double eps = 1e-6) const
        {
            const double res = p[0] * a + p[1] * b + p[2] * c + d;
            if (res > eps)
                return 1;
            else if (res < -eps)
                return -1;
            return 0;
        }
        inline bool IntersectSegment(const vec3d& p1, const vec3d& p2, vec3d &pi, const double eps = 1e-6);
        inline Plane() : pFlag(false) {}
        inline Plane(double _a, double _b, double _c, double _d) : a(_a), b(_b), c(_c), d(_d), pFlag(false) {}
    };

    inline bool SamePointDetect(const vec3d& p0, const vec3d& p1, const double eps = 1e-5)
    {
        return fabs(p0[0] - p1[0]) < eps && fabs(p0[1] - p1[1]) < eps && fabs(p0[2] - p1[2]) < eps;
    }

    inline bool SameVectorDirection(const vec3d& v, const vec3d& w)
    {
        return (v[0] * w[0] + v[1] * w[1] + v[2] * w[2]) > 0.0;
    }

    inline vec3d CrossProduct(const vec3d& v, const vec3d& w)
    {
        return vec3d(v[1] * w[2] - v[2] * w[1], v[2] * w[0] - v[0] * w[2], v[0] * w[1] - v[1] * w[0]);
    }

    inline vec3d CalFaceNormal(const vec3d& p1, const vec3d& p2, const vec3d& p3)
    {
        vec3d v, w;
        v[0] = p2[0] - p1[0];
        v[1] = p2[1] - p1[1];
        v[2] = p2[2] - p1[2];
        w[0] = p3[0] - p1[0];
        w[1] = p3[1] - p1[1];
        w[2] = p3[2] - p1[2];

        const vec3d n = CrossProduct(v, w);
        const double l = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        return vec3d(n[0] / l, n[1] / l, n[2] / l);
    }

    inline double Area(const vec3d& p0, const vec3d& p1, const vec3d& p2)
    {
        const double x = p1[0] * p0[1] - p2[0] * p0[1] - p0[0] * p1[1] + p2[0] * p1[1] + p0[0] * p2[1] - p1[0] * p2[1];
        const double y = p1[0] * p0[2] - p2[0] * p0[2] - p0[0] * p1[2] + p2[0] * p1[2] + p0[0] * p2[2] - p1[0] * p2[2];
        const double z = p1[1] * p0[2] - p2[1] * p0[2] - p0[1] * p1[2] + p2[1] * p1[2] + p0[1] * p2[2] - p1[1] * p2[2];
        return 0.5 * sqrt(x * x + y * y + z * z);
    }

    inline double Volume(const vec3d& p1, const vec3d& p2, const vec3d& p3)
    {
        const double v321 = p3[0] * p2[1] * p1[2];
        const double v231 = p2[0] * p3[1] * p1[2];
        const double v312 = p3[0] * p1[1] * p2[2];
        const double v132 = p1[0] * p3[1] * p2[2];
        const double v213 = p2[0] * p1[1] * p3[2];
        const double v123 = p1[0] * p2[1] * p3[2];
        return (1.0 / 6.0) * (-v321 + v231 + v312 - v132 - v213 + v123);
    }

    void Diagonalize(const array<array<double, 3>, 3>& A, array<array<double, 3>, 3>& Q, array<array<double, 3>, 3>& D);

    inline double pt_norm(const vec3d& p)
    {
        return sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
    }

    inline short Plane::CutSide(vec3d p0, vec3d p1, vec3d p2, Plane plane)
    {
        vec3d normal = CalFaceNormal(p0, p1, p2);
        if (normal[0] * plane.a > 0 || normal[1] * plane.b > 0 || normal[2] * plane.c > 0)
            return -1;
        return 1;
    }

    inline short Plane::BoolSide(const vec3d& p) const
    {
        const double res = p[0] * a + p[1] * b + p[2] * c + d;
        return short((res > 0.0) ? 1 : -1);
    }

    inline bool Plane::IntersectSegment(const vec3d& p1, const vec3d& p2, vec3d &pi, const double eps)
    {
        pi[0] = (p1[0] * b * p2[1] + p1[0] * c * p2[2] + p1[0] * d - p2[0] * b * p1[1] - p2[0] * c * p1[2] - p2[0] * d) / (a * p2[0] - a * p1[0] + b * p2[1] - b * p1[1] + c * p2[2] - c * p1[2]);
        pi[1] = (a * p2[0] * p1[1] + c * p1[1] * p2[2] + p1[1] * d - a * p1[0] * p2[1] - c * p1[2] * p2[1] - p2[1] * d) / (a * p2[0] - a * p1[0] + b * p2[1] - b * p1[1] + c * p2[2] - c * p1[2]);
        pi[2] = (a * p2[0] * p1[2] + b * p2[1] * p1[2] + p1[2] * d - a * p1[0] * p2[2] - b * p1[1] * p2[2] - p2[2] * d) / (a * p2[0] - a * p1[0] + b * p2[1] - b * p1[1] + c * p2[2] - c * p1[2]);

        return min(p1[0] - eps, p2[0] - eps) <= pi[0] && pi[0] <= max(p1[0] + eps, p2[0] + eps) &&
               min(p1[1] - eps, p2[1] - eps) <= pi[1] && pi[1] <= max(p1[1] + eps, p2[1] + eps) &&
               min(p1[2] - eps, p2[2] - eps) <= pi[2] && pi[2] <= max(p1[2] + eps, p2[2] + eps);
    }

    template <typename T>
    struct PointCloud
    {
        struct Point
        {
            T x, y, z;
        };

        std::vector<Point> pts;

        // Must return the number of data points
        inline size_t kdtree_get_point_count() const { return pts.size(); }

        // Returns the dim'th component of the idx'th point in the class:
        // Since this is inlined and the "dim" argument is typically an immediate value, the
        //  "if/else's" are actually solved at compile time.
        inline T kdtree_get_pt(const size_t idx, const size_t dim) const
        {
            const T* ptr = (T* )(pts.data() + idx);
            return *(ptr + dim);
        }

        // Optional bounding-box computation: return false to default to a standard bbox computation loop.
        //   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
        //   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
        template <class BBOX>
        bool kdtree_get_bbox(BBOX & /* bb */) const { return false; }
    };

    template <typename T>
    void vec2PointCloud(PointCloud<T> &point, const vector<vec3d>& V)
    {
        const size_t size = V.size();
        point.pts.resize(size);
        for (size_t i = 0; i < size; i++)
        {
            const vec3d& v = V[i];
            point.pts[i].x = v[0];
            point.pts[i].y = v[1];
            point.pts[i].z = v[2];
        }
    }
}
