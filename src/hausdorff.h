#pragma once
#include <typeinfo>
#include <algorithm>
#include <math.h>
#include <time.h>
#include <thread>
#include <assert.h>

#include <map>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>

#include "nanoflann.hpp"
#include "shape.h"
using namespace std;
using namespace nanoflann;

#define INF numeric_limits<double>::max()

namespace coacd
{

    inline double dist_point2point(const vec3d& p, const vec3d& q)
    {
        const double x = q[0] - p[0];
        const double y = q[1] - p[1];
        const double z = q[2] - p[2];
        return sqrt(x * x + y * y + z * z);
    }

    double dist_point2triangle(const vec3d& p, const vec3d& a, const vec3d& b, const vec3d& c, bool flag = false)
    {
        const vec3d ab(b[0] - a[0], b[1] - a[1], b[2] - a[2]);
        const vec3d ac(c[0] - a[0], c[1] - a[1], c[2] - a[2]);

        // check if p in vertex region outside a
        const vec3d ap(p[0] - a[0], p[1] - a[1], p[2] - a[2]);
        const double d1 = ab[0] * ap[0] + ab[1] * ap[1] + ab[2] * ap[2];
        const double d2 = ac[0] * ap[0] + ac[1] * ap[1] + ac[2] * ap[2];
        if (d1 <= 0.0 && d2 <= 0.0)
        {
            return sqrt(ap[0] * ap[0] + ap[1] * ap[1] + ap[2] * ap[2]);
        }

        // check if p in vertex region outside b
        const vec3d bp(p[0] - b[0], p[1] - b[1], p[2] - b[2]);
        const double d3 = ab[0] * bp[0] + ab[1] * bp[1] + ab[2] * bp[2];
        const double d4 = ac[0] * bp[0] + ac[1] * bp[1] + ac[2] * bp[2];
        if (d3 >= 0.0 && d4 <= d3)
        {
            return sqrt(bp[0] * bp[0] + bp[1] * bp[1] + bp[2] * bp[2]);
        }

        // check if p in vertex region outside c
        const vec3d cp(p[0] - c[0], p[1] - c[1], p[2] - c[2]);
        const double d5 = ab[0] * cp[0] + ab[1] * cp[1] + ab[2] * cp[2];
        const double d6 = ac[0] * cp[0] + ac[1] * cp[1] + ac[2] * cp[2];
        if (d6 >= 0.0 && d5 <= d6)
        {
            return sqrt(cp[0] * cp[0] + cp[1] * cp[1] + cp[2] * cp[2]);
        }

        // check if p in edge region of ab
        const double vc = d1 * d4 - d3 * d2;
        if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
        {
            const double v = d1 / (d1 - d3);
            const double x = a[0] + v * ab[0] - p[0];
            const double y = a[1] + v * ab[1] - p[1];
            const double z = a[2] + v * ab[2] - p[2];
            return sqrt(x * x + y * y + z * z);
        }

        // check if p in edge region of ac
        const double vb = d5 * d2 - d1 * d6;
        if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
        {
            const double w = d2 / (d2 - d6);
            const double x = a[0] + w * ac[0] - p[0];
            const double y = a[1] + w * ac[1] - p[1];
            const double z = a[2] + w * ac[2] - p[2];
            return sqrt(x * x + y * y + z * z);
        }

        // check if p in edge region of bc
        const double va = d3 * d6 - d5 * d4;
        const double d43 = d4 - d3;
        const double d56 = d5 - d6;
        if (va <= 0.0 && d43 >= 0.0 && d56 >= 0.0)
        {
            const double w = d43 / (d43 + d56);
            const double x = b[0] + w * (c[0] - b[0]) - p[0];
            const double y = b[1] + w * (c[1] - b[1]) - p[1];
            const double z = b[2] + w * (c[2] - b[2]) - p[2];
            return sqrt(x * x + y * y + z * z);
        }

        // p inside face region
        const double denom = 1.0 / (va + vb + vc);
        const double v = vb * denom;
        const double w = vc * denom;
        const double pq0 = a[0] + ab[0] * v + ac[0] * w - p[0];
        const double pq1 = a[1] + ab[1] * v + ac[1] * w - p[1];
        const double pq2 = a[2] + ab[2] * v + ac[2] * w - p[2];
        return sqrt(pq0 * pq0 + pq1 * pq1 + pq2 * pq2);
    }

    PointCloud<double>& get_cloud()
    {
        static PointCloud<double> cloud;
        return cloud;
    }

    typedef KDTreeSingleIndexAdaptor<L2_Simple_Adaptor<double, PointCloud<double>>, PointCloud<double>, 3 /* dim */> my_kd_tree_t;

    double face_hausdorff_distance(const Model &meshA, const vector<vec3d> &XA, const vector<int> &idA, const Model &meshB, const vector<vec3d> &XB, const vector<int> &idB, bool flag = false)
    {
        constexpr int n = 7;

        const int nA = XA.size();
        const int nB = XB.size();
        double cmax = 0;

        PointCloud<double>& cloud = get_cloud();
        my_kd_tree_t index1(3 /*dim*/, cloud, KDTreeSingleIndexAdaptorParams(10 /* max leaf */));

        vec2PointCloud(cloud, XA);
        index1.buildIndex();

        size_t ret_index[n];
        double out_dist_sqr[n];

        int seen[n];
        int seen_count = 0;

        for (int i = 0; i < nB; i++)
        {
            const int num_results = (int)index1.knnSearch(XB[i].m, n, ret_index, out_dist_sqr);

            double cmin = INF;
            for (int j = 0; j < num_results; j++)
            {
                const int idx = idA[ret_index[j]];
                for (int k = 0; k < seen_count; ++k)
                {
                    if (seen[k] == idx) { goto next_point; }
                }                

                {
                    const vec3i& triangle = meshA.triangles[idx];
                    const double distance = dist_point2triangle(XB[i], meshA.points[triangle[0]], meshA.points[triangle[1]], meshA.points[triangle[2]]);
                    if (distance < cmin)
                    {
                        cmin = distance;
                        if (cmin < 1e-14) { goto next_hull_point; }
                    }

                    seen[seen_count++] = idx;
                }

            next_point: (void)0;
            }

        next_hull_point:
            if (cmin > 10) { cmin = sqrt(out_dist_sqr[0]); }
            if (cmin > cmax && INF > cmin) { cmax = cmin; }
            seen_count = 0;
        }

        vec2PointCloud(cloud, XB);
        my_kd_tree_t index2(3 /*dim*/, cloud, KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
        index2.buildIndex();

        for (int i = 0; i < nA; i++)
        {
            const int num_results = (int)index2.knnSearch(XA[i].m, n, ret_index, out_dist_sqr);

            double cmin = INF;
            for (int j = 0; j < num_results; j++)
            {
                const int idx = idB[ret_index[j]];
                for (int k = 0; k < seen_count; ++k)
                {
                    if (seen[k] == idx) { goto next_point2; }
                }

                {
                    const vec3i& triangle = meshB.triangles[idx];
                    const double distance = dist_point2triangle(XA[i], meshB.points[triangle[0]], meshB.points[triangle[1]], meshB.points[triangle[2]]);
                    if (distance < cmin)
                    {
                        cmin = distance;
                        if (cmin < 1e-14) { goto next_hull_point2; }
                    }

                    seen[seen_count++] = idx;
                }

            next_point2: (void)0;
            }

        next_hull_point2:
            if (cmin > 10) { cmin = sqrt(out_dist_sqr[0]); }
            if (cmin > cmax && INF > cmin) { cmax = cmin; }
            seen_count = 0;
        }

        return cmax;
    }
}