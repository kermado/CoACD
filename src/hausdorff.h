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

    inline double dist_point2point(const vec3d& pt, const vec3d& p)
    {
        const double dx = pt[0] - p[0];
        const double dy = pt[1] - p[1];
        const double dz = pt[2] - p[2];
        return sqrt(dx * dx + dy * dy + dz * dz);
    }

    double dist_point2segment(const vec3d& pt, const vec3d& s0, const vec3d& s1, bool flag = false)
    {
        assert(!flag);

        // first we build a 3d triangle the compute the height, pt as the top point
        vec3d BA, BC;
        BA[0] = pt[0] - s1[0];
        BA[1] = pt[1] - s1[1];
        BA[2] = pt[2] - s1[2];
        BC[0] = s0[0] - s1[0];
        BC[1] = s0[1] - s1[1];
        BC[2] = s0[2] - s1[2];

        // we calculate the projected vector along the segment
        const double lenBC = sqrt(BC[0] * BC[0] + BC[1] * BC[1] + BC[2] * BC[2]);
        const double proj_dist = (BA[0] * BC[0] + BA[1] * BC[1] + BA[2] * BC[2]) / lenBC;

        // we should make sure the projected point is within the segment, otherwise not consider it
        // if projected distance is negative or bigger than BC, it is out
        const double lenABsq = BA[0] * BA[0] + BA[1] * BA[1] + BA[2] * BA[2];
        if (proj_dist < 0.0 || proj_dist > lenBC) return INF;
        return sqrt(lenABsq - proj_dist * proj_dist);
    }

    double dist_point2triangle(const vec3d& pt, const vec3d& tri_pt0, const vec3d& tri_pt1, const vec3d& tri_pt2, bool flag = false)
    {
        // calculate the funciton of the plane, n = (a, b, c)
        const double _a = (tri_pt1[1] - tri_pt0[1]) * (tri_pt2[2] - tri_pt0[2]) - (tri_pt1[2] - tri_pt0[2]) * (tri_pt2[1] - tri_pt0[1]);
        const double _b = (tri_pt1[2] - tri_pt0[2]) * (tri_pt2[0] - tri_pt0[0]) - (tri_pt1[0] - tri_pt0[0]) * (tri_pt2[2] - tri_pt0[2]);
        const double _c = (tri_pt1[0] - tri_pt0[0]) * (tri_pt2[1] - tri_pt0[1]) - (tri_pt1[1] - tri_pt0[1]) * (tri_pt2[0] - tri_pt0[0]);
        const double _l = sqrt(_a * _a + _b * _b + _c * _c);
        const double a = _a / _l;
        const double b = _b / _l;
        const double c = _c / _l;
        const double d = -(a * tri_pt0[0] + b * tri_pt0[1] + c * tri_pt0[2]);

        // distance can be calculated directly using the function, then we get the projected point as well
        const double dist = fabs(a * pt[0] + b * pt[1] + c * pt[2] + d);
        const Plane p(a, b, c, d);
        const short side = p.Side(pt, 1e-8);

        vec3d proj_pt;
        if (side == 1)
        {
            proj_pt[0] = pt[0] - a * dist;
            proj_pt[1] = pt[1] - b * dist;
            proj_pt[2] = pt[2] - c * dist;
        }
        else if (side == -1)
        {
            proj_pt[0] = pt[0] + a * dist;
            proj_pt[1] = pt[1] + b * dist;
            proj_pt[2] = pt[2] + c * dist;
        }
        else
        {
            proj_pt = pt;
        }

        // judge if the projected point is within the triangle
        // we calculate the cross product of each edge and the line with the point, to see if the results are all along the normal direction
        const vec3d normal = CalFaceNormal(tri_pt0, tri_pt1, tri_pt2);
        vec3d AB, BC, CA, AP, BP, CP;
        AB[0] = tri_pt1[0] - tri_pt0[0];
        AB[1] = tri_pt1[1] - tri_pt0[1];
        AB[2] = tri_pt1[2] - tri_pt0[2];
        BC[0] = tri_pt2[0] - tri_pt1[0];
        BC[1] = tri_pt2[1] - tri_pt1[1];
        BC[2] = tri_pt2[2] - tri_pt1[2];
        CA[0] = tri_pt0[0] - tri_pt2[0];
        CA[1] = tri_pt0[1] - tri_pt2[1];
        CA[2] = tri_pt0[2] - tri_pt2[2];
        AP[0] = proj_pt[0] - tri_pt0[0];
        AP[1] = proj_pt[1] - tri_pt0[1];
        AP[2] = proj_pt[2] - tri_pt0[2];
        BP[0] = proj_pt[0] - tri_pt1[0];
        BP[1] = proj_pt[1] - tri_pt1[1];
        BP[2] = proj_pt[2] - tri_pt1[2];
        CP[0] = proj_pt[0] - tri_pt2[0];
        CP[1] = proj_pt[1] - tri_pt2[1];
        CP[2] = proj_pt[2] - tri_pt2[2];

        const vec3d AB_AP = CrossProduct(AB, AP);
        const vec3d BC_BP = CrossProduct(BC, BP);
        const vec3d CA_CP = CrossProduct(CA, CP);

        // if all the cross product are parallel with normal, then the projected point is within the triangle
        // else, we should calculate the shortest distance to three edges
        if (SameVectorDirection(AB_AP, normal) && SameVectorDirection(BC_BP, normal) && SameVectorDirection(CA_CP, normal))
        {
            return dist;
        }
        
        // if not within the triangle, we calculate the distance to 3 edges and 3 points and use the min
        const double dist_pt2AB = dist_point2segment(pt, tri_pt0, tri_pt1, flag);
        const double dist_pt2BC = dist_point2segment(pt, tri_pt1, tri_pt2, flag);
        const double dist_pt2CA = dist_point2segment(pt, tri_pt2, tri_pt0, flag);
        const double dist_pt2A = dist_point2point(pt, tri_pt0);
        const double dist_pt2B = dist_point2point(pt, tri_pt1);
        const double dist_pt2C = dist_point2point(pt, tri_pt2);

        assert(!flag);
        return min(min(min(dist_pt2AB, dist_pt2BC), dist_pt2CA), min(min(dist_pt2A, dist_pt2B), dist_pt2C));
    }

    PointCloud<double>& get_cloud()
    {
        static PointCloud<double> cloud;
        return cloud;
    }

    typedef KDTreeSingleIndexAdaptor<L2_Simple_Adaptor<double, PointCloud<double>>, PointCloud<double>, 3 /* dim */> my_kd_tree_t;

    double face_hausdorff_distance(const Model &meshA, const vector<vec3d> &XA, const vector<int> &idA, const Model &meshB, const vector<vec3d> &XB, const vector<int> &idB, bool flag = false)
    {
        const int nA = XA.size();
        const int nB = XB.size();
        double cmax = 0;

        PointCloud<double>& cloud = get_cloud();
        my_kd_tree_t index1(3 /*dim*/, cloud, KDTreeSingleIndexAdaptorParams(10 /* max leaf */));

        vec2PointCloud(cloud, XA);
        index1.buildIndex();

        size_t ret_index[10];
        double out_dist_sqr[10];

        for (int i = 0; i < nB; i++)
        {
            double query_pt[3] = {XB[i][0], XB[i][1], XB[i][2]};
            const size_t num_results = index1.knnSearch(&query_pt[0], 10, &ret_index[0], &out_dist_sqr[0]);

            double cmin = INF;
            for (int j = 0; j < (int)num_results; j++)
            {
                double distance;
                distance = dist_point2triangle(XB[i], meshA.points[meshA.triangles[idA[ret_index[j]]][0]], meshA.points[meshA.triangles[idA[ret_index[j]]][1]], meshA.points[meshA.triangles[idA[ret_index[j]]][2]]);
                if (distance < cmin)
                {
                    cmin = distance;
                    if (cmin < 1e-14)
                        break;
                }
            }
            if (cmin > 10)
                cmin = sqrt(out_dist_sqr[0]);
            if (cmin > cmax && INF > cmin)
                cmax = cmin;
        }

        vec2PointCloud(cloud, XB);
        my_kd_tree_t index2(3 /*dim*/, cloud, KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
        index2.buildIndex();

        for (int i = 0; i < nA; i++)
        {
            double query_pt[3] = {XA[i][0], XA[i][1], XA[i][2]};
            const size_t num_results = index2.knnSearch(&query_pt[0], 10, &ret_index[0], &out_dist_sqr[0]);

            double cmin = INF;
            for (int j = 0; j < (int)num_results; j++)
            {
                double distance;
                distance = dist_point2triangle(XA[i], meshB.points[meshB.triangles[idB[ret_index[j]]][0]], meshB.points[meshB.triangles[idB[ret_index[j]]][1]], meshB.points[meshB.triangles[idB[ret_index[j]]][2]]);
                if (distance < cmin)
                {
                    cmin = distance;
                    if (cmin < 1e-14)
                        break;
                }
            }
            if (cmin > 10)
                cmin = sqrt(out_dist_sqr[0]);
            if (cmin > cmax && INF > cmin)
                cmax = cmin;
        }

        return cmax;
    }
}