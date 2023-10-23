#pragma once

#include "model_obj.h"
#include "cache.h"

namespace coacd
{
    bool Clip(const Model &mesh, Model &pos, Model &neg, Plane &plane, double &cut_area, bool foo = false);
    bool CreatePlaneRotationMatrix(const vector<vec3d>& border, vec3d& T, double R[3][3], Plane &plane);
    short Triangulation(vector<vec3d>& border, const vector<pair>& border_edges, vector<vec3i>& border_triangles, Plane& plane);
    void PrintEdgeSet(const vector<pair>& edges);

    inline void addPoint(hash_map<int, int>& vertex_map, vector<vec3d>& border, const vec3d& pt, int id, int &idx)
    {
        if (vertex_map.contains(id) == false)
        {
            int flag = -1;
            const int count = (int)border.size();
            for (int i = 0; i < count; i++)
            {
                const vec3d& br = border[i];
                if (fabs(br[0] - pt[0]) < 1e-4 && fabs(br[1] - pt[1]) < 1e-4 && fabs(br[2] - pt[2]) < 1e-4)
                {
                    flag = i;
                    break;
                }
            }
            if (flag == -1)
            {
                vertex_map[id] = idx;
                border.push_back(pt);
                idx++;
            }
            else
            {
                vertex_map[id] = flag;
            }
        }
    }

    inline void addEdgePoint(hash_map<pair, int, pair_hash>& edge_map, vector<vec3d>& border, const vec3d& pt, int id1, int id2, int &idx)
    {
        const pair edge1(id1, id2);
        const pair edge2(id2, id1);
        if (edge_map.contains(edge1) == false && edge_map.contains(edge2) == false)
        {
            int flag = -1;
            const int count = (int)border.size();
            for (int i = 0; i < count; i++)
            {
                const vec3d& br = border[i];
                if (fabs(br[0] - pt[0]) < 1e-4 && fabs(br[1] - pt[1]) < 1e-4 && fabs(br[2] - pt[2]) < 1e-4)
                {
                    flag = i;
                    break;
                }
            }
            if (flag == -1)
            {
                edge_map[edge1] = idx;
                edge_map[edge2] = idx;
                border.push_back(pt);
                idx++;
            }
            else
            {
                edge_map[edge1] = flag;
                edge_map[edge2] = flag;
            }
        }
    }

    inline bool FaceOverlap(const hash_set<int>& overlap_map, const vec3i& triangle)
    {
        return overlap_map.contains(triangle[0]) || overlap_map.contains(triangle[1]) || overlap_map.contains(triangle[2]);
    }
}
