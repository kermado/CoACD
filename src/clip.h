#pragma once

#include "model_obj.h"
#include <deque>

using std::deque;
using std::endl;

namespace coacd
{
    struct pair_hash
    {
        inline std::size_t operator()(const std::pair<int, int>& p) const
        {
            return p.first * 7 + p.second * 13;
        }
    };

    void SimpleCyclesFromEdges(const vector<pair<int, int>> edges, vector<vector<int>> &simple_cycles);
    void FindCycleDirection(vector<vec3d> border, vector<vector<int>> cycles, Plane plane, map<pair<int, int>, bool> &cycles_dir);
    void RemoveOutlierTriangles(const vector<vec3d>& border, const vector<vec3d>& overlap, const vector<pair<int, int>>& border_edges, int oriN,
                                const vector<vec3i>& border_triangles, vector<vec3i> &final_triangles);
    void RecordIntersectPoint(Model mesh, map<pair<int, int>, int> &edge_map, int i, int ep0, int ep1, int &idx, vector<vec3d> &border, vec3d point);
    bool Clip(const Model &mesh, Model &pos, Model &neg, Plane &plane, double &cut_area, bool foo = false);
    bool CreatePlaneRotationMatrix(vector<vec3d> &border, vec3d &T, double R[3][3], Plane &plane);
    short Triangulation(vector<vec3d> &border, const vector<pair<int, int>>& border_edges, vector<vec3i> &border_triangles, Plane &plane);
    void PrintEdgeSet(vector<pair<int, int>> edges);

    inline void addPoint(std::unordered_map<int, int> &vertex_map, vector<vec3d> &border, vec3d pt, int id, int &idx)
    {
        if (vertex_map.find(id) == vertex_map.end())
        {
            int flag = -1;
            for (int i = 0; i < (int)border.size(); i++)
            {
                if ((fabs(border[i][0] - pt[0])) < 1e-4 && (fabs(border[i][1] - pt[1])) < 1e-4 && (fabs(border[i][2] - pt[2])) < 1e-4)
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
                vertex_map[id] = flag;
        }
    }

    inline void addEdgePoint(std::unordered_map<pair<int, int>, int, pair_hash> &edge_map, vector<vec3d> &border, const vec3d& pt, int id1, int id2, int &idx)
    {
        pair<int, int> edge1 = make_pair(id1, id2);
        pair<int, int> edge2 = make_pair(id2, id1);
        if (edge_map.find(edge1) == edge_map.end() && edge_map.find(edge2) == edge_map.end())
        {
            int flag = -1;
            for (int i = 0; i < (int)border.size(); i++)
            {
                if ((fabs(border[i][0] - pt[0])) < 1e-4 && (fabs(border[i][1] - pt[1])) < 1e-4 && (fabs(border[i][2] - pt[2])) < 1e-4)
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

    inline bool FaceOverlap(const std::unordered_set<int>& overlap_map, const vec3i& triangle)
    {
        const auto end = overlap_map.cend();
        return overlap_map.find(triangle[0]) != end || overlap_map.find(triangle[1]) != end || overlap_map.find(triangle[2]) != end;
    }
}
