#include "clip.h"
#include "process.h"
#include "CDTUtils.h"
#include "CDT.h"

namespace coacd
{
    void Writepoints(vector<array<double, 2>> points, string filename)
    {
        std::ofstream os(filename);
        const int count = (int)points.size();
        for (int n = 0; n < count; n++)
        {
            os << "v " << points[n][0] << " " << points[n][1] << " 0"
                << "\n";
        }
        os.close();
    }

    void PrintEdgeSet(const vector<pair>& edges)
    {
        ofstream of("../edge.txt");
        const int count = (int)edges.size();
        for (int i = 0; i < count; i++)
        {
            of << i + 1 << ' ' << edges[i].m[0] << ' ' << edges[i].m[1] << endl;
        }
        of.close();
    }

    bool CreatePlaneRotationMatrix(const vector<vec3d>& border, vec3d& T, double R[3][3], Plane& plane)
    {
        int idx0 = 0;
        int idx1;
        int idx2;
        bool flag = 0;

        const int count = (int)border.size();
        const vec3d& br_idx0 = border[idx0];
        for (int i = 1; i < count; i++)
        {
            const vec3d& br_i = border[i];
            const double x = br_idx0[0] - br_i[0];
            const double y = br_idx0[1] - br_i[1];
            const double z = br_idx0[2] - br_i[2];
            double dist = sqrt(x * x + y * y + z * z);
            if (dist > 0.01)
            {
                flag = 1;
                idx1 = i;
                break;
            }
        }

        if (!flag) { return false; }
        flag = 0;

        vec3d p0 = border[idx0];
        vec3d p1 = border[idx1];

        for (int i = 2; i < count; i++)
        {
            if (i == idx1) { continue; }
            const vec3d& p2 = border[i];
            const vec3d AB = { p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2] };
            const vec3d BC = { p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2] };

            const double dot = AB[0] * BC[0] + AB[1] * BC[1] + AB[2] * BC[2];
            const double res = dot / (sqrt(AB[0] * AB[0] + AB[1] * AB[1] + AB[2] * AB[2]) * sqrt(BC[0] * BC[0] + BC[1] * BC[1] + BC[2] * BC[2]));
            if (abs(abs(res) - 1) > 1e-6 && abs(res) < INF) // AB not \\ BC, dot product != 1
            {
                flag = 1;
                idx2 = i;
                break;
            }
        }

        if (!flag) { return false; }

        vec3d p2 = border[idx2];
        const vec3d normal = CalFaceNormal(p0, p1, p2);
        if (normal[0] * plane.a > 0 || normal[1] * plane.b > 0 || normal[2] * plane.c > 0.0)
        {
            p0 = p2;
            p2 = border[idx0];
        }

        plane.pFlag = true;
        plane.p0 = p2;
        plane.p1 = p1;
        plane.p2 = p0;

        // translate to origin
        T = p0;

        // rotation matrix
        const double pd0 = p0[0] - p1[0];
        const double pd1 = p0[1] - p1[1];
        const double pd2 = p0[2] - p1[2];

        const double pd0_2 = pd0 * pd0;
        const double pd1_2 = pd1 * pd1;
        const double pd2_2 = pd2 * pd2;

        const double a = sqrt(pd0_2 + pd1_2 + pd2_2);
        R[0][0] = pd0 / a;
        R[0][1] = pd1 / a;
        R[0][2] = pd2 / a;

        const double px = p2[0] - p0[0];
        const double py = p2[1] - p0[1];
        const double pz = p2[2] - p0[2];

        const double rx = R[0][0];
        const double ry = R[0][1];
        const double rz = R[0][2];

        double t0 = pz * ry - py * rz;
        double t1 = px * rz - pz * rx;
        double t2 = py * rx - px * ry;

        const double b = sqrt(t0 * t0 + t1 * t1 + t2 * t2);
        R[2][0] = t0 / b;
        R[2][1] = t1 / b;
        R[2][2] = t2 / b;

        const double ru = R[2][0];
        const double rv = R[2][1];
        const double rw = R[2][2];

        t0 = rw * ry - rv * rz;
        t1 = ru * rz - rw * rx;
        t2 = rv * rx - ru * ry;

        const double c = sqrt(t0 * t0 + t1 * t1 + t2 * t2);
        R[1][0] = t0 / c;
        R[1][1] = t1 / c;
        R[1][2] = t2 / c;

        return true;
    }

    static int edge_start(const pair& p)
    {
        return p.m[0] - 1;
    }

    static int edge_end(const pair& p)
    {
        return p.m[1] - 1;
    }

    short Triangulation(vector<vec3d>& border, const vector<pair>& border_edges, vector<vec3i>& border_triangles, Plane& plane)
    {
        vec3d T;
        double R[3][3];
        if (!CreatePlaneRotationMatrix(border, T, R, plane)) return 1;

        const int count = (int)border.size();

        Cache& cache = Cache::get();
        cache.triangulate_prepare(count);
        vector<CDT::V2d<double>>& points = cache.points;

        const double R00 = R[0][0];
        const double R01 = R[0][1];
        const double R02 = R[0][2];
        const double R10 = R[1][0];
        const double R11 = R[1][1];
        const double R12 = R[1][2];

        const double T0 = T[0];
        const double T1 = T[1];
        const double T2 = T[2];

        double x_min = INF, x_max = -INF, y_min = INF, y_max = -INF;
        for (int i = 0; i < count; i++)
        {
            const vec3d& v = border[i];
            const double x = v[0] - T0;
            const double y = v[1] - T1;
            const double z = v[2] - T2;

            const double px = R00 * x + R01 * y + R02 * z;
            const double py = R10 * x + R11 * y + R12 * z;

            CDT::V2d<double>& p = points[i];
            p.x = px;
            p.y = py;

            x_min = min(x_min, px);
            x_max = max(x_max, px);
            y_min = min(y_min, py);
            y_max = max(y_max, py);
        }
        
        CDT::Triangulation<double> cdt;
        cdt.insertVertices(points);
        cdt.insertEdges(border_edges.begin(), border_edges.end(), edge_start, edge_end);
        cdt.eraseSuperTriangle();

        const CDT::TriangleVec& triangles = cdt.triangles;
        const int triangle_count = (int)triangles.size();
        border_triangles.reserve(triangle_count);
        for (int i = 0; i < triangle_count; i++)
        {
            const CDT::Triangle& triangle = triangles[i];
            const CDT::VerticesArr3& vertices = triangle.vertices;
            border_triangles.emplace_back((int)vertices[0] + 1, (int)vertices[1] + 1, (int)vertices[2] + 1);
        }

        const vector<CDT::V2d<double>>& out_vertices = cdt.vertices;
        const int point_count = (int)points.size();
        border.reserve(point_count);
        for (int i = count; i < point_count; i++)
        {
            const CDT::V2d<double>& vertex = out_vertices[i];
            const double x = R00 * vertex.x + R10 * vertex.y + T0;
            const double y = R01 * vertex.x + R11 * vertex.y + T1;
            const double z = R02 * vertex.x + R12 * vertex.y + T2;
            border.emplace_back(x, y, z);
        }

        return 0;
    }

    void RemoveOutlierTriangles(const vector<vec3d>& border, const vector<vec3d>& overlap, const vector<pair>& border_edges, const vector<vec3i>& border_triangles, int oriN,
        hash_map<int, int>& vertex_map, vector<vec3d>& final_border, vector<vec3i>& final_triangles)
    {
        const int v_length = (int)border.size();
        const int f_length = (int)border_triangles.size();
        Cache& cache = Cache::get();
        cache.outliers_prepare(v_length, f_length);
        std::vector<pair>& BFS_edges = cache.BFS_edges;
        hash_map<pair, pair, pair_hash>& edge_map = cache.tri_edge_map;
        hash_set<pair, pair_hash>& border_map = cache.tri_border_map;
        hash_set<pair, pair_hash>& same_edge_map = cache.same_edge_map;
        hash_set<int>& overlap_map = cache.overlap_map;
        std::vector<bool>& add_vertex = cache.add_vertex;
        std::vector<bool>& remove_map = cache.remove_map;

        BFS_edges.insert(BFS_edges.begin(), border_edges.begin(), border_edges.end());

        for (int i = 0; i < v_length; ++i) { add_vertex[i] = false; }
        for (int i = 0; i < f_length; ++i) { remove_map[i] = false; }

        const int overlap_count = (int)overlap.size();
        for (int i = 0; i < overlap_count; i++)
        {
            const vec3d& p = overlap[i];
            for (int j = 0; j < v_length; j++)
            {
                if (SamePointDetect(p, border[j]))
                {
                    overlap_map.insert(j + 1);
                }
            }
        }

        const int border_edges_count = (int)border_edges.size();
        for (int i = 0; i < border_edges_count; i++)
        {
            same_edge_map.insert(border_edges[i]);
        }

        for (int i = 0; i < border_edges_count; i++)
        {
            const pair& be = border_edges[i];
            const int v0 = be.m[0];
            const int v1 = be.m[1];
            const pair edge(v1, v0);
            if (same_edge_map.contains(edge) == false)
            {
                border_map.insert(pair(v0, v1));
                border_map.insert(edge);
            }
        }

        const int border_triangles_count = (int)border_triangles.size();
        for (int i = 0; i < border_triangles_count; i++)
        {
            const int v0 = border_triangles[i][0];
            const int v1 = border_triangles[i][1];
            const int v2 = border_triangles[i][2];

            if (v0 < 1 || v0 > v_length || v1 < 1 || v1 > v_length || v2 < 1 || v2 > v_length)
            {
                continue; // ignore points added by triangle
            }

            const pair edge01(v0, v1);
            const pair edge10(v1, v0);
            const pair edge12(v1, v2);
            const pair edge21(v2, v1);
            const pair edge20(v2, v0);
            const pair edge02(v0, v2);

            if (!same_edge_map.contains(edge10) || !same_edge_map.contains(edge01))
            {
                if (edge_map.contains(edge10) == false)
                {
                    edge_map[edge01] = pair(i, -1);
                }
                else
                {
                    edge_map[edge10] = pair(edge_map[edge10].m[0], i);
                }
            }

            if (!same_edge_map.contains(edge12) || !same_edge_map.contains(edge21))
            {
                if (edge_map.contains(edge21) == false)
                {
                    edge_map[edge12] = pair(i, -1);
                }
                else
                {
                    edge_map[edge21] = pair(edge_map[edge21].m[0], i);
                }
            }

            if (!same_edge_map.contains(edge02) || !same_edge_map.contains(edge20))
            {
                if (edge_map.contains(edge02) == false)
                {
                    edge_map[edge20] = pair(i, -1);
                }
                else
                {
                    edge_map[edge02] = pair(edge_map[edge02].m[0], i);
                }
            }
        }

        int front = 0;
        int i = 0;
        while (front < BFS_edges.size())
        {
            const pair& item = BFS_edges[front++];
            const int v0 = item.m[0];
            const int v1 = item.m[1];

            int idx;
            const pair edge01(v0, v1);
            const pair edge10(v1, v0);
            if (i < border_edges_count && edge_map.contains(edge10))
            {
                const pair& p = edge_map[edge10];
                idx = p.m[1];
                if (idx != -1) { remove_map[idx] = true; }
                idx = p.m[0];
                if (idx != -1)
                {
                    const vec3i& triangle = border_triangles[idx];
                    if (!remove_map[idx] && !FaceOverlap(overlap_map, triangle))
                    {
                        remove_map[idx] = true;
                        final_triangles.push_back(triangle);

                        const int p0 = triangle[0];
                        const int p1 = triangle[1];
                        const int p2 = triangle[2];

                        add_vertex[p0 - 1] = true;
                        add_vertex[p1 - 1] = true;
                        add_vertex[p2 - 1] = true;

                        if (p2 != v0 && p2 != v1)
                        {
                            const pair pt12(p1, p2);
                            const pair pt20(p2, p0);
                            if (border_map.contains(pt12) == false)
                            {
                                BFS_edges.push_back(pt12);
                            }
                            if (border_map.contains(pt20) == false)
                            {
                                BFS_edges.push_back(pt20);
                            }
                        }
                        else if (p1 != v0 && p1 != v1)
                        {
                            const pair pt12(p1, p2);
                            const pair pt01(p0, p1);
                            if (border_map.contains(pt12) == false)
                            {
                                BFS_edges.push_back(pt12);
                            }
                            if (border_map.contains(pt01) == false)
                            {
                                BFS_edges.push_back(pt01);
                            }
                        }
                        else if (p0 != v0 && p0 != v1)
                        {
                            const pair pt01(p0, p1);
                            const pair pt20(p2, p0);
                            if (border_map.contains(pt01) == false)
                            {
                                BFS_edges.push_back(pt01);
                            }
                            if (border_map.contains(pt20) == false)
                            {
                                BFS_edges.push_back(pt20);
                            }
                        }
                    }
                }
            }
            else if (i < border_edges_count && edge_map.contains(edge01))
            {
                const pair& p = edge_map[edge01];
                idx = p.m[0];
                if (idx != -1) { remove_map[idx] = true; }
                idx = p.m[1];
                if (idx != -1)
                {
                    const vec3i& triangle = border_triangles[idx];
                    if (!remove_map[idx] && !FaceOverlap(overlap_map, triangle))
                    {
                        remove_map[idx] = true;
                        final_triangles.push_back(triangle);

                        const int p0 = triangle[0];
                        const int p1 = triangle[1];
                        const int p2 = triangle[2];

                        add_vertex[p0 - 1] = true;
                        add_vertex[p1 - 1] = true;
                        add_vertex[p2 - 1] = true;

                        if (p2 != v0 && p2 != v1)
                        {
                            const pair pt21(p2, p1);
                            const pair pt02(p0, p2);
                            if (border_map.contains(pt21) == false)
                            {
                                BFS_edges.push_back(pt21);
                            }
                            if (border_map.contains(pt02) == false)
                            {
                                BFS_edges.push_back(pt02);
                            }
                        }
                        else if (p1 != v0 && p1 != v1)
                        {
                            const pair pt21(p2, p1);
                            const pair pt10(p1, p0);
                            if (border_map.contains(pt21) == false)
                            {
                                BFS_edges.push_back(pt21);
                            }
                            if (border_map.contains(pt10) == false)
                            {
                                BFS_edges.push_back(pt10);
                            }
                        }
                        else if (p0 != v0 && p0 != v1)
                        {
                            const pair pt10(p1, p0);
                            const pair pt02(p0, p2);
                            if (border_map.contains(pt10) == false)
                            {
                                BFS_edges.push_back(pt10);
                            }
                            if (border_map.contains(pt02) == false)
                            {
                                BFS_edges.push_back(pt02);
                            }
                        }
                    }
                }
            }
            else if (i >= border_edges_count && (edge_map.contains(edge01) || edge_map.contains(edge10)))
            {
                for (int j = 0; j < 2; j++)
                {
                    if (j == 0)
                    {
                        if (edge_map.contains(edge01))
                        {
                            idx = edge_map[edge01].m[0];
                        }
                        else
                        {
                            idx = edge_map[edge10].m[0];
                        }
                    }
                    else if (edge_map.contains(edge01))
                    {
                        idx = edge_map[edge01].m[1];
                    }
                    else
                    {
                        idx = edge_map[edge10].m[1];
                    }
                    if (idx != -1 && !remove_map[idx])
                    {
                        remove_map[idx] = true;
                        const vec3i& triangle = border_triangles[idx];
                        final_triangles.push_back(triangle);

                        const int p0 = triangle[0];
                        const int p1 = triangle[1];
                        const int p2 = triangle[2];

                        add_vertex[p0 - 1] = true;
                        add_vertex[p1 - 1] = true;
                        add_vertex[p2 - 1] = true;

                        if (p2 != v0 && p2 != v1)
                        {
                            BFS_edges.emplace_back(p1, p2);
                            BFS_edges.emplace_back(p2, p0);
                        }
                        else if (p1 != v0 && p1 != v1)
                        {
                            BFS_edges.emplace_back(p1, p2);
                            BFS_edges.emplace_back(p0, p1);
                        }
                        else if (p0 != v0 && p0 != v1)
                        {
                            BFS_edges.emplace_back(p0, p1);
                            BFS_edges.emplace_back(p2, p0);
                        }
                    }
                }
            }
            i++;
        }

        int index = 0;
        for (int i = 0; i < v_length; i++)
        {
            if (i < oriN || add_vertex[i])
            {
                final_border.push_back(border[i]);
                vertex_map[i + 1] = ++index;
            }
        }
    }

    bool Clip(const Model& mesh, Model& pos, Model& neg, Plane& plane, double& cut_area, bool foo)
    {
        const int N = mesh.points.size();

        Cache& cache = Cache::get();
        cache.clip_prepare(N);
        vector<vec3d>& border = cache.border;
        vector<vec3d>& overlap = cache.overlap;
        vector<vec3i>& border_triangles = cache.border_triangles;
        vector<vec3i>& final_triangles = cache.final_triangles;
        vector<pair>& border_edges = cache.border_edges;
        hash_map<int, int>& border_map = cache.border_map;
        vector<vec3d>& final_border = cache.final_border;
        vector<bool>& pos_map = cache.pos_map;
        vector<bool>& neg_map = cache.neg_map;
        vector<int>& pos_proj = cache.pos_proj;
        vector<int>& neg_proj = cache.neg_proj;
        hash_map<pair, int, pair_hash>& edge_map = cache.edge_map;
        hash_map<int, int>& vertex_map = cache.vertex_map;

        for (int i = 0; i < N; ++i) { pos_map[i] = false; }
        for (int i = 0; i < N; ++i) { neg_map[i] = false; }
        std::memset(pos_proj.data(), 0, sizeof(int) * N);
        std::memset(neg_proj.data(), 0, sizeof(int) * N);

        int idx = 0;

        const int triangle_count = (int)mesh.triangles.size();
        for (int i = 0; i < triangle_count; i++)
        {
            const vec3i& triangle = mesh.triangles[i];
            int id0 = triangle[0];
            int id1 = triangle[1];
            int id2 = triangle[2];

            vec3d p0, p1, p2;
            p0 = mesh.points[id0];
            p1 = mesh.points[id1];
            p2 = mesh.points[id2];

            short s0 = plane.Side(p0);
            short s1 = plane.Side(p1);
            short s2 = plane.Side(p2);
            short sum = s0 + s1 + s2;

            if (s0 == 0 && s1 == 0 && s2 == 0)
            {
                s0 = s1 = s2 = plane.CutSide(p0, p1, p2, plane);
                sum = s0 + s1 + s2;
                overlap.push_back(p0);
                overlap.push_back(p1);
                overlap.push_back(p2);
            }

            if (sum == 3 || sum == 2 || (sum == 1 && ((s0 == 1 && s1 == 0 && s2 == 0) || (s0 == 0 && s1 == 1 && s2 == 0) || (s0 == 0 && s1 == 0 && s2 == 1)))) // pos side
            {
                pos_map[id0] = true;
                pos_map[id1] = true;
                pos_map[id2] = true;
                pos.triangles.push_back(triangle);
                // the plane cross the triangle edge
                if (sum == 1)
                {
                    if (s0 == 1 && s1 == 0 && s2 == 0)
                    {
                        addPoint(vertex_map, border, p1, id1, idx);
                        addPoint(vertex_map, border, p2, id2, idx);
                        const int v1 = vertex_map[id1];
                        const int v2 = vertex_map[id2];
                        if (v1 != v2)
                        {
                            border_edges.emplace_back(v1 + 1, v2 + 1);
                        }
                    }
                    else if (s0 == 0 && s1 == 1 && s2 == 0)
                    {
                        addPoint(vertex_map, border, p2, id2, idx);
                        addPoint(vertex_map, border, p0, id0, idx);
                        const int v2 = vertex_map[id2];
                        const int v0 = vertex_map[id0];
                        if (v2 != v0)
                        {
                            border_edges.emplace_back(v2 + 1, v0 + 1);
                        }
                    }
                    else if (s0 == 0 && s1 == 0 && s2 == 1)
                    {
                        addPoint(vertex_map, border, p0, id0, idx);
                        addPoint(vertex_map, border, p1, id1, idx);
                        const int v0 = vertex_map[id0];
                        const int v1 = vertex_map[id1];
                        if (v0 != v1)
                        {
                            border_edges.emplace_back(v0 + 1, v1 + 1);
                        }
                    }
                }
            }
            else if (sum == -3 || sum == -2 || (sum == -1 && ((s0 == -1 && s1 == 0 && s2 == 0) || (s0 == 0 && s1 == -1 && s2 == 0) || (s0 == 0 && s1 == 0 && s2 == -1)))) // neg side
            {
                neg_map[id0] = true;
                neg_map[id1] = true;
                neg_map[id2] = true;
                neg.triangles.push_back(triangle);
                // the plane cross the triangle edge
                if (sum == -1)
                {
                    if (s0 == -1 && s1 == 0 && s2 == 0)
                    {
                        addPoint(vertex_map, border, p2, id2, idx);
                        addPoint(vertex_map, border, p1, id1, idx);
                        const int v2 = vertex_map[id2];
                        const int v1 = vertex_map[id1];
                        if (v2 != v1)
                        {
                            border_edges.emplace_back(v2 + 1, v1 + 1);
                        }
                    }
                    else if (s0 == 0 && s1 == -1 && s2 == 0)
                    {
                        addPoint(vertex_map, border, p0, id0, idx);
                        addPoint(vertex_map, border, p2, id2, idx);
                        const int v0 = vertex_map[id0];
                        const int v2 = vertex_map[id2];
                        if (v0 != v2)
                        {
                            border_edges.emplace_back(v0 + 1, v2 + 1);
                        }
                    }
                    else if (s0 == 0 && s1 == 0 && s2 == -1)
                    {
                        addPoint(vertex_map, border, p1, id1, idx);
                        addPoint(vertex_map, border, p0, id0, idx);
                        const int v1 = vertex_map[id1];
                        const int v0 = vertex_map[id0];
                        if (v1 != v0)
                        {
                            border_edges.emplace_back(v1 + 1, v0 + 1);
                        }
                    }
                }
            }
            else // different side
            {
                vec3d pi0, pi1, pi2;
                bool f0 = plane.IntersectSegment(p0, p1, pi0);
                bool f1 = plane.IntersectSegment(p1, p2, pi1);
                bool f2 = plane.IntersectSegment(p2, p0, pi2);

                if (f0 && f1 && !f2)
                {
                    // record the points
                    addEdgePoint(edge_map, border, pi0, id0, id1, idx); // f0
                    addEdgePoint(edge_map, border, pi1, id1, id2, idx); // f1

                    // record the edges
                    const int f0_idx = edge_map[pair(id0, id1)];
                    const int f1_idx = edge_map[pair(id1, id2)];
                    const int fd0 = -1 * f0_idx - 1;
                    const int fd1 = -1 * f1_idx - 1;
                    if (s1 == 1)
                    {
                        if (f1_idx != f0_idx)
                        {
                            border_edges.emplace_back(f1_idx + 1, f0_idx + 1); // border
                            pos_map[id1] = true;
                            neg_map[id0] = true;
                            neg_map[id2] = true;
                            pos.triangles.emplace_back(id1, fd1, fd0); // make sure it is not zero
                            neg.triangles.emplace_back(id0, fd0, fd1);
                            neg.triangles.emplace_back(fd1, id2, id0);
                        }
                        else
                        {
                            neg_map[id0] = true;
                            neg_map[id2] = true;
                            neg.triangles.emplace_back(fd1, id2, id0);
                        }
                    }
                    else
                    {
                        if (f0_idx != f1_idx)
                        {
                            border_edges.emplace_back(f0_idx + 1, f1_idx + 1); // border
                            neg_map[id1] = true;
                            pos_map[id0] = true;
                            pos_map[id2] = true;
                            neg.triangles.emplace_back(id1, fd1, fd0);
                            pos.triangles.emplace_back(id0, fd0, fd1);
                            pos.triangles.emplace_back(fd1, id2, id0);
                        }
                        else
                        {
                            pos_map[id0] = true;
                            pos_map[id2] = true;
                            pos.triangles.emplace_back(fd1, id2, id0);
                        }
                    }
                }
                else if (f1 && f2 && !f0)
                {
                    // record the points
                    addEdgePoint(edge_map, border, pi1, id1, id2, idx); // f1
                    addEdgePoint(edge_map, border, pi2, id2, id0, idx); // f2

                    // record the edges
                    const int f1_idx = edge_map[pair(id1, id2)];
                    const int f2_idx = edge_map[pair(id2, id0)];
                    const int fd1 = -1 * f1_idx - 1;
                    const int fd2 = -1 * f2_idx - 1;
                    if (s2 == 1)
                    {
                        if (f2_idx != f1_idx)
                        {
                            border_edges.emplace_back(f2_idx + 1, f1_idx + 1);
                            pos_map[id2] = true;
                            neg_map[id0] = true;
                            neg_map[id1] = true;
                            pos.triangles.emplace_back(id2, fd2, fd1);
                            neg.triangles.emplace_back(id0, fd1, fd2);
                            neg.triangles.emplace_back(fd1, id0, id1);
                        }
                        else
                        {
                            neg_map[id0] = true;
                            neg_map[id1] = true;
                            neg.triangles.emplace_back(fd1, id0, id1);
                        }
                    }
                    else
                    {
                        if (f1_idx != f2_idx)
                        {
                            border_edges.emplace_back(f1_idx + 1, f2_idx + 1);
                            neg_map[id2] = true;
                            pos_map[id0] = true;
                            pos_map[id1] = true;
                            neg.triangles.emplace_back(id2, fd2, fd1);
                            pos.triangles.emplace_back(id0, fd1, fd2);
                            pos.triangles.emplace_back(fd1, id0, id1);
                        }
                        else
                        {
                            pos_map[id0] = true;
                            pos_map[id1] = true;
                            pos.triangles.emplace_back(fd1, id0, id1);
                        }
                    }
                }
                else if (f2 && f0 && !f1)
                {
                    // record the points
                    addEdgePoint(edge_map, border, pi2, id2, id0, idx); // f2
                    addEdgePoint(edge_map, border, pi0, id0, id1, idx); // f0

                    const int f0_idx = edge_map[pair(id0, id1)];
                    const int f2_idx = edge_map[pair(id2, id0)];
                    const int fd0 = -1 * f0_idx - 1;
                    const int fd2 = -1 * f2_idx - 1;
                    if (s0 == 1)
                    {
                        if (f0_idx != f2_idx)
                        {
                            border_edges.emplace_back(f0_idx + 1, f2_idx + 1);
                            pos_map[id0] = true;
                            neg_map[id1] = true;
                            neg_map[id2] = true;
                            pos.triangles.emplace_back(id0, fd0, fd2);
                            neg.triangles.emplace_back(id1, fd2, fd0);
                            neg.triangles.emplace_back(fd2, id1, id2);
                        }
                        else
                        {
                            neg_map[id1] = true;
                            neg_map[id2] = true;
                            neg.triangles.emplace_back(fd2, id1, id2);
                        }
                    }
                    else
                    {
                        if (f2_idx != f0_idx)
                        {
                            border_edges.emplace_back(f2_idx + 1, f0_idx + 1);
                            neg_map[id0] = true;
                            pos_map[id1] = true;
                            pos_map[id2] = true;
                            neg.triangles.emplace_back(id0, fd0, fd2);
                            pos.triangles.emplace_back(id1, fd2, fd0);
                            pos.triangles.emplace_back(fd2, id1, id2);
                        }
                        else
                        {
                            pos_map[id1] = true;
                            pos_map[id2] = true;
                            pos.triangles.emplace_back(fd2, id1, id2);
                        }
                    }
                }
                else if (f0 && f1 && f2)
                {
                    if (s0 == 0 || (s0 != 0 && s1 != 0 && s2 != 0 && SamePointDetect(pi0, pi2))) // intersect at p0
                    {
                        // f2 = f0 = p0
                        addPoint(vertex_map, border, p0, id0, idx);
                        const int v = vertex_map[id0];
                        edge_map[pair(id0, id1)] = v;
                        edge_map[pair(id1, id0)] = v;
                        edge_map[pair(id2, id0)] = v;
                        edge_map[pair(id0, id2)] = v;

                        // record the points
                        addEdgePoint(edge_map, border, pi1, id1, id2, idx); // f1

                        const int f1_idx = edge_map[pair(id1, id2)];
                        const int f0_idx = v;
                        const int fd1 = -1 * f1_idx - 1;
                        const int fd0 = -1 * f0_idx - 1;
                        if (s1 == 1)
                        {
                            if (f1_idx != f0_idx)
                            {
                                border_edges.emplace_back(f1_idx + 1, f0_idx + 1);
                                pos_map[id1] = true;
                                neg_map[id2] = true;
                                pos.triangles.emplace_back(id1, fd1, fd0);
                                neg.triangles.emplace_back(id2, fd0, fd1);
                            }
                        }
                        else
                        {
                            if (f0_idx != f1_idx)
                            {
                                border_edges.emplace_back(f0_idx + 1, f1_idx + 1);
                                neg_map[id1] = true;
                                pos_map[id2] = true;
                                neg.triangles.emplace_back(id1, fd1, fd0);
                                pos.triangles.emplace_back(id2, fd0, fd1);
                            }
                        }
                    }
                    else if (s1 == 0 || (s0 != 0 && s1 != 0 && s2 != 0 && SamePointDetect(pi0, pi1))) // intersect at p1
                    {
                        // f0 = f1 = p1
                        addPoint(vertex_map, border, p1, id1, idx);
                        const int v = vertex_map[id1];
                        edge_map[pair(id0, id1)] = v;
                        edge_map[pair(id1, id0)] = v;
                        edge_map[pair(id1, id2)] = v;
                        edge_map[pair(id2, id1)] = v;

                        // record the points
                        addEdgePoint(edge_map, border, pi2, id2, id0, idx); // f2

                        const int f1_idx = v;
                        const int f2_idx = edge_map[pair(id2, id0)];
                        const int fd1 = -1 * f1_idx - 1;
                        const int fd2 = -1 * f2_idx - 1;
                        if (s0 == 1)
                        {
                            if (f1_idx != f2_idx)
                            {
                                border_edges.emplace_back(f1_idx + 1, f2_idx + 1);
                                pos_map[id0] = true;
                                neg_map[id2] = true;
                                pos.triangles.emplace_back(id0, fd1, fd2);
                                neg.triangles.emplace_back(id2, fd2, fd1);
                            }
                        }
                        else
                        {
                            if (f2_idx != f1_idx)
                            {
                                border_edges.emplace_back(f2_idx + 1, f1_idx + 1);
                                neg_map[id0] = true;
                                pos_map[id2] = true;
                                neg.triangles.emplace_back(id0, fd1, fd2);
                                pos.triangles.emplace_back(id2, fd2, fd1);
                            }
                        }
                    }
                    else if (s2 == 0 || (s0 != 0 && s1 != 0 && s2 != 0 && SamePointDetect(pi1, pi2))) // intersect at p2
                    {
                        // f1 = f2 = p2
                        addPoint(vertex_map, border, p2, id2, idx);
                        const int v = vertex_map[id2];
                        edge_map[pair(id1, id2)] = v;
                        edge_map[pair(id2, id1)] = v;
                        edge_map[pair(id2, id0)] = v;
                        edge_map[pair(id0, id2)] = v;

                        // record the points
                        addEdgePoint(edge_map, border, pi0, id0, id1, idx); // f0

                        const int f0_idx = edge_map[pair(id0, id1)];
                        const int f1_idx = v;
                        const int fd0 = -1 * f0_idx - 1;
                        const int fd1 = -1 * f1_idx - 1;
                        if (s0 == 1)
                        {
                            if (f0_idx != f1_idx)
                            {
                                border_edges.emplace_back(f0_idx + 1, f1_idx + 1);
                                pos_map[id0] = true;
                                neg_map[id1] = true;
                                pos.triangles.emplace_back(id0, fd0, fd1);
                                neg.triangles.emplace_back(id1, fd1, fd0);
                            }
                        }
                        else
                        {
                            if (f1_idx != f0_idx)
                            {
                                border_edges.emplace_back(f1_idx + 1, f0_idx + 1);
                                neg_map[id0] = true;
                                pos_map[id1] = true;
                                neg.triangles.emplace_back(id0, fd0, fd1);
                                pos.triangles.emplace_back(id1, fd1, fd0);
                            }
                        }
                    }
                    else
                    {
                        assert(false && "Intersection error"); // Intersection error. Please report this error to sarahwei0210@gmail.com with your input OBJ and log file.
                        return false;
                    }
                }
            }
        }

        if (border.size() > 2)
        {
            const short flag = Triangulation(border, border_edges, border_triangles, plane);
            if (flag == 0)
            {
                const int oriN = (int)border.size();
                RemoveOutlierTriangles(border, overlap, border_edges, border_triangles, oriN, border_map, final_border, final_triangles);
            }
            else if (flag == 1)
            {
                final_border = border; // remember to fill final_border with border!
            }
            else
            {
                return false; // clip failed
            }
            cut_area = 0;
        }
        else
        {
            final_border = border; // remember to fill final_border with border!
            cut_area = -10;
        }

        // original points in two parts
        double pos_x_min = INF, pos_x_max = -INF, pos_y_min = INF, pos_y_max = -INF, pos_z_min = INF, pos_z_max = -INF;
        double neg_x_min = INF, neg_x_max = -INF, neg_y_min = INF, neg_y_max = -INF, neg_z_min = INF, neg_z_max = -INF;

        int pos_idx = 0;
        int neg_idx = 0;
        for (int i = 0; i < N; i++)
        {
            const vec3d& pi = mesh.points[i];
            const double p0 = pi[0];
            const double p1 = pi[1];
            const double p2 = pi[2];
            if (pos_map[i] == true)
            {
                pos.points.push_back(pi);
                pos_proj[i] = ++pos_idx; // 0 means not exist, so all plus 1

                pos_x_min = min(pos_x_min, p0);
                pos_x_max = max(pos_x_max, p0);
                pos_y_min = min(pos_y_min, p1);
                pos_y_max = max(pos_y_max, p1);
                pos_z_min = min(pos_z_min, p2);
                pos_z_max = max(pos_z_max, p2);
            }
            if (neg_map[i] == true)
            {
                neg.points.push_back(pi);
                neg_proj[i] = ++neg_idx;

                neg_x_min = min(neg_x_min, p0);
                neg_x_max = max(neg_x_max, p0);
                neg_y_min = min(neg_y_min, p1);
                neg_y_max = max(neg_y_max, p1);
                neg_z_min = min(neg_z_min, p2);
                neg_z_max = max(neg_z_max, p2);
            }
        }

        const int pos_N = (int)pos.points.size();
        const int neg_N = (int)neg.points.size();

        // border points & triangles
        const int final_count = (int)final_border.size();
        for (int i = 0; i < final_count; i++)
        {
            const vec3d& pi = final_border[i];
            const double p0 = pi[0];
            const double p1 = pi[1];
            const double p2 = pi[2];

            pos.points.push_back(pi);
            neg.points.push_back(pi);

            pos_x_min = min(pos_x_min, p0);
            pos_x_max = max(pos_x_max, p0);
            pos_y_min = min(pos_y_min, p1);
            pos_y_max = max(pos_y_max, p1);
            pos_z_min = min(pos_z_min, p2);
            pos_z_max = max(pos_z_max, p2);

            neg_x_min = min(neg_x_min, p0);
            neg_x_max = max(neg_x_max, p0);
            neg_y_min = min(neg_y_min, p1);
            neg_y_max = max(neg_y_max, p1);
            neg_z_min = min(neg_z_min, p2);
            neg_z_max = max(neg_z_max, p2);
        }

        pos.bbox[0] = pos_x_min;
        pos.bbox[1] = pos_x_max;
        pos.bbox[2] = pos_y_min;
        pos.bbox[3] = pos_y_max;
        pos.bbox[4] = pos_z_min;
        pos.bbox[5] = pos_z_max;

        neg.bbox[0] = neg_x_min;
        neg.bbox[1] = neg_x_max;
        neg.bbox[2] = neg_y_min;
        neg.bbox[3] = neg_y_max;
        neg.bbox[4] = neg_z_min;
        neg.bbox[5] = neg_z_max;

        // triangles
        for (vec3i& triangle : pos.triangles)
        {
            if (triangle[0] >= 0)
            {
                triangle[0] = pos_proj[triangle[0]] - 1;
            }
            else
            {
                triangle[0] = -triangle[0] + pos_N - 1;
            }
            if (triangle[1] >= 0)
            {
                triangle[1] = pos_proj[triangle[1]] - 1;
            }
            else
            {
                triangle[1] = -triangle[1] + pos_N - 1;
            }
            if (triangle[2] >= 0)
            {
                triangle[2] = pos_proj[triangle[2]] - 1;
            }
            else
            {
                triangle[2] = -triangle[2] + pos_N - 1;
            }
        }

        for (vec3i& triangle : neg.triangles)
        {
            if (triangle[0] >= 0)
            {
                triangle[0] = neg_proj[triangle[0]] - 1;
            }
            else
            {
                triangle[0] = -triangle[0] + neg_N - 1;
            }
            if (triangle[1] >= 0)
            {
                triangle[1] = neg_proj[triangle[1]] - 1;
            }
            else
            {
                triangle[1] = -triangle[1] + neg_N - 1;
            }
            if (triangle[2] >= 0)
            {
                triangle[2] = neg_proj[triangle[2]] - 1;
            }
            else
            {
                triangle[2] = -triangle[2] + neg_N - 1;
            }
        }

        const int final_triangle_count = (int)final_triangles.size();
        for (int i = 0; i < final_triangle_count; i++)
        {
            const vec3i& triangle = final_triangles[i];
            cut_area += Area(final_border[triangle[0] - 1], final_border[triangle[1] - 1], final_border[triangle[2] - 1]);
            const double t0 = border_map[triangle[0]] - 1;
            const double t1 = border_map[triangle[1]] - 1;
            const double t2 = border_map[triangle[2]] - 1;
            pos.triangles.emplace_back(pos_N + t0, pos_N + t1, pos_N + t2);
            neg.triangles.emplace_back(neg_N + t2, neg_N + t1, neg_N + t0);
        }

        return true;
    }
}