#include "vec3.h"
#include "pair.h"

#include "hash/HashSet.h"
#include "hash/HashMap.h"
#include "hash/unordered_dense.h"

#include "CDT.h"

#include <vector>

template <typename KeyT, typename HashT = ankerl::unordered_dense::hash<KeyT>, typename EqT = std::equal_to<KeyT>>
using hash_set = emhash8::HashSet<KeyT, HashT, EqT>;

template <typename KeyT, typename ValueT, typename HashT = ankerl::unordered_dense::hash<KeyT>, typename EqT = std::equal_to<KeyT>>
using hash_map = emhash7::HashMap<KeyT, ValueT, HashT, EqT>;

namespace coacd
{
    struct pair_hash
    {
        inline std::size_t operator()(const pair& p) const
        {
            return hasher(p.i);
        }

    private:
        ankerl::unordered_dense::hash<uint64_t> hasher;
    };

    struct Cache
    {
        std::vector<CDT::V2d<double>> points;

        std::vector<pair> BFS_edges;
        hash_map<pair, pair, pair_hash> tri_edge_map;
        hash_set<pair, pair_hash> tri_border_map;
        hash_set<pair, pair_hash> same_edge_map;
        hash_set<int> overlap_map;
        std::vector<bool> add_vertex;
        std::vector<bool> remove_map;

        std::vector<vec3d> border;
        std::vector<vec3d> overlap;
        std::vector<vec3i> border_triangles;
        std::vector<vec3i> final_triangles;
        std::vector<pair> border_edges;
        hash_map<int, int> border_map;
        std::vector<vec3d> final_border;
        std::vector<bool> pos_map;
        std::vector<bool> neg_map;
        std::vector<int> pos_proj;
        std::vector<int> neg_proj;
        hash_map<pair, int, pair_hash> edge_map;
        hash_map<int, int> vertex_map;

        static Cache& get()
        {
            static thread_local Cache instance;
            return instance;
        }

        void shrink()
        {
            points.clear();
            points.shrink_to_fit();

            tri_edge_map.clear();
            tri_edge_map.shrink_to_fit();

            tri_border_map.clear();
            tri_border_map.shrink_to_fit();

            same_edge_map.clear();
            same_edge_map.shrink_to_fit();

            overlap_map.clear();
            overlap_map.shrink_to_fit();

            add_vertex.clear();
            add_vertex.shrink_to_fit();

            remove_map.clear();
            remove_map.shrink_to_fit();

            border.clear();
            border.shrink_to_fit();

            overlap.clear();
            overlap.shrink_to_fit();

            border_triangles.clear();
            border_triangles.shrink_to_fit();

            final_triangles.clear();
            final_triangles.shrink_to_fit();

            border_edges.clear();
            border_edges.shrink_to_fit();

            border_map.clear();
            border_map.shrink_to_fit();

            final_border.clear();
            final_border.shrink_to_fit();

            pos_map.clear();
            pos_map.shrink_to_fit();

            neg_map.clear();
            neg_map.shrink_to_fit();

            pos_proj.clear();
            pos_proj.shrink_to_fit();

            neg_proj.clear();
            neg_proj.shrink_to_fit();

            edge_map.clear();
            edge_map.shrink_to_fit();

            vertex_map.clear();
            vertex_map.shrink_to_fit();
        }

        void outliers_prepare(const int v_length, const int f_length)
        {
            BFS_edges.clear();
            tri_edge_map.clear();
            tri_border_map.clear();
            same_edge_map.clear();
            overlap_map.clear();

            add_vertex.resize(v_length);
            remove_map.resize(f_length);
        }

        void triangulate_prepare(const unsigned int n)
        {
            points.resize(n);
        }

        void clip_prepare(const unsigned int n)
        {
            border.clear();
            overlap.clear();
            border_triangles.clear();
            final_triangles.clear();
            border_edges.clear();
            border_map.clear();
            final_border.clear();
            edge_map.clear();
            vertex_map.clear();

            pos_map.resize(n);
            neg_map.resize(n);
            pos_proj.resize(n);
            neg_proj.resize(n);
        }
    };
}