#include "./process.h"
#include "mcts.h"
#include "config.h"
#include "./preprocess.h"

#include <emscripten/threading.h>

#include <iostream>
#include <cmath>
#include <algorithm>
#include <mutex>
#include <thread>

namespace coacd
{
    thread_local std::mt19937 random_engine;

    template<class item_t, class worker_t>
    static void parallel_foreach(item_t* start, item_t* end, worker_t worker)
    {
        const int thread_count = emscripten_num_logical_cores() - 1;
        const int item_count = end - start;
        const int items_per_thread = std::ceil((double)item_count / thread_count);
        
        if (thread_count > 1 && item_count > 1)
        {
            int idx_start = 0;
            
            std::vector<std::thread> threads;
            for (int t = 0; t < thread_count; ++t)
            {
                const int idx_end = std::min(item_count, idx_start + items_per_thread);
                if (idx_end > idx_start)
                {
                    std::thread thread([&worker, start, idx_start, idx_end]()
                    {
                        for (int idx = idx_start; idx < idx_end; ++idx)
                        {
                            worker(*(start + idx), idx);
                        }
                    });

                    threads.push_back(std::move(thread));
                    idx_start = idx_end;
                }
            }

            for (std::thread& thread : threads)
            {
                thread.join();
            }
        }
        else
        {
            for (item_t* current = start; current != end; ++current)
            {
                worker(*current, current - start);
            }
        }
    }

    template<class worker_t>
    static void parallel_for(int start, int end, worker_t worker)
    {
        const int thread_count = emscripten_num_logical_cores() - 1;
        const int item_count = end - start;
        const int items_per_thread = std::ceil((double)item_count / thread_count);

        if (thread_count > 1 && item_count > 1)
        {
            int idx_start = start;
            
            std::vector<std::thread> threads;
            for (int t = 0; t < thread_count; ++t)
            {
                const int idx_end = std::min(end, idx_start + items_per_thread);
                if (idx_end > idx_start)
                {
                    std::thread thread([&worker, idx_start, idx_end]()
                    {
                        for (int idx = idx_start; idx < idx_end; ++idx)
                        {
                            worker(idx);
                        }
                    });

                    threads.push_back(std::move(thread));
                    idx_start = idx_end;
                }
            }

            for (std::thread& thread : threads)
            {
                thread.join();
            }
        }
        else
        {
            for (int idx = start; idx < end; ++idx)
            {
                worker(idx);
            }
        }
    }

    void ManifoldPreprocess(Params &params, Model &m)
    {
#ifdef PREPROCESS
        Model tmp = m;
        m.Clear();
        SDFManifold(tmp, m, params.prep_resolution, params.dmc_thres);
#else
        logger::warn("Preprocessing is not available!");
#endif
    }

    void MergeCH(const Model &ch1, const Model &ch2, Model &ch)
    {
        Model merge;
        merge.points.insert(merge.points.end(), ch1.points.begin(), ch1.points.end());
        merge.points.insert(merge.points.end(), ch2.points.begin(), ch2.points.end());
        merge.triangles.reserve(ch1.triangles.size() + ch2.triangles.size());
        merge.triangles.insert(merge.triangles.end(), ch1.triangles.begin(), ch1.triangles.end());
        const int point_count = (int)ch1.points.size();
        const int triangle_count = (int)ch2.triangles.size();
        for (int i = 0; i < triangle_count; i++)
        {
            const vec3i& triangle = ch2.triangles[i];
            merge.triangles.emplace_back(triangle[0] + point_count, triangle[1] + point_count, triangle[2] + point_count);
        }

        merge.ComputeCH(ch);

        for (int i = 0; i < 3; ++i)
        {
            ch.bbox[i * 2] = min(ch1.bbox[i * 2], ch2.bbox[i * 2]);
            ch.bbox[i * 2 + 1] = max(ch1.bbox[i * 2 + 1], ch2.bbox[i * 2 + 1]);
        }
    }

    static void UpdateBoundingBox(Model& m)
    {
        double bbox[6] { 1E12, -1E12, 1E12, -1E12, 1E12, -1E12 };
        const vector<vec3d>& points = m.points;
        for (int i = 0; i < points.size(); ++i)
        {
            const vec3d& p = points[i];
            for (int j = 0; j < 3; ++j)
            {
                bbox[j * 2] = std::min(p[j], bbox[j * 2]);
                bbox[j * 2 + 1] = std::max(p[j], bbox[j * 2 + 1]);
            }
        }

        std::memcpy(m.bbox, bbox, sizeof(double) * 6);
    }

    static double BoundingBoxDistanceSq(Model& m1, Model& m2)
    {
        double intersection[6];
        for (int i = 0; i < 3; ++i)
        {
            intersection[i * 2] = std::max(m1.bbox[i * 2], m2.bbox[i * 2]);
            intersection[i * 2 + 1] = std::min(m1.bbox[i * 2 + 1], m2.bbox[i * 2 + 1]);
        }

        double distancesq = 0;
        for (int i = 0; i < 3; ++i)
        {
            const double d = intersection[i * 2] - intersection[i * 2 + 1];
            if (d > 0)
            {
                distancesq += d * d;
            }
        }

        return distancesq;
    }

    double MergeConvexHulls(Model &m, vector<Model> &meshs, vector<Model> &cvxs, Params &params, const std::atomic<bool>& cancel, double epsilon, double threshold)
    {
        logger::info(" - Merge Convex Hulls");
        const int nConvexHulls = (int)cvxs.size();
        double h = 0;

        if (nConvexHulls > 1)
        {
            for (int i = 0; i < nConvexHulls; ++i)
            {
                UpdateBoundingBox(cvxs[i]);
            }

            const int bound = ((((nConvexHulls - 1) * nConvexHulls)) >> 1);
            // Populate the cost matrix
            vector<double> costMatrix, precostMatrix;
            costMatrix.resize(bound);    // only keeps the top half of the matrix
            precostMatrix.resize(bound); // only keeps the top half of the matrix

#ifdef PARALLEL
            parallel_for(0, bound, [&costMatrix, &precostMatrix, &meshs, &cvxs, &params, &cancel, threshold](int idx) {
#else
            for (int idx = 0; idx < bound; ++idx) {
#endif
                if (cancel)
                {
#ifdef PARALLEL
                    return;
#else
                    return 0.0;
#endif
                }

                int p1 = (int)(sqrt(8 * idx + 1) - 1) >> 1; // compute nearest triangle number index
                const int sum = (p1 * (p1 + 1)) >> 1;       // compute nearest triangle number from index
                const int p2 = idx - sum;                   // modular arithmetic from triangle number
                p1++;

                Model& m1 = cvxs[p1];
                Model& m2 = cvxs[p2];
                const double bboxdistsq = BoundingBoxDistanceSq(m1, m2);
                if (bboxdistsq < threshold * threshold)
                {
                    const double dist = MeshDist(m1, m2);
                    if (dist < threshold)
                    {
                        Model combinedCH;
                        MergeCH(m1, m2, combinedCH);

                        costMatrix[idx] = ComputeHCost(m1, m2, combinedCH, params.rv_k, params.resolution, params.seed);
                        precostMatrix[idx] = max(ComputeHCost(meshs[p1], m1, params.rv_k, 3000, params.seed),
                                                 ComputeHCost(meshs[p2], m2, params.rv_k, 3000, params.seed));
                    }
                    else
                    {
                        costMatrix[idx] = INF;
                    }
                }
                else
                {
                    costMatrix[idx] = INF;
                }
                
#ifdef PARALLEL
            });
#else
            }
#endif

            int costSize = (int)cvxs.size();

            while (true)
            {
                if (cancel) { return 0.0; }

                // Search for lowest cost
                double bestCost = INF;
                const int addr = FindMinimumElement(costMatrix, &bestCost, 0, (int32_t)costMatrix.size());

                if (params.max_convex_hull <= 0)
                {
                    // if dose not set max nConvexHull, stop the merging when bestCost is larger than the threshold
                    if (bestCost > params.threshold)
                        break;
                    if (bestCost > max(params.threshold - precostMatrix[addr], 0.01)) // avoid merging two parts that have already used up the treshold
                    {
                        costMatrix[addr] = INF;
                        continue;
                    }
                }
                else
                {
                    // if set the max nConvexHull, ignore the threshold limitation and stio the merging untill # part reach the constraint
                    if ((int)cvxs.size() <= params.max_convex_hull && bestCost > params.threshold)
                    {
                        if (bestCost > params.threshold + 0.005 && (int)cvxs.size() == params.max_convex_hull)
                            logger::warn("Max concavity {} exceeds the threshold {} due to {} convex hull limitation", bestCost, params.threshold, params.max_convex_hull);
                        break;
                    }
                    if ((int)cvxs.size() <= params.max_convex_hull && bestCost > max(params.threshold - precostMatrix[addr], 0.01)) // avoid merging two parts that have already used up the treshold
                    {
                        costMatrix[addr] = INF;
                        continue;
                    }
                }

                h = max(h, bestCost);
                const int addrI = (static_cast<int32_t>(sqrt(1 + (8 * addr))) - 1) >> 1;
                const int p1 = addrI + 1;
                const int p2 = addr - ((addrI * (addrI + 1)) >> 1);
                // printf("addr %ld, addrI %ld, p1 %ld, p2 %ld\n", addr, addrI, p1, p2);

                // Make the lowest cost row and column into a new hull
                Model cch;
                MergeCH(cvxs[p1], cvxs[p2], cch);
                cvxs[p2] = cch;

                std::swap(cvxs[p1], cvxs[cvxs.size() - 1]);
                cvxs.pop_back();
                costSize--;

                // Calculate costs versus the new hull
                int rowIdx = ((p2 - 1) * p2) >> 1;

#ifdef PARALLEL
                parallel_for(0, p2, [&cvxs, &costMatrix, &precostMatrix, &params, rowIdx, p2, threshold, bestCost](int i) {
#else
                for (int i = 0; i < p2; ++i) {
#endif
                    Model& m1 = cvxs[p2];
                    Model& m2 = cvxs[i];
                    const double bboxdistsq = BoundingBoxDistanceSq(m1, m2);
                    if (bboxdistsq < threshold * threshold)
                    {
                        const double dist = MeshDist(m1, m2);
                        if (dist < threshold)
                        {
                            Model combinedCH;
                            MergeCH(m1, m2, combinedCH);
                            costMatrix[rowIdx + i] = ComputeHCost(m1, m2, combinedCH, params.rv_k, params.resolution, params.seed);
                            precostMatrix[rowIdx + i] = max(precostMatrix[p2] + bestCost, precostMatrix[i]);
                        }
                        else
                        {
                            costMatrix[rowIdx + i] = INF;
                        }
                    }
                    else
                    {
                        costMatrix[rowIdx + i] = INF;
                    }
#ifdef PARALLEL
                });
#else
                }
#endif

                rowIdx += p2 * 2;

#ifdef PARALLEL
                parallel_for(p2 + 1, costSize, [&cvxs, &costMatrix, &precostMatrix, &params, rowIdx, p2, threshold, bestCost](int i) {
#else
                for (int i = p2 + 1; i < costSize; ++i) {
#endif
                    const int j = i - p2 - 1;
                    const int r = rowIdx + ((p2 + 1) * j) + ((j - 1) * j / 2);

                    Model& m1 = cvxs[p2];
                    Model& m2 = cvxs[i];
                    const double bboxdistsq = BoundingBoxDistanceSq(m1, m2);
                    if (bboxdistsq < threshold * threshold)
                    {
                        const double dist = MeshDist(m1, m2);
                        if (dist < threshold)
                        {
                            Model combinedCH;
                            MergeCH(m1, m2, combinedCH);
                            costMatrix[r] = ComputeHCost(m1, m2, combinedCH, params.rv_k, params.resolution, params.seed);
                            precostMatrix[r] = max(precostMatrix[p2] + bestCost, precostMatrix[i]);
                        }
                        else
                        {
                            costMatrix[r] = INF;
                        }
                    }
                    else
                    {
                        costMatrix[r] = INF;
                    }
#ifdef PARALLEL
                });
#else
                }
#endif

                // Move the top column in to replace its space
                const size_t erase_idx = ((costSize - 1) * costSize) >> 1;
                if (p1 < costSize)
                {
                    rowIdx = (addrI * p1) >> 1;
                    size_t top_row = erase_idx;
                    for (size_t i = 0; i < p1; ++i)
                    {
                        if (i != p2)
                        {
                            costMatrix[rowIdx] = costMatrix[top_row];
                            precostMatrix[rowIdx] = precostMatrix[top_row];
                        }
                        ++rowIdx;
                        ++top_row;
                    }

                    ++top_row;
                    rowIdx += p1;
                    for (size_t i = p1 + 1; i < (costSize + 1); ++i)
                    {
                        costMatrix[rowIdx] = costMatrix[top_row];
                        precostMatrix[rowIdx] = precostMatrix[top_row++];
                        rowIdx += i;
                        assert(rowIdx >= 0);
                    }
                }
                costMatrix.resize(erase_idx);
                precostMatrix.resize(erase_idx);
            }
        }

        return h;
    }
    
    vector<Model> Compute(Model &mesh, Params &params, const std::atomic<bool>& cancel, std::atomic<uint32_t>& progress)
    {
        random_engine.seed(params.seed);

        vector<Model> input_parts { mesh };
        vector<Model> next_input_parts;
        vector<Model> parts, pmeshs;

        logger::info("# Points: {}", mesh.points.size());
        logger::info("# Triangles: {}", mesh.triangles.size());
        logger::info(" - Decomposition (MCTS)");

        clock_t start, end;
        start = clock();

#ifdef PARALLEL
        std::mutex mutex;
#endif

        size_t iter = 0;
        std::atomic<float> concavity(0);
        while ((int)input_parts.size() > 0)
        {
            concavity = 0.0F;

            const int start_sort_idx = (int)parts.size();
            next_input_parts.clear();
            logger::info("iter {} ---- waiting pool: {}", iter, input_parts.size());

#ifdef PARALLEL
            parallel_foreach(input_parts.data(), input_parts.data() + input_parts.size(), [&mutex, &cancel, &concavity, &params, &pmeshs, &parts, &next_input_parts](Model& input_part, int idx) {
                random_engine.seed(params.seed);
#else
            for (int idx = 0; idx < (int)input_parts.size(); ++idx) {
                Model& input_part = input_parts[idx];
#endif
                if (cancel)
                {
#ifdef PARALLEL
                    return;
#else
                    Cache::get().shrink();
                    return vector<Model>();
#endif
                }

                Model ch_part;
                Plane bestplane;
                input_part.ComputeVCH(ch_part);
                const double h = ComputeHCost(input_part, ch_part, params.rv_k, params.resolution, params.seed, 0.0001, false);
                concavity = std::max(concavity.load(), (float)h);

                if (h > params.threshold)
                {
                    vector<Plane> best_path;

                    // MCTS for cutting plane
                    Node *node = new Node(params);
                    State state(params, input_part);
                    node->set_state(state);
                    Node *best_next_node = MonteCarloTreeSearch(params, node, best_path, cancel);
                    
                    if (cancel)
                    {
                        free_tree(node, 0);
                        Cache::get().shrink();

#ifdef PARALLEL
                        return;
#else
                        Cache::get().shrink();
                        return vector<Model>();
#endif
                    }

                    if (best_next_node == NULL)
                    {
#ifdef PARALLEL
                        mutex.lock();
#endif
                        parts.push_back(ch_part);
                        pmeshs.push_back(input_part);
                        parts[parts.size() -1].idx = idx;
                        pmeshs[pmeshs.size() -1].idx = idx;
#ifdef PARALLEL
                        mutex.unlock();
#endif

                        free_tree(node, 0);
                    }
                    else
                    {
                        bestplane = best_next_node->state->current_value.first;
                        TernaryMCTS(input_part, params, bestplane, best_path, best_next_node->quality_value); // using Rv to Ternary refine
                        free_tree(node, 0);

                        Model pos, neg;
                        double cut_area;
                        const bool clipf = Clip(input_part, pos, neg, bestplane, cut_area);

                        if (!clipf)
                        {
                            logger::error("Wrong clip proposal!");
                            Cache::get().shrink();
#ifdef PARALLEL
                            return;
#else
                            Cache::get().shrink();
                            return vector<Model>();
#endif
                        }
#ifdef PARALLEL
                        mutex.lock();
#endif
                        if (pos.triangles.empty() == false)
                        {
                            next_input_parts.push_back(pos);
                            next_input_parts[next_input_parts.size() - 1].idx = idx;
                        }

                        if (neg.triangles.empty() == false)
                        {
                            next_input_parts.push_back(neg);
                            next_input_parts[next_input_parts.size() - 1].idx = idx;
                        }
#ifdef PARALLEL
                        mutex.unlock();
#endif
                    }
                }
                else
                {
#ifdef PARALLEL
                    mutex.lock();
#endif
                    parts.push_back(ch_part);
                    pmeshs.push_back(input_part);
                    parts[parts.size() -1].idx = idx;
                    pmeshs[pmeshs.size() -1].idx = idx;
#ifdef PARALLEL
                    mutex.unlock();
#endif
                }

#ifdef PARALLEL
            });
#else
            }  
#endif

            if (cancel)
            {
#ifndef PARALLEL
                Cache::get().shrink();
#endif
                return vector<Model>();
            }

            progress = *reinterpret_cast<uint32_t*>(&concavity);

            std::stable_sort(parts.begin() + start_sort_idx, parts.end(), [](const Model& p1, const Model& p2) { return p1.idx < p2.idx; });
            std::stable_sort(pmeshs.begin() + start_sort_idx, pmeshs.end(), [](const Model& p1, const Model& p2) { return p1.idx < p2.idx; });
            std::stable_sort(next_input_parts.begin(), next_input_parts.end(), [](const Model& p1, const Model& p2) { return p1.idx < p2.idx; });

            std::swap(next_input_parts, input_parts);
            iter++;
        }

        if (params.merge)
        {
            MergeConvexHulls(mesh, pmeshs, parts, params, cancel);
        }

        if (cancel)
        {
#ifndef PARALLEL
            Cache::get().shrink();
#endif
            return vector<Model>();
        }

        end = clock();
        logger::info("Compute Time: {}s", double(end - start) / CLOCKS_PER_SEC);
        logger::info("# Convex Hulls: {}", (int)parts.size());

        #ifndef PARALLEL
            Cache::get().shrink();
        #endif

        return parts;
    }
}
