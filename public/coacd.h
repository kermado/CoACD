#pragma once
#include "../src/vec3.h"

#include <array>
#include <string>
#include <string_view>
#include <vector>
#include <atomic>

namespace coacd {

#if _WIN32
#define COACD_API __declspec(dllexport)
#else
#define COACD_API
#endif

struct Mesh {
  std::vector<vec3d> vertices;
  std::vector<vec3i> indices;
};

std::vector<Mesh> CoACD(Mesh const &input, std::atomic<bool>& cancel, std::atomic<uint32_t>& progress, double threshold = 0.05,
                        int max_convex_hull = -1, std::string preprocess = "auto",
                        int prep_resolution = 50, int sample_resolution = 2000,
                        int mcts_nodes = 20, int mcts_iteration = 150,
                        int mcts_max_depth = 3, bool pca = false,
                        bool merge = true, unsigned int seed = 0);
void set_log_level(std::string_view level);

} // namespace coacd
