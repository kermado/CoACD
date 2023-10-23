#pragma once

#include <iostream>
#include <string>
#include <random>
#include <fstream>
#include <vector>
#include <math.h>
#include <limits>
#include <typeinfo>
#include <algorithm>
#include <assert.h>
#include <regex>

#include "./io.h"
#include "clip.h"
#include "config.h"
#include "model_obj.h"
#include "cost.h"

namespace coacd
{
  extern thread_local std::mt19937 random_engine;

  void ManifoldPreprocess(Params &params, Model &m);
  void MergeCH(Model &ch1, Model &ch2, Model &ch);
  double MergeConvexHulls(Model &m, vector<Model> &meshs, vector<Model> &cvxs, Params &params, const std::atomic<bool>& cancel, double epsilon = 0.02, double threshold = 0.01);
  vector<Model> Compute(Model &mesh, Params &params, const std::atomic<bool>& cancel, std::atomic<uint32_t>& progress);

  inline int32_t FindMinimumElement(const vector<double>& d, double *const m, const int32_t begin, const int32_t end)
  {
    int32_t idx = -1;
    double min = (std::numeric_limits<double>::max)();
    for (size_t i = begin; i < size_t(end); ++i)
    {
      if (d[i] < min)
      {
        idx = i;
        min = d[i];
      }
    }

    *m = min;
    return idx;
  }
}