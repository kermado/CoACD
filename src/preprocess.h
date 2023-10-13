#pragma once
#if PREPROCESS
#include <openvdb/Exceptions.h>
#include <openvdb/openvdb.h>
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/VolumeToMesh.h>
#include <openvdb/util/Util.h>
#endif
#include <vector>
#include <cstdio>
#include <string>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <algorithm>

#include "model_obj.h"
#include "logger.h"
#include "bvh.h"

#if PREPROCESS
using namespace openvdb;
#endif

namespace coacd
{
#ifdef PREPROCESS
    void SDFManifold(Model &input, Model &output, double scale = 50.0f, double level_set = 0.55f);
#endif

    bool IsManifold(Model &input);
}
