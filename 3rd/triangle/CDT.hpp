/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

/**
 * @file
 * Public API - implementation
 */

#include "CDT.h"

#include <algorithm>
#include <deque>
#include <limits>
#include <stdexcept>

namespace CDT
{

CDT_INLINE_IF_HEADER_ONLY VerticesTriangles calculateTrianglesByVertex(const TriangleVec& triangles, const VertInd verticesSize)
{
    VerticesTriangles vertTris(verticesSize);
    for(TriInd iT(0); iT < triangles.size(); ++iT)
    {
        const VerticesArr3& vv = triangles[iT].vertices;
        for(VerticesArr3::const_iterator v = vv.begin(); v != vv.end(); ++v)
        {
            vertTris[*v].push_back(iT);
        }
    }
    return vertTris;
}

template <typename T>
inline DuplicatesInfo RemoveDuplicates(std::vector<V2d<T>>& vertices)
{
    const DuplicatesInfo di = FindDuplicates<T>(vertices.begin(), vertices.end(), getX_V2d<T>, getY_V2d<T>);
    RemoveDuplicates(vertices, di.duplicates);
    return di;
}

CDT_INLINE_IF_HEADER_ONLY void
inline RemapEdges(std::vector<Edge>& edges, const std::vector<std::size_t>& mapping)
{
    RemapEdges(
        edges.begin(),
        edges.end(),
        mapping,
        edge_get_v1,
        edge_get_v2,
        edge_make);
}

template <typename T>
inline DuplicatesInfo RemoveDuplicatesAndRemapEdges(std::vector<V2d<T>>& vertices, std::vector<Edge>& edges)
{
    return RemoveDuplicatesAndRemapEdges<T>(
        vertices,
        getX_V2d<T>,
        getY_V2d<T>,
        edges.begin(),
        edges.end(),
        edge_get_v1,
        edge_get_v2,
        edge_make);
}

CDT_INLINE_IF_HEADER_ONLY EdgeUSet
extractEdgesFromTriangles(const TriangleVec& triangles)
{
    EdgeUSet edges;
    typedef TriangleVec::const_iterator CIt;
    for(CIt t = triangles.begin(); t != triangles.end(); ++t)
    {
        const VertInd v0(t->vertices[0]);
        const VertInd v1(t->vertices[1]);
        const VertInd v2(t->vertices[2]);
        edges.insert(Edge(v0, v1));
        edges.insert(Edge(v1, v2));
        edges.insert(Edge(v2, v0));
    }
    return edges;
}

CDT_INLINE_IF_HEADER_ONLY emhash7::HashMap<Edge, EdgeVec, ankerl::unordered_dense::hash<Edge>>
EdgeToPiecesMapping(const emhash7::HashMap<Edge, EdgeVec, ankerl::unordered_dense::hash<Edge>>& pieceToOriginals)
{
    emhash7::HashMap<Edge, EdgeVec, ankerl::unordered_dense::hash<Edge>> originalToPieces;
    typedef emhash7::HashMap<Edge, EdgeVec, ankerl::unordered_dense::hash<Edge>>::const_iterator Cit;
    for(Cit ptoIt = pieceToOriginals.begin(); ptoIt != pieceToOriginals.end(); ++ptoIt)
    {
        const Edge piece = ptoIt->first;
        const EdgeVec& originals = ptoIt->second;
        for(EdgeVec::const_iterator origIt = originals.begin(); origIt != originals.end(); ++origIt)
        {
            originalToPieces[*origIt].push_back(piece);
        }
    }
    return originalToPieces;
}

} // namespace CDT